//! Bounded process-local ownership. IDs never repeat during a context lifetime.
//! Source quota follows shared native ownership, including after parent release.
use std::{cell::Cell, collections::BTreeMap, rc::Rc};

#[derive(Clone, Copy, Debug, serde::Deserialize, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Limits {
    pub sources: usize,
    pub consumers: usize,
    pub retained_bytes: usize,
    pub commands: u64,
    pub input_line_bytes: usize,
    pub output_line_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            sources: 16,
            consumers: 16,
            retained_bytes: 1 << 30,
            commands: 4096,
            input_line_bytes: 16 << 20,
            output_line_bytes: 64 << 20,
        }
    }
}
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum Error {
    Quota,
    UnknownHandle,
    WrongKind,
    Identity,
    IdExhausted,
    OutputLine,
}

struct Charge {
    bytes: usize,
    ledger: Rc<Cell<usize>>,
    source_count: Option<Rc<Cell<usize>>>,
}
// Reserved before source acquisition/native preparation. The dispatcher maps
// this remaining envelope into native preparation preflight, then converts it
// to the owning type's actual retained charge without a quota gap.
pub(crate) struct Peak {
    charge: Charge,
}
impl Peak {
    pub fn bytes(&self) -> usize {
        self.charge.bytes
    }
}
impl Drop for Charge {
    fn drop(&mut self) {
        self.ledger.set(self.ledger.get() - self.bytes);
        if let Some(count) = &self.source_count {
            count.set(count.get() - 1);
        }
    }
}
pub(crate) struct Source<S> {
    pub value: S,
    pub preparation_digest: String,
    _charge: Charge,
}
pub(crate) struct Consumer<S, C> {
    pub value: C,
    pub _source: Option<Rc<Source<S>>>,
    pub preparation_digest: String,
    _charge: Charge,
}
enum Entry<S, C> {
    Source(Rc<Source<S>>),
    Consumer(Consumer<S, C>),
}
pub(crate) struct Context<S, C> {
    limits: Limits,
    entries: BTreeMap<u64, Entry<S, C>>,
    ledger: Rc<Cell<usize>>,
    source_count: Rc<Cell<usize>>,
    next: u64,
}
impl<S, C> Context<S, C> {
    pub fn new(limits: Limits) -> Self {
        Self {
            limits,
            entries: BTreeMap::new(),
            ledger: Rc::new(Cell::new(0)),
            source_count: Rc::new(Cell::new(0)),
            next: 1,
        }
    }
    #[cfg(test)]
    pub fn check_reply(&self, bytes: usize) -> Result<(), Error> {
        if bytes > self.limits.output_line_bytes {
            Err(Error::OutputLine)
        } else {
            Ok(())
        }
    }
    #[cfg(test)]
    pub fn retained_bytes(&self) -> usize {
        self.ledger.get()
    }
    pub fn remaining_bytes(&self) -> usize {
        self.limits.retained_bytes - self.ledger.get()
    }
    pub fn insert_independent_consumer_reserved(
        &mut self,
        value: C,
        digest: String,
        bytes: usize,
        peak: Peak,
    ) -> Result<u64, Error> {
        if !Self::digest(&digest) {
            return Err(Error::Identity);
        }
        if self
            .entries
            .values()
            .filter(|x| matches!(x, Entry::Consumer(_)))
            .count()
            >= self.limits.consumers
        {
            return Err(Error::Quota);
        }
        let charge = self.retain_peak(peak, bytes)?;
        let id = self.id()?;
        self.entries.insert(
            id,
            Entry::Consumer(Consumer {
                value,
                _source: None,
                preparation_digest: digest,
                _charge: charge,
            }),
        );
        Ok(id)
    }
    pub fn reserve_peak(&self, bytes: usize) -> Result<Peak, Error> {
        Ok(Peak {
            charge: self.charge(bytes)?,
        })
    }
    fn retain_peak(&self, mut peak: Peak, bytes: usize) -> Result<Charge, Error> {
        if !Rc::ptr_eq(&peak.charge.ledger, &self.ledger) {
            return Err(Error::WrongKind);
        }
        if bytes > peak.charge.bytes {
            return Err(Error::Quota);
        }
        self.ledger
            .set(self.ledger.get() - (peak.charge.bytes - bytes));
        peak.charge.bytes = bytes;
        Ok(peak.charge)
    }
    pub fn insert_source_reserved(
        &mut self,
        value: S,
        digest: String,
        bytes: usize,
        peak: Peak,
    ) -> Result<u64, Error> {
        if !Self::digest(&digest) {
            return Err(Error::Identity);
        }
        if self.source_count.get() >= self.limits.sources {
            return Err(Error::Quota);
        }
        let mut charge = self.retain_peak(peak, bytes)?;
        let id = self.id()?;
        self.source_count.set(self.source_count.get() + 1);
        charge.source_count = Some(self.source_count.clone());
        self.entries.insert(
            id,
            Entry::Source(Rc::new(Source {
                value,
                preparation_digest: digest,
                _charge: charge,
            })),
        );
        Ok(id)
    }
    pub fn insert_consumer_reserved(
        &mut self,
        source: Rc<Source<S>>,
        value: C,
        digest: String,
        bytes: usize,
        peak: Peak,
    ) -> Result<u64, Error> {
        if !Self::digest(&digest) {
            return Err(Error::Identity);
        }
        if !Rc::ptr_eq(&source._charge.ledger, &self.ledger) {
            return Err(Error::WrongKind);
        }
        if self
            .entries
            .values()
            .filter(|x| matches!(x, Entry::Consumer(_)))
            .count()
            >= self.limits.consumers
        {
            return Err(Error::Quota);
        }
        let charge = self.retain_peak(peak, bytes)?;
        let id = self.id()?;
        self.entries.insert(
            id,
            Entry::Consumer(Consumer {
                value,
                _source: Some(source),
                preparation_digest: digest,
                _charge: charge,
            }),
        );
        Ok(id)
    }
    fn charge(&self, bytes: usize) -> Result<Charge, Error> {
        let total = self.ledger.get().checked_add(bytes).ok_or(Error::Quota)?;
        if total > self.limits.retained_bytes {
            return Err(Error::Quota);
        }
        self.ledger.set(total);
        Ok(Charge {
            bytes,
            ledger: self.ledger.clone(),
            source_count: None,
        })
    }
    fn id(&mut self) -> Result<u64, Error> {
        let id = self.next;
        self.next = self.next.checked_add(1).ok_or(Error::IdExhausted)?;
        Ok(id)
    }
    fn digest(value: &str) -> bool {
        value.len() == 64
            && value
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    }
    #[cfg(test)]
    pub fn insert_source(&mut self, value: S, digest: String, bytes: usize) -> Result<u64, Error> {
        if !Self::digest(&digest) {
            return Err(Error::Identity);
        }
        if self.source_count.get() >= self.limits.sources {
            return Err(Error::Quota);
        }
        let mut charge = self.charge(bytes)?;
        let id = self.id()?;
        self.source_count.set(self.source_count.get() + 1);
        charge.source_count = Some(self.source_count.clone());
        self.entries.insert(
            id,
            Entry::Source(Rc::new(Source {
                value,
                preparation_digest: digest,
                _charge: charge,
            })),
        );
        Ok(id)
    }
    pub fn source(&self, id: u64) -> Result<Rc<Source<S>>, Error> {
        match self.entries.get(&id) {
            Some(Entry::Source(source)) => Ok(source.clone()),
            Some(_) => Err(Error::WrongKind),
            None => Err(Error::UnknownHandle),
        }
    }
    // Source is obtained before native preparation. A failed construction
    // drops its lease without publishing a handle or changing source charge.
    #[cfg(test)]
    pub fn insert_consumer(
        &mut self,
        source: Rc<Source<S>>,
        value: C,
        digest: String,
        bytes: usize,
    ) -> Result<u64, Error> {
        if !Self::digest(&digest) {
            return Err(Error::Identity);
        }
        if !Rc::ptr_eq(&source._charge.ledger, &self.ledger) {
            return Err(Error::WrongKind);
        }
        if self
            .entries
            .values()
            .filter(|x| matches!(x, Entry::Consumer(_)))
            .count()
            >= self.limits.consumers
        {
            return Err(Error::Quota);
        }
        let charge = self.charge(bytes)?;
        let id = self.id()?;
        self.entries.insert(
            id,
            Entry::Consumer(Consumer {
                value,
                _source: Some(source),
                preparation_digest: digest,
                _charge: charge,
            }),
        );
        Ok(id)
    }
    pub fn consumer(&self, id: u64) -> Result<&Consumer<S, C>, Error> {
        match self.entries.get(&id) {
            Some(Entry::Consumer(value)) => Ok(value),
            Some(_) => Err(Error::WrongKind),
            None => Err(Error::UnknownHandle),
        }
    }
    pub fn release(&mut self, id: u64) -> Result<(), Error> {
        self.entries
            .remove(&id)
            .map(drop)
            .ok_or(Error::UnknownHandle)
    }
    // A newly inserted handle is not published until its complete reply has
    // fitted the output envelope. Failed encoding rolls ownership and quota
    // back, while the monotonic ID remains consumed.
    pub fn publish_reply<T>(
        &mut self,
        id: u64,
        encode: impl FnOnce(u64) -> Result<T, Error>,
    ) -> Result<T, Error> {
        if !self.entries.contains_key(&id) {
            return Err(Error::UnknownHandle);
        }
        match encode(id) {
            Ok(reply) => Ok(reply),
            Err(error) => {
                self.entries.remove(&id);
                Err(error)
            }
        }
    }
    // Encode the truthful release acknowledgement before changing ownership.
    // The caller writes the already complete frame; an I/O failure terminates
    // the session, whose remaining owners are then dropped.
    pub fn release_reply<T>(
        &mut self,
        id: u64,
        encode: impl FnOnce(u64) -> Result<T, Error>,
    ) -> Result<T, Error> {
        if !self.entries.contains_key(&id) {
            return Err(Error::UnknownHandle);
        }
        let reply = encode(id)?;
        self.release(id)?;
        Ok(reply)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn digest() -> String {
        "a".repeat(64)
    }
    #[test]
    fn parent_release_keeps_source_and_charge_until_last_consumer() {
        let mut context = Context::new(Limits {
            retained_bytes: 110,
            ..Limits::default()
        });
        let source = context.insert_source(vec![1, 2], digest(), 80).unwrap();
        let a = context
            .insert_consumer(context.source(source).unwrap(), 3, digest(), 10)
            .unwrap();
        let b = context
            .insert_consumer(context.source(source).unwrap(), 4, digest(), 20)
            .unwrap();
        context.release(source).unwrap();
        assert_eq!(context.retained_bytes(), 110);
        assert_eq!(
            context.consumer(a).unwrap()._source.as_ref().unwrap().value,
            [1, 2]
        );
        assert_eq!(context.consumer(a).unwrap().value, 3);
        assert_eq!(context.consumer(b).unwrap().preparation_digest, digest());
        assert_eq!(
            context
                .consumer(a)
                .unwrap()
                ._source
                .as_ref()
                .unwrap()
                .preparation_digest,
            digest()
        );
        assert_eq!(
            context.insert_source(vec![], digest(), 1),
            Err(Error::Quota)
        );
        context.release(a).unwrap();
        assert_eq!(context.retained_bytes(), 100);
        context.release(b).unwrap();
        assert_eq!(context.retained_bytes(), 0);
        assert_eq!(context.release(source), Err(Error::UnknownHandle));
        let next = context.insert_source(vec![], digest(), 110).unwrap();
        assert!(next > b);
    }
    #[test]
    fn failed_prepare_kind_and_stream_limits_do_not_leak_or_reuse() {
        let mut context: Context<(), ()> = Context::new(Limits {
            sources: 1,
            consumers: 1,
            retained_bytes: 9,
            commands: 2,
            input_line_bytes: 3,
            output_line_bytes: 4,
        });
        assert_eq!(
            context.insert_source((), "bad".into(), 1),
            Err(Error::Identity)
        );
        let source = context.insert_source((), digest(), 8).unwrap();
        let lease = context.source(source).unwrap();
        assert_eq!(
            context.insert_consumer(lease.clone(), (), digest(), 2),
            Err(Error::Quota)
        );
        assert_eq!(context.retained_bytes(), 8);
        let consumer = context.insert_consumer(lease, (), digest(), 1).unwrap();
        assert!(matches!(context.source(consumer), Err(Error::WrongKind)));
        assert_eq!(context.check_reply(5), Err(Error::OutputLine));
        assert_eq!(context.check_reply(4), Ok(()));
    }
    #[test]
    fn peak_reservation_precedes_allocation_and_handover_has_no_quota_gap() {
        let mut context: Context<(), ()> = Context::new(Limits {
            retained_bytes: 100,
            ..Limits::default()
        });
        let peak = context.reserve_peak(80).unwrap();
        assert_eq!(peak.bytes(), 80);
        assert!(matches!(context.reserve_peak(21), Err(Error::Quota)));
        let source = context
            .insert_source_reserved((), digest(), 60, peak)
            .unwrap();
        assert_eq!(context.remaining_bytes(), 40);
        {
            let failed = context.reserve_peak(40).unwrap();
            assert_eq!(context.remaining_bytes(), 0);
            drop(failed);
        }
        assert_eq!(context.remaining_bytes(), 40);
        let peak = context.reserve_peak(40).unwrap();
        let consumer = context
            .insert_consumer_reserved(context.source(source).unwrap(), (), digest(), 20, peak)
            .unwrap();
        context.release(source).unwrap();
        assert_eq!(context.remaining_bytes(), 20);
        context.release(consumer).unwrap();
        assert_eq!(context.remaining_bytes(), 100);
    }
    #[test]
    fn unseen_prepare_rolls_back_and_unencoded_release_preserves_owner() {
        let mut context: Context<(), ()> = Context::new(Limits::default());
        let unseen = context.insert_source((), digest(), 80).unwrap();
        assert_eq!(
            context.publish_reply::<Vec<u8>>(unseen, |_| Err(Error::OutputLine)),
            Err(Error::OutputLine)
        );
        assert_eq!(context.retained_bytes(), 0);
        assert!(matches!(context.source(unseen), Err(Error::UnknownHandle)));
        let source = context.insert_source((), digest(), 80).unwrap();
        assert!(source > unseen);
        assert_eq!(
            context.release_reply::<Vec<u8>>(source, |_| Err(Error::OutputLine)),
            Err(Error::OutputLine)
        );
        assert_eq!(context.retained_bytes(), 80);
        context.source(source).unwrap();
        context.release_reply(source, |_| Ok(vec![b'\n'])).unwrap();
        assert_eq!(context.retained_bytes(), 0);
    }
}
