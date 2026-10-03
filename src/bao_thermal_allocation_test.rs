//! Root-admitted test-only requested-byte controls. Unexecuted; no supported profile.
//! Original requested Layout bytes are distinct from this tracker/header and RSS.
use std::alloc::{GlobalAlloc, Layout, System};
use std::{cell::Cell, ptr};
#[global_allocator]
static ALLOCATOR: RequestedBytes = RequestedBytes;
thread_local! {
    static OWNER_EPOCH: Cell<u64> = const { Cell::new(0) };
    static CONTROL_TOKEN: Cell<u8> = const { Cell::new(0) };
}
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering::SeqCst};

pub(crate) struct RequestedBytes;
#[repr(C)]
#[derive(Clone, Copy)]
struct Header {
    epoch: u64,
    requested: usize,
}
static ACTIVE: AtomicBool = AtomicBool::new(false);
static EPOCH: AtomicU64 = AtomicU64::new(0);
static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
static REALLOC_OVERLAP_PEAK: AtomicUsize = AtomicUsize::new(0);
static CALLS: AtomicUsize = AtomicUsize::new(0);
static MAX_REQUEST: AtomicUsize = AtomicUsize::new(0);
static MAX_ALIGN: AtomicUsize = AtomicUsize::new(0);
static TRACKER_LIVE: AtomicUsize = AtomicUsize::new(0);
static TRACKER_PEAK: AtomicUsize = AtomicUsize::new(0);
static ADOPTED_GROWTH: AtomicUsize = AtomicUsize::new(0);
static INVALID: AtomicBool = AtomicBool::new(false);
static PHASE: AtomicUsize = AtomicUsize::new(0);
static PHASE_PEAK: [AtomicUsize; 4] = [const { AtomicUsize::new(0) }; 4];
static PHASE_OVERLAP: [AtomicUsize; 4] = [const { AtomicUsize::new(0) }; 4];

fn layout(l: Layout) -> Option<(Layout, usize)> {
    let a = l.align().max(std::mem::align_of::<Header>());
    let offset = std::mem::size_of::<Header>().checked_add(a - 1)? & !(a - 1);
    Some((
        Layout::from_size_align(offset.checked_add(l.size().max(1))?, a).ok()?,
        offset,
    ))
}
fn observe(live: usize) {
    PEAK.fetch_max(live, SeqCst);
    PHASE_PEAK[PHASE.load(SeqCst)].fetch_max(live, SeqCst);
}
fn tracking_add(n: usize) {
    let old = TRACKER_LIVE.fetch_add(n, SeqCst);
    if let Some(x) = old.checked_add(n) {
        TRACKER_PEAK.fetch_max(x, SeqCst);
    } else {
        INVALID.store(true, SeqCst);
    }
}
fn tracking_remove(n: usize) {
    if TRACKER_LIVE.fetch_sub(n, SeqCst) < n {
        INVALID.store(true, SeqCst);
    }
}
fn add(n: usize) {
    let old = LIVE.fetch_add(n, SeqCst);
    match old.checked_add(n) {
        Some(live) => observe(live),
        None => {
            INVALID.store(true, SeqCst);
        }
    }
}
fn remove(n: usize) {
    let old = LIVE.fetch_sub(n, SeqCst);
    if old < n {
        INVALID.store(true, SeqCst);
    }
}
unsafe impl GlobalAlloc for RequestedBytes {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let Some((actual, offset)) = layout(l) else {
            return ptr::null_mut();
        };
        let p = unsafe { System.alloc(actual) };
        if p.is_null() {
            return p;
        }
        let epoch = OWNER_EPOCH.try_with(Cell::get).unwrap_or(0);
        unsafe {
            p.cast::<Header>().write(Header {
                epoch,
                requested: l.size(),
            })
        };
        if epoch != 0 {
            CALLS.fetch_add(1, SeqCst);
            MAX_REQUEST.fetch_max(l.size(), SeqCst);
            MAX_ALIGN.fetch_max(l.align(), SeqCst);
            add(l.size());
            tracking_add(actual.size() - l.size());
        }
        unsafe { p.add(offset) }
    }
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        let p = unsafe { self.alloc(l) };
        if !p.is_null() {
            unsafe { ptr::write_bytes(p, 0, l.size()) }
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        let Some((actual, offset)) = layout(l) else {
            INVALID.store(true, SeqCst);
            return;
        };
        let base = unsafe { p.sub(offset) };
        let h = unsafe { base.cast::<Header>().read() };
        if h.requested != l.size() {
            INVALID.store(true, SeqCst);
        }
        if h.epoch != 0 && h.epoch == EPOCH.load(SeqCst) {
            remove(h.requested);
            tracking_remove(actual.size() - l.size());
        }
        unsafe { System.dealloc(base, actual) };
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        let Some((old, offset)) = layout(l) else {
            return ptr::null_mut();
        };
        let Ok(new_l) = Layout::from_size_align(n, l.align()) else {
            return ptr::null_mut();
        };
        let Some((new, new_offset)) = layout(new_l) else {
            return ptr::null_mut();
        };
        if offset != new_offset {
            INVALID.store(true, SeqCst);
            return ptr::null_mut();
        }
        let base = unsafe { p.sub(offset) };
        let h = unsafe { base.cast::<Header>().read() };
        let counted = h.epoch != 0 && h.epoch == EPOCH.load(SeqCst);
        let current = OWNER_EPOCH.try_with(Cell::get).unwrap_or(0);
        let adopted = !counted && current != 0;
        if counted || adopted {
            CALLS.fetch_add(1, SeqCst);
            MAX_REQUEST.fetch_max(n, SeqCst);
            MAX_ALIGN.fetch_max(l.align(), SeqCst);
            // This is a conservative request-overlap envelope, not a claim
            // that System.realloc actually retained both allocations.
            let old_charge = if adopted { h.requested } else { 0 };
            if let Some(overlap) = LIVE
                .load(SeqCst)
                .checked_add(old_charge)
                .and_then(|x| x.checked_add(n))
            {
                REALLOC_OVERLAP_PEAK.fetch_max(overlap, SeqCst);
                PHASE_OVERLAP[PHASE.load(SeqCst)].fetch_max(overlap, SeqCst);
            } else {
                INVALID.store(true, SeqCst);
            }
        }
        let q = unsafe { System.realloc(base, old, new.size()) };
        if q.is_null() {
            return q;
        }
        unsafe {
            q.cast::<Header>().write(Header {
                epoch: if adopted { current } else { h.epoch },
                requested: n,
            })
        };
        if counted {
            remove(h.requested);
            tracking_remove(old.size() - l.size());
            add(n);
            tracking_add(new.size() - n);
        } else if adopted {
            ADOPTED_GROWTH.fetch_add(1, SeqCst);
            add(n);
            tracking_add(new.size() - n);
        }
        unsafe { q.add(offset) }
    }
}

#[derive(Clone, Copy, Debug, serde::Serialize)]
pub(crate) struct Snapshot {
    pub live: usize,
    pub peak: usize,
    pub realloc_overlap_envelope: usize,
    pub calls: usize,
    pub maximum_request: usize,
    pub maximum_alignment: usize,
    pub phase_peak: [usize; 4],
    pub phase_realloc_overlap: [usize; 4],
    pub invalid: bool,
    pub tracker_overhead_live: usize,
    pub tracker_overhead_peak: usize,
    pub adopted_preexisting_reallocations: usize,
}
pub(crate) fn snapshot() -> Snapshot {
    Snapshot {
        live: LIVE.load(SeqCst),
        peak: PEAK.load(SeqCst),
        realloc_overlap_envelope: REALLOC_OVERLAP_PEAK.load(SeqCst),
        calls: CALLS.load(SeqCst),
        maximum_request: MAX_REQUEST.load(SeqCst),
        maximum_alignment: MAX_ALIGN.load(SeqCst),
        phase_peak: std::array::from_fn(|i| PHASE_PEAK[i].load(SeqCst)),
        phase_realloc_overlap: std::array::from_fn(|i| PHASE_OVERLAP[i].load(SeqCst)),
        invalid: INVALID.load(SeqCst),
        tracker_overhead_live: TRACKER_LIVE.load(SeqCst),
        tracker_overhead_peak: TRACKER_PEAK.load(SeqCst),
        adopted_preexisting_reallocations: ADOPTED_GROWTH.load(SeqCst),
    }
}
pub(crate) fn begin() {
    assert!(!ACTIVE.load(SeqCst));
    assert_eq!(
        LIVE.load(SeqCst),
        0,
        "prior epoch still retains a derived owner"
    );
    assert_eq!(TRACKER_LIVE.load(SeqCst), 0);
    assert_ne!(EPOCH.load(SeqCst), u64::MAX);
    let epoch = EPOCH.fetch_add(1, SeqCst) + 1;
    for x in [
        &PEAK,
        &REALLOC_OVERLAP_PEAK,
        &CALLS,
        &MAX_REQUEST,
        &MAX_ALIGN,
        &TRACKER_PEAK,
        &ADOPTED_GROWTH,
    ] {
        x.store(0, SeqCst);
    }
    for x in PHASE_PEAK.iter().chain(PHASE_OVERLAP.iter()) {
        x.store(0, SeqCst);
    }
    INVALID.store(false, SeqCst);
    PHASE.store(0, SeqCst);
    OWNER_EPOCH.with(|x| x.set(epoch));
    ACTIVE.store(true, SeqCst);
}
pub(crate) fn phase(i: usize) {
    assert!(i < 4);
    PHASE.store(i, SeqCst);
    observe(LIVE.load(SeqCst));
}
pub(crate) fn end() -> Snapshot {
    OWNER_EPOCH.with(|x| x.set(0));
    ACTIVE.store(false, SeqCst);
    snapshot()
}

pub(crate) fn native_admitted() -> bool {
    CONTROL_TOKEN.try_with(Cell::get).unwrap_or(0) == 2
}
pub(crate) fn decoder_if_admitted() {
    if CONTROL_TOKEN.try_with(Cell::get).unwrap_or(0) != 0 {
        begin();
    }
}
pub(crate) fn phase_if_active(i: usize) {
    if OWNER_EPOCH.try_with(Cell::get).unwrap_or(0) != 0 {
        phase(i);
    }
}
// Only the explicitly named private unit allocation control can create this
// token. Neither the released binary nor a runtime environment has a bypass.
struct Admission;
impl Admission {
    fn new(native: bool) -> Self {
        CONTROL_TOKEN.with(|x| {
            assert_eq!(x.get(), 0);
            x.set(if native { 2 } else { 1 });
        });
        Self
    }
}
impl Drop for Admission {
    fn drop(&mut self) {
        CONTROL_TOKEN.with(|x| x.set(0));
    }
}

#[derive(Default, Debug)]
struct Census {
    atoms: usize,
    entries: usize,
    maps: usize,
    arrays: usize,
    text: usize,
}
impl Census {
    fn visit(&mut self, v: &serde_json::Value) {
        self.atoms += 1;
        match v {
            serde_json::Value::Object(m) => {
                self.maps += 1;
                self.entries += m.len();
                for (k, v) in m {
                    self.text += k.len();
                    self.visit(v);
                }
            }
            serde_json::Value::Array(a) => {
                self.arrays += 1;
                for v in a {
                    self.visit(v);
                }
            }
            serde_json::Value::String(s) => {
                self.text += s.len();
            }
            _ => {}
        }
    }
    fn within(&self, e: &crate::bao_thermal_run::Envelope) {
        assert!(self.atoms <= e.atoms, "atoms {self:?}");
        assert!(self.entries <= e.entries, "entries {self:?}");
        assert!(self.maps <= e.maps, "maps {self:?}");
        assert!(self.arrays <= e.arrays, "arrays {self:?}");
        assert!(self.text <= e.text, "text {self:?}");
    }
}
fn released(s: Snapshot) {
    assert!(!s.invalid, "invalid tracking {s:?}");
    assert_eq!(s.live, 0, "retained operation owner {s:?}");
    assert_eq!(s.tracker_overhead_live, 0, "retained tracker header {s:?}");
}
fn structural_controls() {
    use serde_json::{Map, Value};
    use std::mem::{align_of, size_of};
    assert_eq!(size_of::<usize>(), 8);
    assert_eq!(size_of::<Option<std::ptr::NonNull<u8>>>(), 8);
    assert!(
        size_of::<String>() <= 32
            && size_of::<Value>() <= 64
            && size_of::<Map<String, Value>>() <= 64
    );
    assert!(
        align_of::<String>() <= 8
            && align_of::<Value>() <= 8
            && align_of::<Map<String, Value>>() <= 8
    );
    println!(
        "THERMAL_RUST_LAYOUT pointer={} string={}/{} value={}/{} map={}/{} request={} wire={} tracker_header={}",
        size_of::<usize>(),
        size_of::<String>(),
        align_of::<String>(),
        size_of::<Value>(),
        align_of::<Value>(),
        size_of::<Map<String, Value>>(),
        align_of::<Map<String, Value>>(),
        size_of::<crate::bao_thermal_run::Request>(),
        crate::bridge::bao_thermal::wire_size(),
        size_of::<Header>()
    );
    for count in [1usize, 11, 12, 143, 144, 1024, 4096] {
        for descending in [false, true] {
            begin();
            let mut m = Map::new();
            for i in 0..count {
                let k = if descending { count - 1 - i } else { i };
                m.insert(format!("{k:04}"), Value::Null);
                // The allowance uses final entries after the pending insert,
                // including all nodes allocated by its transient split.
                assert!(snapshot().peak <= 512 + 1024 * (i + 1));
            }
            let clone = m.clone();
            let peak = snapshot();
            drop(clone);
            drop(m);
            let s = end();
            println!(
                "THERMAL_MAP count={count} descending={descending} peak={peak:?} released={s:?}"
            );
            released(s);
            assert!(peak.maximum_request <= 1280 && peak.maximum_alignment <= 8);
            assert!(peak.peak <= 2 * (512 + 1024 * count));
        }
    }
    // An old owner grown during the operation must be adopted/charged. The
    // old+new overlap is explicitly conservative, including that old capacity.
    let mut preexisting = Vec::<u8>::with_capacity(8);
    preexisting.extend_from_slice(b"12345678");
    begin();
    preexisting.reserve(1000);
    let adopted = snapshot();
    assert_eq!(adopted.adopted_preexisting_reallocations, 1);
    assert!(adopted.live >= 1008 && adopted.realloc_overlap_envelope >= adopted.live + 8);
    drop(preexisting);
    let s = end();
    println!("THERMAL_ADOPTED_GROWTH {adopted:?} released={s:?}");
    released(s);
    // Actual byte/String/Value Vec growth, with old/new request overlap.
    begin();
    let mut bytes = Vec::new();
    let mut text = String::new();
    let mut values = Vec::new();
    for _ in 0..4096 {
        bytes.push(b'x');
        text.push('x');
        values.push(Value::Null);
    }
    let growth = snapshot();
    assert!(growth.realloc_overlap_envelope <= 4 * (4096 * 2 + 4096 * size_of::<Value>()) + 4096);
    drop(bytes);
    drop(text);
    drop(values);
    let s = end();
    println!("THERMAL_VEC_GROWTH {growth:?} released={s:?}");
    released(s);
}
fn inventory(root: &std::path::Path, base: &std::path::Path, entries: &mut Vec<serde_json::Value>) {
    let mut paths: Vec<_> = std::fs::read_dir(root)
        .unwrap()
        .map(|x| x.unwrap().path())
        .collect();
    paths.sort();
    for p in paths {
        if p.is_dir() {
            inventory(&p, base, entries);
        } else {
            let b = std::fs::read(&p).unwrap();
            entries.push(serde_json::json!({"path":p.strip_prefix(base).unwrap().to_str(),"bytes":b.len(),"sha256":crate::records::hash(&b)}));
        }
    }
}
fn actual_case(
    base: &std::path::Path,
    name: &str,
    raw: &[u8],
    expected: Option<&str>,
    e: &crate::bao_thermal_run::Envelope,
    native: bool,
) {
    let path = base.join(name);
    std::fs::create_dir_all(&path).unwrap();
    let request = path.join("request.json");
    let store = path.join("store");
    std::fs::write(&request, raw).unwrap();
    let args = vec![
        "irred-unit-profile".into(),
        "run".into(),
        request.to_str().unwrap().into(),
        store.to_str().unwrap().into(),
    ];
    // Warm the shared stdout context before the operation's allocation epoch.
    println!(
        "THERMAL_PROFILE_CASE name={name} input_sha256={}",
        crate::records::hash(raw)
    );
    let token = Admission::new(native);
    let attempted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        crate::cli::execute_args(&args)
    }));
    let result = match attempted {
        Ok(result) => result,
        Err(payload) => {
            drop(payload);
            Err("PROFILE_RUNNER_PANIC".to_owned())
        }
    };
    if let Err(error) = &result {
        eprintln!("{}", serde_json::json!({"kind":"failure","error_id":error}));
    }
    drop(token);
    let before_return_drop = end();
    // Controller metadata is copied after the epoch; the real returned error
    // owner is then released, rather than falsely declaring it a leak.
    let captured_result = result.clone();
    drop(result);
    let s = snapshot();
    let mut entries = Vec::new();
    inventory(&path, &path, &mut entries);
    // Preserve the real recorder/store, original raw bytes and measurements
    // before any expected result, lifetime, phase or census assertion.
    std::fs::write(path.join("profile.json"),serde_json::to_vec_pretty(&serde_json::json!({"case":name,"result_error":captured_result.as_ref().err(),"expected_error":expected,
        "snapshot":s,"before_return_drop":before_return_drop,"requested_payload_metric":"original Layout bytes; separate conservative realloc envelope","tracker_overhead_excluded":true,
        "native_payload_not_measured_by_rust_tracker":true,"inventory":entries})).unwrap()).unwrap();
    println!("THERMAL_PROFILE_RESULT {name} {s:?} {captured_result:?}");
    assert_eq!(captured_result.as_ref().err().map(String::as_str), expected);
    released(s);
    for (i, bound) in [e.decoder, e.bridge, e.completion, e.recorder]
        .into_iter()
        .enumerate()
    {
        assert!(
            s.phase_peak[i] <= bound && s.phase_realloc_overlap[i] <= bound,
            "phase{i} {s:?}"
        );
    }
    assert!(s.peak <= e.peak && s.realloc_overlap_envelope <= e.peak);
    let attempts = store.join("attempts");
    let complete = std::fs::read_dir(&attempts)
        .unwrap()
        .map(|x| x.unwrap().path())
        .find(|p| p.to_string_lossy().ends_with(".complete.json"))
        .unwrap();
    let record: serde_json::Value =
        serde_json::from_slice(&std::fs::read(complete).unwrap()).unwrap();
    let output: serde_json::Value = serde_json::from_slice(
        &std::fs::read(
            store
                .join("objects")
                .join(record["output_digest"].as_str().unwrap()),
        )
        .unwrap(),
    )
    .unwrap();
    if native {
        assert_eq!(output["thermal_batch_call_attempted"], true);
        assert_eq!(output["source_factor_completed"], true);
    }
    let mut census = Census::default();
    if let Some(spec) = record["scientific_specification_digest"].as_str() {
        let v: serde_json::Value =
            serde_json::from_slice(&std::fs::read(store.join("objects").join(spec)).unwrap())
                .unwrap();
        census.visit(&v);
    }
    census.visit(&output);
    for v in [
        &record["method"],
        &record["precision"],
        &record["resource_budget"]["operation"],
        &record["outputs"],
    ] {
        census.visit(v);
    }
    println!(
        "THERMAL_PROFILE_CENSUS {name} {census:?} dom_bound={} serialized_bound={}",
        e.dom, e.serialized
    );
    census.within(e);
    if expected == Some("NUMERICAL_QUALIFICATION_REQUIRED") && !native {
        assert!(s.phase_peak[3] <= 2 << 20 && s.phase_realloc_overlap[3] <= 2 << 20);
        assert!(
            census.atoms <= 256 && census.entries <= 64 && census.maps <= 16 && census.text <= 4096,
            "short census {census:?}"
        );
        assert_eq!(output["native_payload_absent"], true);
    }
}

#[test]
#[ignore = "root source/profile/runtime lease required; public allocation profile remains false"]
fn allocation_profile_whole_cli() {
    // One explicit test-only token owner, selected serially by its exact name.
    // This witnesses the actual recorder and one unchanged original target;
    // it does not establish independent numerical or physical qualification.
    let manifest: serde_json::Value =
        serde_json::from_str(include_str!(env!("IRRED_BUILD_MANIFEST"))).unwrap();
    assert!(
        manifest["target"]
            .as_str()
            .unwrap()
            .contains("release: 1.98.0")
    );
    assert!(
        manifest["target"]
            .as_str()
            .unwrap()
            .contains("commit-hash: 88d9e12ae")
    );
    structural_controls();
    let raw = include_bytes!("../tests/fixtures/bao-thermal.json");
    assert_eq!(
        crate::records::hash(raw),
        "8730f3bfa2d84b6428b4bea60114eebb50e13b714e8e0fc3ceefea4bc8cb9208"
    );
    let r: Box<crate::bao_thermal_run::Request> = Box::new(serde_json::from_slice(raw).unwrap());
    let chars = r.validate().unwrap();
    let e = crate::bao_thermal_run::envelope(
        raw.len(),
        r.observation.queries.len(),
        r.models.len(),
        r.models
            .iter()
            .map(|x| x.physical_model.species.len())
            .sum(),
        chars,
    )
    .unwrap();
    drop(r);
    let root = std::env::var_os("IRRED_TEST_ARTIFACT_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::path::PathBuf::from("build/test-artifacts/thermal-bao-profile"));
    let base = root.join(format!(
        "{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&base).unwrap();
    actual_case(&base, "original-unchanged", raw, None, &e, true);
    actual_case(
        &base,
        "unsupported-profile-short",
        raw,
        Some("NUMERICAL_QUALIFICATION_REQUIRED"),
        &e,
        false,
    );
    let escaped = std::str::from_utf8(raw)
        .unwrap()
        .replace("synthetic_controls", r"\u0073ynthetic_controls");
    actual_case(
        &base,
        "escaped-decoder-short",
        escaped.as_bytes(),
        Some("NUMERICAL_QUALIFICATION_REQUIRED"),
        &e,
        false,
    );
    // A bounded syntactically valid incomplete object reaches the actual D/C
    // error path, unlike a syntax error refused by the outer strict scanner.
    actual_case(
        &base,
        "decoder-prefix",
        br#"{"operation":"bao.thermal_density","schema_version":2}"#,
        Some("THERMAL_BAO_REQUEST"),
        &e,
        false,
    );
}
