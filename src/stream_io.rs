//! Bounded complete JSONL frames. Never emits truncated JSON to stdout.
use std::io::{self, BufRead, Write};

pub(crate) enum Line {
    Complete(Vec<u8>),
    TooLarge,
    End,
}
pub(crate) fn read_line<R: BufRead>(reader: &mut R, limit: usize) -> io::Result<Line> {
    let mut data = Vec::new();
    let mut large = false;
    let mut seen = false;
    loop {
        let chunk = reader.fill_buf()?;
        if chunk.is_empty() {
            return Ok(if !seen {
                Line::End
            } else if large {
                Line::TooLarge
            } else {
                Line::Complete(data)
            });
        }
        seen = true;
        let end = chunk.iter().position(|byte| *byte == b'\n');
        let count = end.map_or(chunk.len(), |n| n + 1);
        if !large {
            if count > limit.saturating_sub(data.len()) {
                large = true;
                data.clear();
            } else {
                data.extend_from_slice(&chunk[..count]);
            }
        }
        reader.consume(count);
        if end.is_some() {
            return Ok(if large {
                Line::TooLarge
            } else {
                Line::Complete(data)
            });
        }
    }
}
struct Capped {
    data: Vec<u8>,
    limit: usize,
}
impl Write for Capped {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if bytes.len() > self.limit.saturating_sub(self.data.len()) {
            return Err(io::Error::other("OUTPUT_LINE_LIMIT"));
        }
        self.data.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
pub(crate) fn encode<T: serde::Serialize>(value: &T, limit: usize) -> Result<Vec<u8>, String> {
    let payload_limit = limit.checked_sub(1).ok_or("OUTPUT_LINE_LIMIT")?;
    let mut writer = Capped {
        data: Vec::new(),
        limit: payload_limit,
    };
    serde_json::to_writer(&mut writer, value).map_err(|_| "OUTPUT_LINE_LIMIT")?;
    writer.data.push(b'\n');
    Ok(writer.data)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn oversized_input_discards_exactly_one_frame_and_recovers() {
        let mut input = io::BufReader::with_capacity(3, &b"0123456789\n{}\n"[..]);
        assert!(matches!(read_line(&mut input, 4).unwrap(), Line::TooLarge));
        match read_line(&mut input, 4).unwrap() {
            Line::Complete(line) => assert_eq!(line, b"{}\n"),
            _ => panic!("lost next frame"),
        };
        assert!(matches!(read_line(&mut input, 4).unwrap(), Line::End));
    }
    #[test]
    fn reply_is_complete_or_absent_at_exact_boundary() {
        let value = serde_json::json!({"id":1});
        let bytes = encode(&value, 9).unwrap();
        assert_eq!(bytes, b"{\"id\":1}\n");
        assert_eq!(encode(&value, 8), Err("OUTPUT_LINE_LIMIT".into()));
    }
}
