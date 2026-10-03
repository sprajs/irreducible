"""Read a preserved native JSON corpus and check the independent angular route.

This checker never launches the native binary. Importing/running it requires
the root's separate source/profile execution lease.
"""
import argparse
import json
import sys

from supplied_shell_projection_reference import (
    SCHEMA, MAX_SOURCE_BYTES, Payload, Refusal, check_corpus)


def nesting_guard(text):
    depth = 0
    quoted = False
    escaped = False
    for char in text:
        if quoted:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                quoted = False
        elif char == '"':
            quoted = True
        elif char in "[{":
            depth += 1
            if depth > 8:
                raise Refusal("JSON nesting cap")
        elif char in "]}":
            depth -= 1
            if depth < 0:
                raise Refusal("JSON nesting mismatch")
    if quoted or depth:
        raise Refusal("unterminated JSON")


def reject_float(_):
    raise Refusal("JSON floating scalars forbidden: use IEEE64 bits")


def bounded_integer(text):
    if len(text.lstrip("-")) > 8:
        raise Refusal("JSON integer literal cap")
    return int(text)


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise Refusal("duplicate JSON object key")
        result[key] = value
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("corpus", nargs="?", help="preserved JSON file; stdin if omitted")
    args = parser.parse_args()
    try:
        if args.corpus:
            with open(args.corpus, "rb") as stream:
                raw = stream.read(MAX_SOURCE_BYTES + 1)
        else:
            raw = sys.stdin.buffer.read(MAX_SOURCE_BYTES + 1)
        if len(raw) > MAX_SOURCE_BYTES:
            raise Refusal("raw corpus byte cap")
        # Reserve the whole-live parser graph BEFORE JSON allocation.
        Payload(len(raw))
        text = raw.decode("utf-8", errors="strict")
        nesting_guard(text)
        corpus = json.loads(text, parse_float=reject_float,
                            parse_int=bounded_integer,
                            parse_constant=reject_float,
                            object_pairs_hook=unique_object)
        report = check_corpus(corpus, len(raw))
    except (Refusal, UnicodeError, json.JSONDecodeError, OSError, MemoryError,
            ValueError) as error:
        report = {"schema": SCHEMA, "status": "refused", "reason": str(error),
                  "outputs": [], "execution_scope": "saved corpus only; no binary launch"}
    json.dump(report, sys.stdout, separators=(",", ":"), allow_nan=False)
    sys.stdout.write("\n")
    return 0 if report["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
