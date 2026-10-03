"""Preserve timestamps of byte-identical generated build inputs."""
from pathlib import Path


def write_if_changed(path, text):
    path = Path(path)
    content = text.encode("utf-8")
    try:
        if path.read_bytes() == content:
            return
    except FileNotFoundError:
        pass
    path.write_bytes(content)
