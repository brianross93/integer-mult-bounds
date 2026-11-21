"""Bind the production files used for capture to the submitted Git subtree."""

import argparse
import hashlib
import json
from pathlib import Path


HEAD = "9585c1fa2ed5f2316ac09e04b357cb720badfcb8"
PACKAGE = Path("research/triangular-native-supplier")
# GitHub Git Trees API, GET /git/trees/<HEAD>:research/triangular-native-supplier
TREE = "2802c261f08c6d42d0f371a1b2adc822e4372d4d"


def object_id(kind, data):
    return hashlib.sha1(kind + b" " + str(len(data)).encode() + b"\0" + data).digest()


def tree_id(entries):
    rows = []
    for name, value in entries.items():
        directory = isinstance(value, dict)
        mode = b"40000" if directory else b"100644"
        identifier = tree_id(value) if directory else value
        basename = name.encode("utf-8")
        key = basename + (b"/" if directory else b"")
        rows.append((key, mode + b" " + basename + b"\0" + identifier))
    return object_id(b"tree", b"".join(row for _, row in sorted(rows)))


def verify(source):
    package = Path(source) / PACKAGE
    manifest = json.loads(Path(__file__).with_name("source-blobs.json").read_text("utf-8"))
    entries = {}
    for name, expected in manifest.items():
        path = Path(name)
        actual = object_id(b"blob", (package / path).read_bytes())
        if actual.hex() != expected:
            raise ValueError(f"source bytes differ from Git at {name}")
        node = entries
        for component in path.parts[:-1]:
            node = node.setdefault(component, {})
        node[path.name] = actual
    if tree_id(entries).hex() != TREE:
        raise ValueError("source manifest differs from pinned Git subtree")
    print(f"PASS_SOURCE {len(manifest)} Git blobs; tree {TREE}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="root of the pinned source checkout")
    verify(parser.parse_args().source)
