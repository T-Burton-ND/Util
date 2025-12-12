#!/usr/bin/env python3
"""
pp_pickle.py
------------

Recursively inspect a pickle file and print every contained object/value.

Usage
-----
python pp_pickle.py data.pkl
python pp_pickle.py data.pkl --max-list 10
python pp_pickle.py data.pkl --output dump.txt
"""

import argparse
import pickle
import sys
from collections.abc import Mapping
from typing import Any, TextIO


def short_desc(obj: Any) -> str:
    """Compact type/size description for container headers."""
    t = type(obj).__name__
    if isinstance(obj, Mapping):
        return f"{t}(len={len(obj)})"
    if isinstance(obj, (list, tuple, set)):
        return f"{t}(len={len(obj)})"
    if hasattr(obj, "shape"):
        return f"{t}(shape={obj.shape})"
    return t


def format_value(obj: Any, max_len: int = 120) -> str:
    """Readable repr that trims very long outputs."""
    try:
        text = repr(obj)
    except Exception:
        text = f"<unreprable {type(obj).__name__}>"

    if len(text) > max_len:
        return text[: max_len - 3] + "..."
    return text


def is_leaf(obj: Any) -> bool:
    """Whether the object should be treated as a value instead of a container."""
    if isinstance(obj, (str, bytes, bytearray)):
        return True
    if isinstance(obj, (int, float, complex, bool, type(None))):
        return True
    return False


def walk(
    obj: Any,
    path: str = "",
    max_list: int = None,
    indent: int = 0,
    visited=None,
    out: TextIO = sys.stdout,
):
    """Recursively print every object/value inside the pickle."""
    if visited is None:
        visited = set()

    prefix = "  " * indent
    path_label = path or "<root>"

    # Primitive/leaf-like values can be printed directly and reused safely.
    if is_leaf(obj):
        print(f"{prefix}{path_label}: {format_value(obj)}", file=out)
        return

    obj_id = id(obj)
    if obj_id in visited:
        print(f"{prefix}{path_label}: <recursion {short_desc(obj)}>", file=out)
        return
    visited.add(obj_id)

    # Containers
    if isinstance(obj, Mapping):
        print(f"{prefix}{path_label}: {short_desc(obj)}", file=out)
        for k, v in obj.items():
            key_label = repr(k)
            walk(v, f"{path_label}.{key_label}", max_list, indent + 1, visited, out)
        return

    if isinstance(obj, (list, tuple)):
        print(f"{prefix}{path_label}: {short_desc(obj)}", file=out)
        seq = obj if max_list is None else obj[:max_list]
        for i, v in enumerate(seq):
            walk(v, f"{path_label}[{i}]", max_list, indent + 1, visited, out)
        if max_list is not None and len(obj) > max_list:
            print(f"{prefix}  ... ({len(obj) - max_list} more elements)", file=out)
        return

    if isinstance(obj, set):
        print(f"{prefix}{path_label}: {short_desc(obj)}", file=out)
        seq = list(obj) if max_list is None else list(obj)[:max_list]
        for i, v in enumerate(seq):
            walk(v, f"{path_label}{{{i}}}", max_list, indent + 1, visited, out)
        if max_list is not None and len(obj) > max_list:
            print(f"{prefix}  ... ({len(obj) - max_list} more elements)", file=out)
        return

    # Anything else: print type header then inspect __dict__ if present
    print(f"{prefix}{path_label}: {short_desc(obj)} -> {format_value(obj)}", file=out)
    if hasattr(obj, "__dict__"):
        print(f"{prefix}  attributes ({len(obj.__dict__)}):", file=out)
        for k, v in obj.__dict__.items():
            walk(v, f"{path_label}.{k}", max_list, indent + 2, visited, out)


def main():
    parser = argparse.ArgumentParser(
        description="Recursively show the full contents of a pickle file"
    )
    parser.add_argument("pickle_file")
    parser.add_argument(
        "--max-list",
        type=int,
        default=None,
        help="Max list/tuple/set elements to descend into (omit for all)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Write output to a file instead of stdout"
    )
    args = parser.parse_args()

    try:
        with open(args.pickle_file, "rb") as f:
            obj = pickle.load(f)
    except Exception as e:
        print(f"[ERROR] Failed to load pickle: {e}", file=sys.stderr)
        sys.exit(1)

    out_handle = sys.stdout
    try:
        if args.output:
            out_handle = open(args.output, "w", encoding="utf-8")
        walk(obj, max_list=args.max_list, out=out_handle)
    finally:
        if args.output and out_handle is not sys.stdout:
            out_handle.close()


if __name__ == "__main__":
    main()
