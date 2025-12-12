#!/usr/bin/env python3
"""
callgraph.py — Build a function call graph and render a cleaned, layered SVG.

Pipeline (default):
1) Build edges (AST static or dynamic trace via harness)
2) Write raw DOT with node fill colors by BFS layer from the root
3) Prune noisy nodes (decorators/argmaps/backends; configurable)
4) Transitive reduction via `tred`
5) Layout with `sfdp` to SVG

If Graphviz CLI tools (`tred`, `sfdp`, `dot`) are not available, the script
still writes DOTs and prints exact fallbacks.

Usage (dir mode, local NetworkX checkout):
  python callgraph.py /path/to/networkx \
      networkx.algorithms.similarity.graph_edit_distance \
      --mode trace --harness harness.py \
      --out ged_callgraph.svg \
      --max-depth 6

Usage (single file):
  python callgraph.py /path/to/networkx/networkx/algorithms/similarity.py \
      similarity.graph_edit_distance \
      --mode ast --out ged_ast.svg
"""

from __future__ import annotations
import argparse, ast, importlib, importlib.util, inspect, os, runpy, sys, types, re, shutil, subprocess
from typing import Iterable, Optional, Set, Tuple, Dict, List

try:
    from graphviz import Digraph  # pip install graphviz (CLI graphviz recommended)
except Exception:
    Digraph = None


# ------------------------------ Import helpers --------------------------------

class _TmpPath:
    """Temporarily make a local repo or module file importable."""
    def __init__(self, path: str, module_name_for_file: str = "local_module"):
        self.path = os.path.abspath(path)
        self.mode = "dir" if os.path.isdir(self.path) else "file"
        self.module_name = module_name_for_file
        self.mod = None
        self._old_sys_path = None
    def __enter__(self):
        if self.mode == "dir":
            self._old_sys_path = list(sys.path)
            sys.path.insert(0, self.path)
            return ("dir", None)
        if not (os.path.isfile(self.path) and self.path.endswith(".py")):
            raise FileNotFoundError(f"{self.path} is not a directory or .py file")
        spec = importlib.util.spec_from_file_location(self.module_name, self.path)
        if not spec or not spec.loader:
            raise ImportError(f"Cannot import from {self.path}")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[self.module_name] = mod
        spec.loader.exec_module(mod)  # type: ignore[attr-defined]
        self.mod = mod
        return ("file", mod)
    def __exit__(self, *exc):
        if self.mode == "dir" and self._old_sys_path is not None:
            sys.path[:] = self._old_sys_path


def _infer_prefix(func_path: str) -> Optional[str]:
    parts = [p for p in func_path.split(".") if p]
    return parts[0] if parts else None


def _resolve_from_file_module(root_module: types.ModuleType, dotted: str):
    parts = [p for p in dotted.split(".") if p]
    if parts and parts[0] == root_module.__name__:
        parts = parts[1:]
    obj = root_module
    for p in parts:
        obj = getattr(obj, p)
    return obj


def _resolve_function(code_mode: str, root_obj: Optional[types.ModuleType], func_path: str) -> Tuple[str, object]:
    """Return (qualified_name, function_object) for the requested function."""
    if code_mode == "dir":
        parts = [p for p in func_path.split(".") if p]
        if not parts:
            raise ValueError("Empty function path")
        for cut in range(len(parts), 0, -1):
            mod_name = ".".join(parts[:cut])
            try:
                mod = importlib.import_module(mod_name)
            except Exception:
                continue
            remainder = parts[cut:]
            if not remainder:
                raise ValueError(f"'{func_path}' points to a module, not a function")
            obj = mod
            for p in remainder:
                obj = getattr(obj, p)
            return (f"{mod.__name__}.{'.'.join(remainder)}", obj)
        raise ImportError(f"Cannot import any module from '{func_path}'")
    else:
        if root_obj is None:
            raise RuntimeError("Missing root module")
        obj = _resolve_from_file_module(root_obj, func_path)
        return (f"{root_obj.__name__}.{func_path.split('.')[-1]}", obj)


# ------------------------------ Call extraction --------------------------------

class _CallCollector(ast.NodeVisitor):
    """Collect simple function/method names from Call nodes."""
    def __init__(self): self.calls: Set[str] = set()
    def _qual(self, node: ast.AST) -> Optional[str]:
        parts: List[str] = []
        cur = node
        while isinstance(cur, ast.Attribute):
            parts.append(cur.attr); cur = cur.value  # type: ignore
        if isinstance(cur, ast.Name): parts.append(cur.id)
        parts.reverse(); return ".".join(parts) if parts else None
    def visit_Call(self, node: ast.Call):
        name=None
        if isinstance(node.func, ast.Name): name=node.func.id
        elif isinstance(node.func, ast.Attribute): name=self._qual(node.func)
        if name: self.calls.add(name)
        self.generic_visit(node)


def _static_calls(func_obj: object) -> Set[str]:
    """Return set of callee names from the function's AST (single-file scope)."""
    src = inspect.getsource(func_obj)
    tree = ast.parse(src)
    c = _CallCollector(); c.visit(tree)
    return c.calls


def _dynamic_edges(root_prefix: str, harness_path: str) -> Set[Tuple[str,str]]:
    """Trace calls during harness execution; keep edges wholly under root_prefix."""
    edges:set[Tuple[str,str]] = set()
    def q(frame):
        mod = frame.f_globals.get("__name__","")
        fn = frame.f_code.co_name
        return f"{mod}.{fn}" if mod else fn
    def profiler(frame, event, arg):
        if event == "call":
            caller = frame.f_back
            if caller:
                a = q(caller); b = q(frame)
                if a.startswith(root_prefix) and b.startswith(root_prefix):
                    edges.add((a,b))
        return profiler
    sys.setprofile(profiler)
    try: runpy.run_path(harness_path, run_name="__main__")
    finally: sys.setprofile(None)
    return edges


# ------------------------------ Layering & coloring ----------------------------

def _build_adj(edges: Iterable[Tuple[str,str]]) -> Dict[str, Set[str]]:
    adj: Dict[str, Set[str]] = {}
    for a, b in edges:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set())
    return adj

def _layer_by_bfs(edges: Iterable[Tuple[str,str]], root: str) -> Dict[str, int]:
    """Compute BFS layers (distance from root) on directed edges."""
    adj = _build_adj(edges)
    layer: Dict[str, int] = {}
    if root not in adj:  # root might be the qualified caller
        # Sometimes the exact qualified name used in edges differs slightly.
        # Fall back to the first node that startswith root if present.
        for n in adj:
            if n.endswith(root) or n.startswith(root):
                root = n
                break
        else:
            return layer
    from collections import deque
    dq = deque([(root, 0)])
    layer[root] = 0
    while dq:
        u, d = dq.popleft()
        for v in adj.get(u, ()):
            if v not in layer:
                layer[v] = d + 1
                dq.append((v, d + 1))
    return layer

def _default_layer_palette() -> List[str]:
    # Gentle pastels (blue, green, yellow, red, purple, teal, gray)
    return ["#DBEAFE", "#DCFCE7", "#FDE68A", "#FECACA", "#E9D5FF", "#CFFAFE", "#E5E7EB"]


# ------------------------------ DOT & rendering --------------------------------

def _write_raw_dot(edges: Iterable[Tuple[str,str]], title: str, dot_path: str,
                   node_attrs: Dict[str, Dict[str, str]] | None = None):
    """Write a DOT file, optionally with per-node attributes (e.g., fillcolor)."""
    nodes:set[str] = set()
    for a, b in edges: 
        nodes.add(a); nodes.add(b)

    lines = [
        "digraph G {",
        f'labelloc="t"; label="{title}"; fontsize=20;',
        'node [style="filled", shape="box", penwidth=1];'
    ]

    node_attrs = node_attrs or {}
    for n in sorted(nodes):
        attrs = node_attrs.get(n, {})
        # build attribute string safely
        attr_pairs = []
        for k, v in attrs.items():
            safe_val = str(v).replace('"', '\\"')
            attr_pairs.append(f'{k}="{safe_val}"')
        safe_name = n.replace('"', '\\"')
        if attr_pairs:
            lines.append(f"\"{safe_name}\" [{', '.join(attr_pairs)}];")
        else:
            lines.append(f"\"{safe_name}\";")

    for a, b in sorted(edges):
        sa = a.replace('"','\\"'); sb=b.replace('"','\\"')
        lines.append(f"\"{sa}\" -> \"{sb}\";")

    lines.append("}")
    with open(dot_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))



def _prune_dot(in_path: str, out_path: str, patterns: List[str]):
    """Remove node/edge lines mentioning any regex pattern. Keep braces."""
    regex = re.compile("|".join(patterns)) if patterns else None
    with open(in_path, "r", encoding="utf-8") as f, open(out_path, "w", encoding="utf-8") as g:
        for line in f:
            if line.strip() in {"digraph G {", "}"}:
                g.write(line); continue
            if not regex or not regex.search(line):
                g.write(line)


def _which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

def _run(cmd: List[str]):
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def _render_clean_svg(dot_in: str, svg_out: str) -> bool:
    """Preferred render: tred -> sfdp (SVG). Requires Graphviz CLI tools."""
    tred = _which("tred"); sfdp = _which("sfdp")
    if not (tred and sfdp):
        return False
    tmp_tr = dot_in.replace(".dot", "_tr.dot")
    # tred
    tr = subprocess.run([tred, dot_in], check=True, stdout=subprocess.PIPE)
    with open(tmp_tr, "wb") as f: f.write(tr.stdout)
    # sfdp
    sp = subprocess.run([sfdp, "-x", "-Goverlap=scale", "-Tsvg", tmp_tr], check=True, stdout=subprocess.PIPE)
    with open(svg_out, "wb") as f: f.write(sp.stdout)
    try: os.remove(tmp_tr)
    except OSError: pass
    return True

def _render_fallback(dot_in: str, out_path: str) -> bool:
    """Fallback render via `dot` CLI (or leave DOT if unavailable)."""
    dot_bin = _which("dot")
    if not dot_bin:
        return False
    base, ext = os.path.splitext(out_path)
    fmt = (ext.lstrip(".") or "svg").lower()
    _run([dot_bin, f"-T{fmt}", dot_in, "-o", out_path])
    return True


# ------------------------------------ CLI -------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Build and render a cleaned, layered call-graph SVG.")
    ap.add_argument("code_path", help="Path to code (directory to add to sys.path, or a single .py file).")
    ap.add_argument("function", help="Function path. Dir mode: pkg.mod.func  |  File mode: module.func or func")
    ap.add_argument("--mode", choices=["ast","trace"], default="ast", help="Static AST or dynamic tracing via harness.")
    ap.add_argument("--harness", help="Path to a script that calls the target (required for --mode trace).")
    ap.add_argument("--limit-prefix", default=None, help="Only keep edges where both ends start with this prefix (default inferred).")
    ap.add_argument("--out", default="callgraph_clean.svg", help="Output SVG filename (default: callgraph_clean.svg).")
    ap.add_argument("--no-prune", action="store_true", help="Disable node/edge pruning before layout.")
    ap.add_argument("--prune-pattern", action="append", default=[
        r"networkx\.utils\.decorators",
        r"argmap_",
        r"utils\.backends",
        r"^\"__main__\.",   # helper frames from harness
    ], help="Regex (can repeat) to drop noisy lines from DOT.")
    ap.add_argument("--layer-colors", nargs="+", default=None,
                    help="List of fill colors per BFS layer (wraps if deeper). Defaults to a pastel palette.")
    ap.add_argument("--max-depth", type=int, default=None,
                    help="Optional max BFS layer to keep (drop nodes/edges deeper than this).")
    args = ap.parse_args()

    default_prefix = _infer_prefix(args.function)
    limit_prefix = args.limit_prefix or default_prefix
    palette = args.layer_colors or _default_layer_palette()

    with _TmpPath(args.code_path) as (mode, root_mod):
        qual, func_obj = _resolve_function(mode, root_mod, args.function)

        # Build edges (AST vs trace)
        if args.mode == "ast":
            calls = _static_calls(func_obj)
            caller = qual
            edges: Set[Tuple[str,str]] = set()
            for callee in calls:
                if "." in callee:
                    dotted = callee
                else:
                    mod = getattr(func_obj, "__module__", "")
                    dotted = f"{mod}.{callee}" if mod else callee
                edges.add((caller, dotted))
        else:
            if not args.harness:
                raise SystemExit("--mode trace requires --harness")
            prefix = limit_prefix or getattr(func_obj, "__module__", "").split(".", 1)[0]
            edges = _dynamic_edges(prefix, args.harness)

        # Namespace pruning by prefix (focus the graph)
        if limit_prefix:
            edges = {e for e in edges if e[0].startswith(limit_prefix) and e[1].startswith(limit_prefix)}

        if not edges:
            print("No edges found. Try --mode trace with a harness, or relax --limit-prefix.")
            return

        # Compute BFS layers from the root (the caller node is `qual`)
        root_node = None
        # Prefer exact match; otherwise find a node whose name endswith the qualified target
        node_set = {n for ab in edges for n in ab}
        if qual in node_set:
            root_node = qual
        else:
            for n in node_set:
                if n.endswith(qual):
                    root_node = n
                    break
        if root_node is None:
            # fallback: choose a source-only node as root
            sources = {a for (a, _) in edges} - {b for (_, b) in edges}
            root_node = next(iter(sources), qual)

        layers = _layer_by_bfs(edges, root_node)

        # Optional depth trimming to keep huge graphs readable
        if args.max_depth is not None:
            keep_nodes = {n for n, d in layers.items() if d <= args.max_depth}
            edges = {(a, b) for (a, b) in edges if a in keep_nodes and b in keep_nodes}
            # recompute after trimming
            layers = _layer_by_bfs(edges, root_node)

        # Color nodes by layer
        node_attrs: Dict[str, Dict[str, str]] = {}
        for n in {x for ab in edges for x in ab}:
            d = layers.get(n, 0)
            node_attrs[n] = {"fillcolor": palette[d % len(palette)]}

        # Paths
        base = os.path.splitext(args.out)[0]
        raw_dot = base + "_raw.dot"
        prn_dot = base + ".dot"  # pruned copy we render from

        # 1) write raw DOT (colored)
        _write_raw_dot(edges, f"Call graph: {qual}", raw_dot, node_attrs=node_attrs)

        # 2) prune noisy lines (unless disabled)
        if args.no_prune:
            with open(raw_dot, "r", encoding="utf-8") as f, open(prn_dot, "w", encoding="utf-8") as g:
                g.write(f.read())
        else:
            _prune_dot(raw_dot, prn_dot, args.prune_pattern)

        # 3) render clean SVG via tred+sfdp if available, otherwise fallback
        svg_out = args.out if args.out.lower().endswith(".svg") else base + ".svg"
        rendered = _render_clean_svg(prn_dot, svg_out)
        if not rendered:
            fallback_ok = _render_fallback(prn_dot, svg_out)
            if not fallback_ok:
                print("Rendered SVG not created (missing Graphviz tools).")
                print("Render manually with:")
                print(f"  tred {prn_dot} | sfdp -x -Goverlap=scale -Tsvg > {svg_out}")
                print(f"Or basic layout (no tred) with:")
                print(f"  dot -Tsvg {prn_dot} -o {svg_out}")
                print(f"Raw DOT: {raw_dot}")
                print(f"Pruned DOT: {prn_dot}")
                return

        print(f"SVG written to: {svg_out}")
        print(f"Raw DOT:        {raw_dot}")
        print(f"Pruned DOT:     {prn_dot}")


if __name__ == "__main__":
    main()
