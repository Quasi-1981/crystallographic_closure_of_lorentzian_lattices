#!/usr/bin/env python3
# author: B (lane-B chair a7296aa8), S1765.  Typeset build for the preprint-7 deposit.
#
# WHY THIS FILE EXISTS.  active-v10.2/src/tools/md_to_tex.py carries the inline layer
# (the glyph map, the paragraph assembly, the glyph gate) but handles neither markdown
# TABLES nor fenced CODE BLOCKS -- this paper has ten of the first and many of the
# second.  The preprint-6 build wrapper that did handle them was lost with its
# scratchpad, which is why that deposit now carries hand-mirrored edits.  This wrapper
# is kept in the repository so the same thing does not happen twice.
#
# Source of truth: the Markdown.  The .tex is generated; never edit it by hand.
# Handles: 0.
# RUN LINE:  python S1765_build_tex.py --md <paper.md> --out <paper.tex>
import argparse, os, re, sys
from pathlib import Path

# The inline layer (glyph map, paragraph assembly, glyph gate) lives in md_to_tex.py.
# Look for it BESIDE this file first -- that is how the deposit ships, so the package runs
# from a fresh clone with nothing else present -- and only then fall back to the working
# tree of the programme.  A hard-coded absolute path here is what would make the shipped
# package unrunnable for anyone but its author.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from md_to_tex import _inline, gate_glyph                  # inline layer + glyph gate
except ImportError:
    sys.exit("md_to_tex.py not found: it must sit beside this file (it ships in src/).")

PREAMBLE = r"""%% GENERATED from %(src)s -- DO NOT EDIT BY HAND.
%% Source of truth is the Markdown; regenerate with src/S1765_build_tex.py.
\documentclass[11pt,a4paper,onecolumn]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb,amsthm}
\usepackage[margin=2.4cm]{geometry}
\usepackage{microtype}
\usepackage{longtable,booktabs,array}
\usepackage[hidelinks]{hyperref}
\newcommand{\bookmark}{\textbf{[cited]}}
\setlength{\parskip}{0.55em}
\setlength{\parindent}{0pt}
\setlength{\LTpre}{0.6em}
\setlength{\LTpost}{0.6em}
\title{%(title)s}
\author{%(author)s}
\date{%(date)s}
\begin{document}
\maketitle
"""


# Glyphs this paper uses that the shared inline layer does not map.  Kept HERE rather
# than added to the shared tool: the shared map serves other preprints, and a local gap
# is repaired locally.  The glyph gate at the end of the run is what found every one of
# these -- an unmapped character does not raise, it vanishes silently from the PDF.
MATH = {
    "⟨": r"\langle ", "⟩": r"\rangle ", "ℤ": r"\mathbb{Z}", "ℚ": r"\mathbb{Q}",
    "⊆": r"\subseteq ", "⊇": r"\supseteq ", "⊗": r"\otimes ", "∩": r"\cap ",
    "∉": r"\notin ", "∨": r"{\vee}", "↦": r"\mapsto ", "↪": r"\hookrightarrow ",
    "≇": r"\ncong ", "≢": r"\not\equiv ", "⟸": r"\Longleftarrow ",
    "Γ": r"\Gamma ", "Λ": r"\Lambda ", "Σ": r"\Sigma ", "σ": r"\sigma ",
    "ι": r"\iota ", "′": r"'", "⁶": r"^{6}", "₄": r"_{4}", "₈": r"_{8}",
    "¼": r"\tfrac{1}{4}", "½": r"\tfrac{1}{2}", "¾": r"\tfrac{3}{4}",
    "ȳ": r"\bar{y}", "̄": "",
    "ℓ": r"\ell ", "₅": r"_{5}", "✓": r"\checkmark ",
    "δ": r"\delta ", "₆": r"_{6}",
    "📖": r"\mbox{\bookmark{}}",
    "—": r"\mbox{---}",
}

# Glyphs that are TEXT, not mathematics: wrapping them in $...$ would be wrong even
# outside a math span.  Handled by their own branch in _extra below.
TEXT_ONLY = {"📖": r"\bookmark{}", "—": "---", "̄": ""}


def _extra(s):
    """Apply the local glyph map, MATH-AWARE.

       The shared inline layer already wraps code spans in $...$, so a naive replacement
       inserting its own $ produced `$$\\Lambda$_{W}$` and pdfTeX stopped on "Display math
       should end with $$".  Splitting on $ and treating odd segments as math is what makes
       one map serve both prose and formulae."""
    parts = s.split("$")
    for idx in range(len(parts)):
        inside = (idx % 2 == 1)
        for k, v in MATH.items():
            if k not in parts[idx]:
                continue
            if inside:
                parts[idx] = parts[idx].replace(k, v)
            elif k in TEXT_ONLY:
                parts[idx] = parts[idx].replace(k, TEXT_ONLY[k])
            else:
                parts[idx] = parts[idx].replace(k, "$" + v.strip() + "$")
    return "$".join(parts)


def cell(s):
    """one table cell: unescape the markdown pipe guard, then run the inline layer."""
    return _extra(_inline(s.strip().replace(r"\|", "|")))


def is_rule_row(ln):
    return bool(re.match(r"^\|[\s:\-|]+\|$", ln.strip()))


def emit_table(rows):
    """markdown table -> longtable; first row is the header, second is the alignment rule."""
    header = [c for c in rows[0].strip().strip("|").split("|")]
    ncol = len(header)
    aligns = []
    if len(rows) > 1 and is_rule_row(rows[1]):
        for a in rows[1].strip().strip("|").split("|"):
            a = a.strip()
            aligns.append("r" if a.endswith(":") and not a.startswith(":") else "l")
        body_rows = rows[2:]
    else:
        aligns = ["l"] * ncol
        body_rows = rows[1:]
    while len(aligns) < ncol:
        aligns.append("l")
    # p-columns keep long prose cells inside the text block instead of overflowing
    spec = " ".join("p{%.2f\\linewidth}" % (0.94 / ncol) for _ in range(ncol))
    out = [r"\begin{longtable}{" + spec + "}", r"\toprule"]
    out.append(" & ".join(r"\textbf{%s}" % cell(c) for c in header) + r" \\")
    out.append(r"\midrule\endhead")
    for r in body_rows:
        cells = r.strip().strip("|").split("|")
        cells = (cells + [""] * ncol)[:ncol]
        out.append(" & ".join(cell(c) for c in cells) + r" \\")
    out.append(r"\bottomrule")
    out.append(r"\end{longtable}")
    return out


def emit_verbatim(lines):
    """fenced block -> verbatim, with non-ASCII rendered through the inline layer is NOT
       possible inside verbatim; such lines are emitted as a centred \\texttt paragraph
       instead so that no glyph is silently dropped."""
    if all(ord(ch) < 128 for ln in lines for ch in ln):
        return [r"\begin{verbatim}"] + list(lines) + [r"\end{verbatim}"]
    out = [r"\begin{quote}\ttfamily\small\obeylines"]
    for ln in lines:
        out.append(_extra(_inline(ln)) if ln.strip() else r"~")
    out.append(r"\end{quote}")
    return out


def convert(md_text, src_name, title, author, date):
    md_text = re.sub(r"<!--.*?-->", "", md_text, flags=re.S)
    lines = md_text.split("\n")
    body, i, in_list = [], 0, False

    def close_list():
        nonlocal in_list
        if in_list:
            body.append(r"\end{itemize}")
            in_list = False

    while i < len(lines):
        ln = lines[i].rstrip()
        if not ln.strip():
            close_list()
            body.append("")
            i += 1
            continue
        if ln.startswith("```"):                                   # fenced block
            close_list()
            i += 1
            blk = []
            while i < len(lines) and not lines[i].startswith("```"):
                blk.append(lines[i].rstrip())
                i += 1
            i += 1
            body.extend(emit_verbatim(blk))
            continue
        if ln.lstrip().startswith("|"):                             # table
            close_list()
            rows = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                rows.append(lines[i].rstrip())
                i += 1
            body.extend(emit_table(rows))
            continue
        if ln.startswith("# ") or ln.strip() == "---":
            i += 1
            continue
        if ln.startswith("#### "):
            close_list(); body.append(r"\subsubsection*{" + _extra(_inline(ln[5:])) + "}")
        elif ln.startswith("### "):
            close_list(); body.append(r"\subsection*{" + _extra(_inline(ln[4:])) + "}")
        elif ln.startswith("## "):
            close_list(); body.append(r"\section*{" + _extra(_inline(ln[3:])) + "}")
        elif ln.startswith(">"):
            close_list()
            blk = []
            while i < len(lines) and lines[i].startswith(">"):
                blk.append(lines[i].lstrip(">").strip())
                i += 1
            body.append(r"\begin{quote}" + _extra(_inline(" ".join(blk))) + r"\end{quote}")
            continue
        elif re.match(r"^\d+\.\s", ln):                             # numbered item
            close_list()
            items = []
            while i < len(lines) and (re.match(r"^\d+\.\s", lines[i]) or
                                      (lines[i].startswith("   ") and lines[i].strip())):
                if re.match(r"^\d+\.\s", lines[i]):
                    items.append([re.sub(r"^\d+\.\s", "", lines[i]).rstrip()])
                else:
                    items[-1].append(lines[i].strip())
                i += 1
            body.append(r"\begin{enumerate}")
            for it in items:
                body.append(r"\item " + _extra(_inline(" ".join(it))))
            body.append(r"\end{enumerate}")
            continue
        elif ln.startswith("- "):
            if not in_list:
                body.append(r"\begin{itemize}"); in_list = True
            item = [ln[2:]]
            i += 1
            while i < len(lines) and lines[i].startswith("  ") and lines[i].strip():
                item.append(lines[i].strip()); i += 1
            body.append(r"\item " + _extra(_inline(" ".join(item))))
            continue
        else:
            para = [ln]
            i += 1
            while i < len(lines):
                nxt = lines[i].rstrip()
                if (not nxt.strip() or nxt.startswith(("#", ">", "- ", "```"))
                        or nxt.lstrip().startswith("|") or nxt.strip() == "---"
                        or re.match(r"^\d+\.\s", nxt)):
                    break
                para.append(nxt)
                i += 1
            body.append(_extra(_inline(" ".join(para))))
            continue
        i += 1
    close_list()
    head = PREAMBLE % {"src": src_name, "title": title, "author": author, "date": date}
    return head + "\n".join(body) + "\n\\end{document}\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--author", required=True)
    ap.add_argument("--date", required=True)
    args = ap.parse_args()
    src = Path(args.md)
    tex = convert(src.read_text(encoding="utf-8"), src.name,
                  args.title, args.author, args.date)
    Path(args.out).write_text(tex, encoding="utf-8")
    print("written: %s  (%d lines)" % (args.out, tex.count("\n")))
    print("\n-- GLYPH GATE: any character the inline layer did not map --")
    gate_glyph(tex, src.name)


if __name__ == "__main__":
    main()
