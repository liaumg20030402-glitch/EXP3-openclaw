#!/usr/bin/env python
# /// script
# requires-python = ">=3.8"
# dependencies = []   # matplotlib is OPTIONAL: used for PNG charts, else SVG fallback
# ///

"""
stoichiometry_helper.py

A (mostly) dependency-free CLI for reaction analysis. Core balancing and
stoichiometry use only the Python standard library; chart rendering uses
matplotlib if available and otherwise falls back to hand-written SVG.

Subcommands
-----------
balance    Parse an (un)balanced equation, balance with exact rational linear
           algebra, and report integer coefficients + molar masses.
stoich     Balance, then from given reactant amounts compute the limiting
           reagent, theoretical yields and leftovers. Optionally:
             --actual  -> percent yield of a product
             --plot    -> render a yield bar chart (PNG via matplotlib, else SVG)
           Always reports atom economy and a theoretical E-factor.
empirical  Determine empirical (and molecular) formula from combustion data
           (masses of CO2/H2O + sample mass) or from element mass percentages.
thermo     Balance, then compute the standard reaction enthalpy dH_rxn by
           Hess's law from a built-in standard formation-enthalpy table.

Why this is more than a chat LLM can do reliably:
- Exact rational balancing (no hallucinated coefficients on hard redox cases).
- Deterministic, auditable numbers.
- Produces real artifacts (PNG/SVG chart, JSON file) that a text model cannot.

Examples
--------
  python stoichiometry_helper.py balance --equation "C2H6 + O2 -> CO2 + H2O"
  python stoichiometry_helper.py stoich  --equation "C2H6 + O2 -> CO2 + H2O" \\
      --given "C2H6=10g, O2=40g" --actual "CO2=25g" --plot out/ethane.png
  python stoichiometry_helper.py empirical --combustion "CO2=0.881g, H2O=0.180g" \\
      --sample 0.240g --molar-mass 180.16
  python stoichiometry_helper.py thermo --equation "CH4 + O2 -> CO2 + H2O"
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Standard atomic weights (IUPAC conventional values). Z = 1..103.
# ---------------------------------------------------------------------------

ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.91,
    "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91,
    "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
    "Pm": 145.0, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93,
    "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05,
    "Lu": 174.97, "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21,
    "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59,
    "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209.0, "At": 210.0,
    "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0, "Th": 232.04,
    "Pa": 231.04, "U": 238.03, "Np": 237.0, "Pu": 244.0, "Am": 243.0,
    "Cm": 247.0, "Bk": 247.0, "Cf": 251.0, "Es": 252.0, "Fm": 257.0,
    "Md": 258.0, "No": 259.0, "Lr": 262.0,
}

# Standard molar enthalpies of formation dHf° at 298.15 K, in kJ/mol.
# Keyed by formula in the most common standard state (noted in STATE_NOTE).
# Elements in their reference state are 0 by definition. Override per-run with
# --dhf if a different state/value is needed.
FORMATION_ENTHALPY: Dict[str, float] = {
    # elements (reference state)
    "O2": 0.0, "H2": 0.0, "N2": 0.0, "Cl2": 0.0, "Br2": 0.0, "F2": 0.0,
    "C": 0.0, "Fe": 0.0, "Al": 0.0, "Na": 0.0, "K": 0.0, "Ca": 0.0,
    "Mg": 0.0, "Cu": 0.0, "Zn": 0.0, "S": 0.0,
    # oxides of carbon / water
    "CO2": -393.5, "CO": -110.5, "H2O": -285.8, "H2O2": -187.8,
    # hydrocarbons / alcohols
    "CH4": -74.6, "C2H6": -84.0, "C2H4": 52.4, "C2H2": 227.4,
    "C3H8": -103.8, "C4H10": -125.6, "C6H6": 49.1,
    "CH3OH": -238.6, "C2H5OH": -277.6, "C6H12O6": -1273.3,
    # nitrogen / sulfur compounds
    "NH3": -45.9, "NO": 91.3, "NO2": 33.2, "N2O": 81.6, "N2O4": 9.2,
    "HNO3": -174.1, "SO2": -296.8, "SO3": -395.7, "H2SO4": -814.0, "H2S": -20.6,
    # halides / common salts / bases
    "HCl": -92.3, "HBr": -36.3, "HF": -273.3, "NaCl": -411.2, "NaOH": -425.6,
    "KCl": -436.7, "KOH": -424.8, "CaO": -634.9, "Ca(OH)2": -985.2,
    "CaCO3": -1207.6, "MgO": -601.6, "Al2O3": -1675.7, "Fe2O3": -824.2,
    "Fe3O4": -1118.4, "CuO": -157.3, "ZnO": -350.5, "CaCl2": -795.4,
}
STATE_NOTE = ("dHf° assume conventional standard states (H2O as liquid; "
              "metals/oxides/salts solid; gases as written). Use --dhf to override.")

MASS_UNITS = {"kg": 1000.0, "g": 1.0, "mg": 1e-3, "ug": 1e-6, "µg": 1e-6}
MOLE_UNITS = {"mol": 1.0, "mmol": 1e-3, "umol": 1e-6, "µmol": 1e-6}


# ---------------------------------------------------------------------------
# Formula parsing
# ---------------------------------------------------------------------------

_STATE_RE = re.compile(r"\((?:s|l|g|aq|sln)\)", re.IGNORECASE)
_CHARGE_RE = re.compile(r"\^?\d*[+-]+$")
_HYDRATE_SEP = re.compile(r"[·•⋅*.]")


def parse_formula(formula: str) -> Dict[str, int]:
    """Parse a chemical formula into {element: count}. Supports nested ()/[],
    multipliers, and hydrate notation (CuSO4·5H2O). State symbols and a trailing
    charge are ignored. Raises ValueError on malformed input."""
    raw = _STATE_RE.sub("", formula.strip()).strip()
    if not raw:
        raise ValueError("Empty formula.")
    segments = _HYDRATE_SEP.split(raw)
    totals: Dict[str, int] = {}
    for seg in segments:
        seg = _CHARGE_RE.sub("", seg.strip()).strip()
        if not seg:
            continue
        m = re.match(r"^(\d+)(.*)$", seg)
        mult, body = 1, seg
        if m and ((m.group(2) and m.group(2)[0].isalpha()) or m.group(2).startswith("(")):
            mult, body = int(m.group(1)), m.group(2)
        for el, n in _parse_group(body).items():
            totals[el] = totals.get(el, 0) + n * mult
    if not totals:
        raise ValueError(f"No atoms parsed from formula: {formula!r}")
    for el in totals:
        if el not in ATOMIC_WEIGHTS:
            raise ValueError(f"Unknown element symbol: {el!r} in {formula!r}")
    return totals


def _parse_group(s: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    i, n = 0, len(s)
    stack: List[Dict[str, int]] = [counts]
    while i < n:
        ch = s[i]
        if ch in "([{":
            stack.append({}); i += 1
        elif ch in ")]}":
            i += 1
            num, i = _read_int(s, i)
            grp = stack.pop()
            if not stack:
                raise ValueError(f"Unbalanced brackets in {s!r}")
            for el, c in grp.items():
                stack[-1][el] = stack[-1].get(el, 0) + c * num
        elif ch.isalpha():
            m = re.match(r"[A-Z][a-z]?", s[i:])
            if not m:
                raise ValueError(f"Cannot parse element near {s[i:]!r}")
            el = m.group(0); i += len(el)
            num, i = _read_int(s, i)
            stack[-1][el] = stack[-1].get(el, 0) + num
        elif ch.isspace():
            i += 1
        else:
            raise ValueError(f"Unexpected character {ch!r} in {s!r}")
    if len(stack) != 1:
        raise ValueError(f"Unbalanced brackets in {s!r}")
    return counts


def _read_int(s: str, i: int) -> Tuple[int, int]:
    m = re.match(r"\d+", s[i:])
    if m:
        return int(m.group(0)), i + len(m.group(0))
    return 1, i


def molar_mass(formula: str) -> float:
    return sum(ATOMIC_WEIGHTS[el] * n for el, n in parse_formula(formula).items())


# ---------------------------------------------------------------------------
# Equation parsing and balancing
# ---------------------------------------------------------------------------

def split_equation(equation: str) -> Tuple[List[str], List[str]]:
    norm = equation.replace("→", "->").replace("⟶", "->").replace("=>", "->")
    parts = re.split(r"->|=", norm, maxsplit=1)
    if len(parts) != 2:
        raise ValueError("Equation must contain exactly one '->' or '=' "
                         "separating reactants and products.")

    def species(side: str) -> List[str]:
        out = []
        for it in side.split("+"):
            it = it.strip()
            if not it:
                continue
            it = re.sub(r"^\s*\d+(\.\d+)?\s*", "", it)  # drop user coefficient
            out.append(it.strip())
        return out

    reactants, products = species(parts[0]), species(parts[1])
    if not reactants or not products:
        raise ValueError("Both sides of the equation must list at least one species.")
    return reactants, products


def _rref(matrix: List[List[Fraction]]) -> List[List[Fraction]]:
    M = [row[:] for row in matrix]
    rows = len(M)
    cols = len(M[0]) if rows else 0
    pr = 0
    for col in range(cols):
        sel = next((r for r in range(pr, rows) if M[r][col] != 0), None)
        if sel is None:
            continue
        M[pr], M[sel] = M[sel], M[pr]
        pv = M[pr][col]
        M[pr] = [x / pv for x in M[pr]]
        for r in range(rows):
            if r != pr and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f * b for a, b in zip(M[r], M[pr])]
        pr += 1
        if pr == rows:
            break
    return M


def _nullspace_basis(matrix: List[List[Fraction]], ncols: int) -> List[List[Fraction]]:
    if not matrix:
        return [[Fraction(1) if j == k else Fraction(0) for j in range(ncols)]
                for k in range(ncols)]
    R = _rref(matrix)
    pivot_cols = []
    for row in R:
        for j in range(ncols):
            if row[j] == 1 and all(row[k] == 0 for k in range(j)):
                pivot_cols.append(j); break
    pivot_set = set(pivot_cols)
    free_cols = [j for j in range(ncols) if j not in pivot_set]
    basis = []
    for free in free_cols:
        vec = [Fraction(0)] * ncols
        vec[free] = Fraction(1)
        for ri, pc in enumerate(pivot_cols):
            if ri < len(R):
                vec[pc] = -R[ri][free]
        basis.append(vec)
    return basis


def balance_equation(reactants, products) -> Tuple[List[int], List[int]]:
    species = reactants + products
    nr = len(reactants)
    parsed = [parse_formula(s) for s in species]
    elements = sorted({el for p in parsed for el in p})
    matrix = []
    for el in elements:
        matrix.append([Fraction(p.get(el, 0) if i < nr else -p.get(el, 0))
                       for i, p in enumerate(parsed)])
    basis = _nullspace_basis(matrix, len(species))
    if len(basis) == 0:
        raise ValueError("Equation cannot be balanced (no non-trivial solution).")
    if len(basis) > 1:
        raise ValueError("Equation is underdetermined ({} independent reactions). "
                         "Provide a single net reaction.".format(len(basis)))
    vec = basis[0]
    lcm = 1
    for f in vec:
        lcm = lcm * f.denominator // math.gcd(lcm, f.denominator)
    ints = [int(f * lcm) for f in vec]
    if ints and ints[0] < 0:
        ints = [-x for x in ints]
    g = 0
    for x in ints:
        g = math.gcd(g, abs(x))
    if g > 1:
        ints = [x // g for x in ints]
    if any(x <= 0 for x in ints):
        raise ValueError("Balancing produced non-positive coefficients; the "
                         "reaction may be invalid or not a single net reaction.")
    return ints[:nr], ints[nr:]


def format_balanced(reactants, products, rc, pc) -> str:
    def side(sp, co):
        return " + ".join(s if c == 1 else f"{c} {s}" for s, c in zip(sp, co))
    return f"{side(reactants, rc)} -> {side(products, pc)}"


# ---------------------------------------------------------------------------
# Amount parsing
# ---------------------------------------------------------------------------

def parse_amounts(spec: str) -> Dict[str, Tuple[float, str]]:
    out: Dict[str, Tuple[float, str]] = {}
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError(f"Bad amount entry (need species=value unit): {chunk!r}")
        name, rhs = chunk.split("=", 1)
        m = re.match(r"^([0-9]*\.?[0-9]+)\s*([A-Za-zµ%]+)$", rhs.strip())
        if not m:
            raise ValueError(f"Bad amount value (e.g. 10g, 0.5mol): {rhs!r}")
        out[name.strip()] = (float(m.group(1)), m.group(2))
    return out


def parse_value_unit(s: str) -> Tuple[float, str]:
    m = re.match(r"^([0-9]*\.?[0-9]+)\s*([A-Za-zµ]+)$", s.strip())
    if not m:
        raise ValueError(f"Bad value (e.g. 0.240g): {s!r}")
    return float(m.group(1)), m.group(2)


def to_moles(value: float, unit: str, mw: float) -> float:
    if unit in MOLE_UNITS:
        return value * MOLE_UNITS[unit]
    if unit in MASS_UNITS:
        return value * MASS_UNITS[unit] / mw
    raise ValueError(f"Unknown unit {unit!r}. Use: "
                     f"{', '.join(list(MASS_UNITS) + list(MOLE_UNITS))}")


def to_grams(value: float, unit: str, mw: float) -> float:
    if unit in MASS_UNITS:
        return value * MASS_UNITS[unit]
    if unit in MOLE_UNITS:
        return value * MOLE_UNITS[unit] * mw
    raise ValueError(f"Unknown unit {unit!r}.")


# ---------------------------------------------------------------------------
# Chart rendering: matplotlib PNG if available, else hand-written SVG
# ---------------------------------------------------------------------------

def render_bar_chart(title: str, labels: List[str],
                     series: Dict[str, List[float]], out_path: Path,
                     ylabel: str = "mass (g)") -> str:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        names = list(series.keys())
        x = np.arange(len(labels))
        width = 0.8 / max(1, len(names))
        fig, ax = plt.subplots(figsize=(max(4, 1.2 * len(labels) + 1), 4))
        for i, name in enumerate(names):
            ax.bar(x + i * width, series[name], width, label=name)
        ax.set_xticks(x + width * (len(names) - 1) / 2)
        ax.set_xticklabels(labels)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if len(names) > 1:
            ax.legend()
        fig.tight_layout()
        target = out_path if out_path.suffix.lower() == ".png" else out_path.with_suffix(".png")
        fig.savefig(str(target), dpi=140)
        plt.close(fig)
        return str(target.resolve())
    except Exception as e:  # noqa: BLE001 - fall back to SVG
        sys.stderr.write(f"[WARN] matplotlib unavailable ({e}); writing SVG.\n")
        target = out_path.with_suffix(".svg")
        target.write_text(_svg_bar_chart(title, labels, series, ylabel),
                          encoding="utf-8")
        return str(target.resolve())


def _svg_bar_chart(title, labels, series, ylabel) -> str:
    W, H = 640, 400
    ml, mr, mt, mb = 60, 20, 40, 60
    pw, ph = W - ml - mr, H - mt - mb
    names = list(series.keys())
    allv = [v for vals in series.values() for v in vals] or [1.0]
    vmax = max(allv) or 1.0
    palette = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759"]
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
             f'font-family="sans-serif" font-size="12">',
             f'<rect width="{W}" height="{H}" fill="white"/>',
             f'<text x="{W/2}" y="22" text-anchor="middle" font-size="15" '
             f'font-weight="bold">{_esc(title)}</text>',
             f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#333"/>',
             f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#333"/>',
             f'<text x="14" y="{mt+ph/2}" transform="rotate(-90 14 {mt+ph/2})" '
             f'text-anchor="middle">{_esc(ylabel)}</text>']
    n_groups = len(labels)
    gw = pw / max(1, n_groups)
    bw = gw * 0.8 / max(1, len(names))
    for gi, lab in enumerate(labels):
        gx = ml + gi * gw + gw * 0.1
        for si, name in enumerate(names):
            v = series[name][gi]
            bh = (v / vmax) * ph if vmax else 0
            bx = gx + si * bw
            by = mt + ph - bh
            parts.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bw:.1f}" '
                         f'height="{bh:.1f}" fill="{palette[si % len(palette)]}"/>')
            parts.append(f'<text x="{bx+bw/2:.1f}" y="{by-3:.1f}" '
                         f'text-anchor="middle" font-size="10">{v:.2f}</text>')
        parts.append(f'<text x="{ml+gi*gw+gw/2:.1f}" y="{mt+ph+16}" '
                     f'text-anchor="middle">{_esc(lab)}</text>')
    for li, name in enumerate(names):
        lx = ml + li * 110
        parts.append(f'<rect x="{lx}" y="{H-20}" width="12" height="12" '
                     f'fill="{palette[li % len(palette)]}"/>')
        parts.append(f'<text x="{lx+16}" y="{H-10}">{_esc(name)}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def _esc(s: str) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def ascii_bar_chart(labels: List[str], series: Dict[str, List[float]],
                    width: int = 28) -> str:
    """A text bar chart that renders anywhere (no image media needed).

    Robust to chat platforms / models that cannot deliver image attachments.
    """
    allv = [v for vals in series.values() for v in vals] or [1.0]
    vmax = max(allv) or 1.0
    names = list(series.keys())
    lw = max((len(l) for l in labels), default=1)
    nw = max((len(n) for n in names), default=1)
    lines = []
    for gi, lab in enumerate(labels):
        for si, name in enumerate(names):
            v = series[name][gi]
            bars = int(round((v / vmax) * width)) if vmax else 0
            head = lab if si == 0 else ""
            lines.append(f"{head:<{lw}}  {name:<{nw}} |{'#' * bars} {v:.2f}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Empirical / molecular formula helpers
# ---------------------------------------------------------------------------

def integer_ratios(moles: Dict[str, float], max_mult: int = 8,
                   tol: float = 0.08) -> Dict[str, int]:
    """Convert mole amounts to smallest integer subscripts."""
    items = [(el, m) for el, m in moles.items() if m > 1e-9]
    if not items:
        raise ValueError("No positive mole amounts to form a formula.")
    smallest = min(m for _, m in items)
    base = {el: m / smallest for el, m in items}
    for mult in range(1, max_mult + 1):
        scaled = {el: v * mult for el, v in base.items()}
        if all(abs(v - round(v)) <= tol for v in scaled.values()):
            return {el: int(round(v)) for el, v in scaled.items()}
    return {el: int(round(v)) for el, v in base.items()}


def formula_string(subscripts: Dict[str, int]) -> str:
    order = sorted(subscripts, key=lambda e: (e != "C", e != "H", e))
    out = ""
    for el in order:
        n = subscripts[el]
        out += el + (str(n) if n != 1 else "")
    return out


# ---------------------------------------------------------------------------
# Output helper
# ---------------------------------------------------------------------------

def _emit(result: dict, output: Optional[str]) -> None:
    print("[JSON] " + json.dumps(result, ensure_ascii=False))
    if output:
        p = Path(output).expanduser().resolve()
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"[RESULT] json_file={p}")


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_balance(args) -> int:
    reactants, products = split_equation(args.equation)
    rc, pc = balance_equation(reactants, products)
    balanced = format_balanced(reactants, products, rc, pc)
    species = [{"formula": n, "coefficient": c, "molar_mass": round(molar_mass(n), 4)}
               for n, c in list(zip(reactants, rc)) + list(zip(products, pc))]
    print(f"[RESULT] balanced={balanced}")
    _emit({"mode": "balance", "input": args.equation, "balanced": balanced,
           "reactants": reactants, "products": products,
           "reactant_coefficients": rc, "product_coefficients": pc,
           "species": species}, args.output)
    return 0


def cmd_stoich(args) -> int:
    reactants, products = split_equation(args.equation)
    rc, pc = balance_equation(reactants, products)
    balanced = format_balanced(reactants, products, rc, pc)
    print(f"[RESULT] balanced={balanced}")

    coeff, mw = {}, {}
    for n, c in list(zip(reactants, rc)) + list(zip(products, pc)):
        coeff[n], mw[n] = c, molar_mass(n)

    given = parse_amounts(args.given)
    for n in given:
        if n not in coeff:
            raise ValueError(f"Given species {n!r} not in equation. Known: {list(coeff)}")
    if not any(n in reactants for n in given):
        raise ValueError("Provide amounts for at least one reactant.")

    given_moles, extents = {}, {}
    for n, (v, u) in given.items():
        given_moles[n] = to_moles(v, u, mw[n])
        if n in reactants:
            extents[n] = given_moles[n] / coeff[n]
    xi = min(extents.values())
    limiting = [n for n, e in extents.items() if abs(e - xi) < 1e-12]

    reactant_rows = []
    for n in reactants:
        consumed = coeff[n] * xi
        row = {"formula": n, "coefficient": coeff[n], "molar_mass": round(mw[n], 4),
               "is_limiting": n in limiting}
        if n in given_moles:
            sup = given_moles[n]
            row.update({"supplied_mol": round(sup, 6), "supplied_g": round(sup*mw[n], 4),
                        "consumed_mol": round(consumed, 6), "consumed_g": round(consumed*mw[n], 4),
                        "leftover_mol": round(sup-consumed, 6),
                        "leftover_g": round((sup-consumed)*mw[n], 4)})
        else:
            row.update({"required_mol": round(consumed, 6),
                        "required_g": round(consumed*mw[n], 4), "note": "amount not supplied"})
        reactant_rows.append(row)

    target = args.target if args.target in products else products[0]
    product_rows = []
    for n in products:
        produced = coeff[n] * xi
        product_rows.append({"formula": n, "coefficient": coeff[n],
                             "molar_mass": round(mw[n], 4),
                             "theoretical_mol": round(produced, 6),
                             "theoretical_g": round(produced*mw[n], 4)})

    # Green metrics (reaction-intrinsic atom economy; theoretical E-factor).
    reactant_mass_unit = sum(coeff[n]*mw[n] for n in reactants)
    atom_economy = 100.0 * (coeff[target]*mw[target]) / reactant_mass_unit
    consumed_reactant_g = sum(coeff[n]*xi*mw[n] for n in reactants)
    target_g = coeff[target]*xi*mw[target]
    e_factor = (consumed_reactant_g - target_g) / target_g if target_g else float("nan")

    # Optional percent yield from actual obtained product.
    yield_info = None
    if args.actual:
        amap = parse_amounts(args.actual)
        ydata = {}
        for n, (v, u) in amap.items():
            if n not in products:
                raise ValueError(f"--actual species {n!r} is not a product.")
            actual_g = to_grams(v, u, mw[n])
            theo_g = coeff[n]*xi*mw[n]
            pct = 100.0*actual_g/theo_g if theo_g else float("nan")
            ydata[n] = {"actual_g": round(actual_g, 4),
                        "theoretical_g": round(theo_g, 4),
                        "percent_yield": round(pct, 2)}
            print(f"[RESULT] percent_yield[{n}]={round(pct,2)}%")
        yield_info = ydata

    print(f"[RESULT] limiting_reagent={','.join(limiting)}")
    print(f"[RESULT] extent_mol={round(xi,6)}")
    print(f"[RESULT] atom_economy={round(atom_economy,2)}%")

    # Build the yield series once; used by both the text and image charts.
    labels = products
    theo = [coeff[n] * xi * mw[n] for n in products]
    series = {"theoretical (g)": [round(t, 4) for t in theo]}
    if yield_info:
        series["actual (g)"] = [round(yield_info.get(n, {}).get("actual_g", 0.0), 4)
                                for n in products]

    # Always emit a text bar chart: renders in any chat / model, no media needed.
    ascii_chart = ascii_bar_chart(labels, series)
    print("[CHART]")
    for line in ascii_chart.splitlines():
        print(line)

    # Optional PNG/SVG artifact (saved to disk; useful for reports). Sending it
    # as a chat image is best-effort and not required -- the text chart suffices.
    chart_path = None
    if args.plot:
        chart_path = render_bar_chart(
            f"Theoretical yield: {balanced}", labels, series, Path(args.plot),
            ylabel="mass (g)")
        print(f"[RESULT] chart_file={chart_path}")

    _emit({"mode": "stoich", "input": args.equation, "balanced": balanced,
           "extent_of_reaction_mol": round(xi, 6), "limiting_reagent": limiting,
           "target_product": target, "atom_economy_percent": round(atom_economy, 2),
           "e_factor": round(e_factor, 4), "percent_yield": yield_info,
           "reactants": reactant_rows, "products": product_rows,
           "ascii_chart": ascii_chart, "chart_file": chart_path}, args.output)
    return 0


def cmd_empirical(args) -> int:
    moles: Dict[str, float] = {}
    detail = {}
    if args.combustion:
        if not args.sample:
            raise ValueError("--combustion requires --sample (sample mass).")
        amap = parse_amounts(args.combustion)
        sample_g = to_grams(*parse_value_unit(args.sample), mw=1.0)  # mass directly
        mass_el: Dict[str, float] = {}
        # C from CO2, H from H2O.
        for prod, (v, u) in amap.items():
            g = to_grams(v, u, molar_mass(prod))
            comp = parse_formula(prod)
            if prod.upper() in ("CO2",) or comp == {"C": 1, "O": 2}:
                molC = g / molar_mass("CO2")
                mass_el["C"] = mass_el.get("C", 0.0) + molC*ATOMIC_WEIGHTS["C"]
            elif comp == {"H": 2, "O": 1}:
                molH = 2*g/molar_mass("H2O")
                mass_el["H"] = mass_el.get("H", 0.0) + molH*ATOMIC_WEIGHTS["H"]
            else:
                raise ValueError(f"Unsupported combustion product {prod!r}; "
                                 f"use CO2 and H2O.")
        accounted = sum(mass_el.values())
        o_mass = sample_g - accounted
        if o_mass > 1e-4:
            mass_el["O"] = o_mass
        for el, mg in mass_el.items():
            moles[el] = mg / ATOMIC_WEIGHTS[el]
        detail = {"sample_g": round(sample_g, 4),
                  "element_mass_g": {k: round(v, 4) for k, v in mass_el.items()}}
    elif args.percent:
        # Mass percent, e.g. "C=40.0, H=6.7, O=53.3"; assume a 100 g basis.
        pct = {}
        for chunk in args.percent.split(","):
            chunk = chunk.strip()
            if not chunk:
                continue
            if "=" not in chunk:
                raise ValueError(f"Bad --percent entry (need El=value): {chunk!r}")
            el, val = chunk.split("=", 1)
            el = el.strip()
            if el not in ATOMIC_WEIGHTS:
                raise ValueError(f"Unknown element {el!r} in --percent.")
            num = float(val.strip().rstrip("%").strip())
            pct[el] = num
            moles[el] = num / ATOMIC_WEIGHTS[el]  # 100 g basis -> grams = percent
        detail = {"basis": "100 g", "mass_percent": pct}
    else:
        raise ValueError("Provide either --combustion (+--sample) or --percent.")

    subs = integer_ratios(moles)
    empirical = formula_string(subs)
    emp_mass = sum(ATOMIC_WEIGHTS[e]*n for e, n in subs.items())
    result = {"mode": "empirical", "moles": {k: round(v, 6) for k, v in moles.items()},
              "empirical_formula": empirical, "empirical_formula_mass": round(emp_mass, 4),
              **detail}
    print(f"[RESULT] empirical_formula={empirical}")
    print(f"[RESULT] empirical_formula_mass={round(emp_mass,4)}")
    if args.molar_mass:
        mult = round(args.molar_mass / emp_mass)
        mult = max(1, mult)
        molec = formula_string({e: n*mult for e, n in subs.items()})
        result.update({"given_molar_mass": args.molar_mass, "multiplier": mult,
                       "molecular_formula": molec})
        print(f"[RESULT] molecular_formula={molec}")
    _emit(result, args.output)
    return 0


def cmd_thermo(args) -> int:
    reactants, products = split_equation(args.equation)
    rc, pc = balance_equation(reactants, products)
    balanced = format_balanced(reactants, products, rc, pc)
    print(f"[RESULT] balanced={balanced}")

    overrides = {}
    if args.dhf:
        for chunk in args.dhf.split(","):
            chunk = chunk.strip()
            if not chunk:
                continue
            k, val = chunk.split("=")
            overrides[k.strip()] = float(val)

    def dhf(n: str) -> float:
        if n in overrides:
            return overrides[n]
        if n in FORMATION_ENTHALPY:
            return FORMATION_ENTHALPY[n]
        raise ValueError(f"No standard formation enthalpy for {n!r}. "
                         f"Provide it via --dhf \"{n}=<kJ/mol>\".")

    terms = []
    h_prod = 0.0
    for n, c in zip(products, pc):
        val = dhf(n); h_prod += c*val
        terms.append({"species": n, "role": "product", "coefficient": c, "dHf_kJ_mol": val})
    h_react = 0.0
    for n, c in zip(reactants, rc):
        val = dhf(n); h_react += c*val
        terms.append({"species": n, "role": "reactant", "coefficient": c, "dHf_kJ_mol": val})
    dh_rxn = h_prod - h_react
    nature = "放热 (exothermic)" if dh_rxn < 0 else "吸热 (endothermic)"
    print(f"[RESULT] dH_rxn_kJ={round(dh_rxn,2)}")
    print(f"[RESULT] reaction_nature={nature}")
    _emit({"mode": "thermo", "input": args.equation, "balanced": balanced,
           "dH_rxn_kJ_per_mol_rxn": round(dh_rxn, 2), "nature": nature,
           "sum_products_kJ": round(h_prod, 2), "sum_reactants_kJ": round(h_react, 2),
           "terms": terms, "assumptions": STATE_NOTE}, args.output)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="stoichiometry_helper.py",
        description="Balance equations, run stoichiometry & green metrics, "
                    "determine empirical formulas, and compute reaction enthalpy.")
    sub = p.add_subparsers(dest="cmd", required=True)

    pb = sub.add_parser("balance", help="Balance an equation; report molar masses.")
    pb.add_argument("--equation", required=True)
    pb.add_argument("--output", default=None)
    pb.set_defaults(func=cmd_balance)

    ps = sub.add_parser("stoich", help="Stoichiometry + green metrics + optional chart.")
    ps.add_argument("--equation", required=True)
    ps.add_argument("--given", required=True, help="e.g. 'C2H6=10g, O2=40g'")
    ps.add_argument("--actual", default=None, help="Actual product mass, e.g. 'CO2=25g'")
    ps.add_argument("--target", default=None, help="Desired product for atom economy.")
    ps.add_argument("--plot", default=None, help="Path for yield chart (PNG or SVG).")
    ps.add_argument("--output", default=None)
    ps.set_defaults(func=cmd_stoich)

    pe = sub.add_parser("empirical", help="Empirical/molecular formula from data.")
    pe.add_argument("--combustion", default=None, help="e.g. 'CO2=0.881g, H2O=0.180g'")
    pe.add_argument("--sample", default=None, help="Sample mass, e.g. '0.240g'")
    pe.add_argument("--percent", default=None, help="Mass %, e.g. 'C=40.0, H=6.7, O=53.3'")
    pe.add_argument("--molar-mass", type=float, default=None,
                    help="Known molar mass -> molecular formula.")
    pe.add_argument("--output", default=None)
    pe.set_defaults(func=cmd_empirical)

    pt = sub.add_parser("thermo", help="Reaction enthalpy via Hess's law.")
    pt.add_argument("--equation", required=True)
    pt.add_argument("--dhf", default=None,
                    help="Override/add dHf°, e.g. 'C2H6=-84.0, CO2=-393.5'")
    pt.add_argument("--output", default=None)
    pt.set_defaults(func=cmd_thermo)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    try:
        args = build_parser().parse_args()
        raise SystemExit(int(args.func(args)))
    except KeyboardInterrupt:
        print("[ERROR] Interrupted.", file=sys.stderr); raise
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] {type(e).__name__}: {e}", file=sys.stderr)
        raise SystemExit(2)
