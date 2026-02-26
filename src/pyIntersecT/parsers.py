# parsers.py
# Enhanced parser with element extraction and stoichiometric coefficient handling

import re
from pathlib import Path
import pandas as pd
import logging
from typing import Sequence, Dict, List, Any, Optional, Tuple

logger = logging.getLogger(__name__)
logging.getLogger("pyIntersecT.parsers").setLevel(logging.INFO)


def _find_case_insensitive_key(columns, target: str) -> Optional[str]:
    """Find a column name matching target case-insensitively."""
    target_lower = target.strip().lower()
    for c in columns:
        if isinstance(c, str) and c.strip().lower() == target_lower:
            return c
    return None


def extract_element_info(column_name: str) -> Tuple[str, float]:
    """
    Extract base element name and stoichiometric coefficient from a column name.
    
    Handles patterns like:
    - "Al2O3,mol,pfu" -> ("Al", 2.0)
    - "Na2O" -> ("Na", 2.0)
    - "FeO,mol,pfu" -> ("Fe", 1.0)  # Note: FeO has coef 1, not 2
    - "MgO,mol,pfu" -> ("Mg", 1.0)
    - "Mg[apfu]" -> ("Mg", 1.0)
    - "K2O[wt%]" -> ("K", 2.0)
    - "Fe" -> ("Fe", 1.0)
    - "CaO" -> ("Ca", 1.0)
    
    Returns:
        Tuple of (element_symbol, stoichiometric_coefficient)
    """
    # Clean the column name
    col = column_name.strip()
    
    # Pattern 1: Element with numeric coefficient followed by O (oxide notation)
    # Examples: Al2O3, Na2O, K2O, Fe2O3
    # Note: FeO, MgO, CaO etc. have implicit coefficient of 1
    pattern_oxide = r'^([A-Z][a-z]?)(\d*)O\d*'
    match = re.match(pattern_oxide, col, re.IGNORECASE)
    if match:
        element = match.group(1)
        coef_str = match.group(2)
        # If no number before O, coefficient is 1 (e.g., FeO, MgO, CaO)
        # If there's a number, that's the coefficient (e.g., Al2O3, Na2O, K2O)
        coef = float(coef_str) if coef_str else 1.0
        return (element, coef)
    
    # Pattern 2: Element in brackets or with special notation
    # Examples: Mg[apfu], Fe[wt%], Si(apfu)
    pattern_bracket = r'^([A-Z][a-z]?)\[.*\]|^([A-Z][a-z]?)\(.*\)'
    match = re.match(pattern_bracket, col)
    if match:
        element = match.group(1) or match.group(2)
        return (element, 1.0)
    
    # Pattern 3: Just element symbol (possibly followed by non-alphabetic chars)
    # Examples: Mg, Al, Fe, Si
    pattern_simple = r'^([A-Z][a-z]?)(?:[^a-zA-Z]|$)'
    match = re.match(pattern_simple, col)
    if match:
        element = match.group(1)
        return (element, 1.0)
    
    # Pattern 4: Element symbol at start (fallback)
    # This catches remaining cases where element is clearly at the start
    if len(col) >= 1 and col[0].isupper():
        # Take first capital letter plus optional lowercase
        element = col[0]
        if len(col) > 1 and col[1].islower():
            element += col[1]
        return (element, 1.0)
    
    # No match found - return the original column and coefficient 1
    return (col, 1.0)


def normalize_element_name(element: str) -> str:
    """
    Normalize element names to standard format.
    Handles common variations and returns capitalized element symbol.
    """
    # Common element name variations
    element_map = {
        'si': 'Si', 'al': 'Al', 'fe': 'Fe', 'mg': 'Mg', 'ca': 'Ca',
        'na': 'Na', 'k': 'K', 'ti': 'Ti', 'mn': 'Mn', 'cr': 'Cr',
        'ni': 'Ni', 'p': 'P', 'o': 'O', 'h': 'H', 'c': 'C',
        'cl': 'Cl', 'f': 'F', 's': 'S', 'ba': 'Ba', 'sr': 'Sr',
        'zn': 'Zn', 'cu': 'Cu', 'co': 'Co', 'v': 'V'
    }
    
    elem_lower = element.strip().lower()
    return element_map.get(elem_lower, element.capitalize())


def _auto_detect_coords_from_columns(columns) -> List[str]:
    """Auto-detect coordinate columns from available column names."""
    temp_candidates = ["t(k)", "t", "t[°c]", "temperature", "temp"]
    pres_candidates = ["p(bar)", "p[kbar]", "p", "pressure", "press"]
    x_candidates = ["x", "x[0.0-1.0]", "X(C1)", "x(c1)", "Y(CO2)", "y(co2)"]

    def pick(cands):
        for cand in cands:
            for c in columns:
                if isinstance(c, str) and c.strip().lower() == cand.strip().lower():
                    return c
        return None

    t = pick(temp_candidates)
    p = pick(pres_candidates)
    x = pick(x_candidates)

    # If all three coordinates are present (MAGEMin X-P-T case), return all three
    if t and p and x:
        return [x, p, t]  # Order: composition, pressure, temperature
    
    if t and p:
        return [t, p]
    if x and p:
        return [x, p]
    if x and t:
        return [x, t]

    # fallback: first two columns (excluding obvious non-coords)
    cols = [c for c in columns if isinstance(c, str)]
    exclusions = {"phase", "name", "point", "counter"}
    cols = [c for c in cols if c.strip().lower() not in exclusions]
    if len(cols) >= 2:
        return [cols[0], cols[1]]
    # last resort: return raw first two original columns
    return list(columns)[:2]


def parse_table(path: str | Path,
                coord_columns: Optional[Sequence[str]] = None,
                phase_column: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    Robust parse_table supporting MAGEMin (CSV) and PerpleX (column table).
    coord_columns=None -> automatic detection.
    
    Returns list of calculation blocks, each containing:
    - coords: dict of coordinate values
    - phases: list of phase dicts with composition data
    - system: 'magemin' or 'perplex'
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    # read lines with latin-1 to avoid encoding issues
    text = path.read_text(encoding="latin-1", errors="ignore")
    lines = text.splitlines()

    # find first non-empty line
    first_idx = next((i for i, L in enumerate(lines) if L.strip()), 0)
    first_line = lines[first_idx] if lines else ""

    # Heuristic: if first non-empty line contains many commas -> treat as CSV (MAGEMin)
    is_csv_like = first_line.count(",") >= 2

    df = None
    tried_sep = None

    if is_csv_like:
        # MAGEMin style: try comma-separated
        try:
            try:
                df = pd.read_csv(path, sep=None, engine="python", encoding="utf-8")
            except:
                df = pd.read_csv(path, sep=None, engine="python", encoding="latin1")

        except Exception as e:
            logger.warning("Failed to read as comma CSV: %s; falling back to whitespace", e)
            df = pd.read_csv(path, sep=r"\s+", engine="python", encoding="latin-1", comment="#")
            tried_sep = "whitespace_fallback"
    else:
        # Try to detect PerpleX header line
        header_idx = None
        for i, L in enumerate(lines):
            if "Name" in L and ("T(" in L or "Counter" in L or "P(" in L or "T[K" in L):
                header_idx = i
                break

        if header_idx is not None:
            try:
                df = pd.read_csv(path, sep=r"\s+", engine="python", header=header_idx, 
                               comment="#", encoding="latin-1")
                tried_sep = f"whitespace_with_header@{header_idx}"
            except Exception as e:
                logger.warning("Failed to read PerpleX-style with header at %d: %s", header_idx, e)
                df = pd.read_csv(path, sep=r"\s+", engine="python", encoding="latin-1", comment="#")
                tried_sep = "whitespace_fallback2"
        else:
            # Unknown: try pandas sniff then whitespace fallback
            try:
                df = pd.read_csv(path, sep=None, engine="python", comment="#", encoding="latin-1")
                tried_sep = "sniff"
            except Exception:
                df = pd.read_csv(path, sep=r"\s+", engine="python", comment="#", encoding="latin-1")
                tried_sep = "whitespace_final"

    logger.info("parse_table: used sep=%s, columns=%s", tried_sep, list(df.columns)[:20])

    # Phase column resolution (case-insensitive)
    inferred_phase_col = None
    if phase_column:
        inferred_phase_col = _find_case_insensitive_key(df.columns, phase_column)
    else:
        inferred_phase_col = (_find_case_insensitive_key(df.columns, "phase") or 
                             _find_case_insensitive_key(df.columns, "Name"))

    if inferred_phase_col is None:
        preview_headers = list(df.columns)
        sample_lines = "\n".join(lines[first_idx:first_idx+5])
        raise KeyError(
            f"Phase column not found. Available columns: {preview_headers}\n"
            f"First non-empty lines of file:\n{sample_lines}"
        )

    # Coordinate columns detection or mapping
    if coord_columns is None:
        coord_cols = _auto_detect_coords_from_columns(df.columns)
        logger.info("Auto-detected coordinate columns: %s", coord_cols)
    else:
        coord_cols = []
        for c in coord_columns:
            actual = _find_case_insensitive_key(df.columns, c)
            if actual is None:
                raise KeyError(f"Coordinate column '{c}' not found. Available: {list(df.columns)}")
            coord_cols.append(actual)

    # Group rows into calculation blocks
    blocks: Dict[Tuple[Tuple[str, float], ...], Dict[str, Any]] = {}
    
    for idx, row in df.iterrows():
        coords: Dict[str, float] = {}
        skip_row = False
        
        for col_name in coord_cols:
            try:
                coords[col_name] = float(row[col_name])
            except Exception:
                skip_row = True
                break
        
        if skip_row:
            continue

        key = tuple(sorted((c, coords[c]) for c in coords))
        if key not in blocks:
            system = "magemin" if is_csv_like else "perplex"
            blocks[key] = {
                "coords": coords.copy(), 
                "phases": [], 
                "system": system
            }

        phase_entry = {"phase": row[inferred_phase_col]}
        
        # Process columns and filter for composition data only
        for col in df.columns:
            if col == inferred_phase_col or col in coord_cols:
                continue
                
            # Filter: only process columns that represent chemical compositions
            # These are typically oxide columns ending with O, O2, O3, etc. followed by modifiers
            # Examples: Na2O,mol,pfu, MgO,mol,pfu, Al2O3,mol,pfu, SiO2,mol,pfu, FeO,mol,pfu
            # Skip thermodynamic properties like alpha,1/K, beta,1/bar, cp,J/K/mol, etc.
            
            col_lower = col.lower()
            
            # Explicitly skip known non-composition columns
            skip_patterns = [
                'alpha', 'beta', 'gamma', 'delta',  # Greek letters for properties
                'cp,', 'cv,', 'entropy', 's,j/',     # Heat capacity and entropy
                'gruneisen', 'ks,', 'gs,', 'bulk',   # Elastic properties
                'v,j/', 'h,j/', 'g,j/',              # Thermodynamic potentials
                'vp,', 'vs,', 'v0,', 'vp/', 'rho,',  # Seismic and density
                'counter', 'n,mol', 'n,g',           # Counters and totals
                'mu[', 'wt,%', 'vol,%', 'mol,%',     # Chemical potentials and percentages (not apfu)
                'nom_ox'                              # Nominal oxide (not composition)
            ]
            
            should_skip = any(pattern in col_lower for pattern in skip_patterns)
            if should_skip:
                continue
            
            # Only process columns that look like oxide compositions
            # Must contain 'O' (for oxide) and ideally 'mol' or 'pfu' or 'apfu'
            is_oxide = 'o' in col_lower and ('mol' in col_lower or 'pfu' in col_lower or 'apfu' in col_lower)
            
            # Also allow simple element names followed by [apfu] or similar
            is_element_apfu = '[apfu]' in col_lower or '(apfu)' in col_lower
            
            if not (is_oxide or is_element_apfu):
                continue
            
            val = row[col]
            try:
                num = float(val)
                # Extract element info from column name
                element, coef = extract_element_info(col)
                element_normalized = normalize_element_name(element)
                
                # Store both raw column name and processed element
                phase_entry[col] = num  # Keep original column
                
                # Also store with normalized element name and applied coefficient
                if coef != 1.0:
                    # Apply stoichiometric coefficient for oxides
                    phase_entry[f"{element_normalized}_apfu"] = num * coef
                else:
                    phase_entry[f"{element_normalized}_apfu"] = num
                    
            except Exception:
                # Ignore non-numeric fields
                continue

        blocks[key]["phases"].append(phase_entry)

    logger.info("Parsed %d calculation points with %d total phases", 
                len(blocks), sum(len(b["phases"]) for b in blocks.values()))
    
    return list(blocks.values())