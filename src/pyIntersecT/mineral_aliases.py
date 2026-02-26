"""
mineral_aliases.py

Robust loader / resolver for mineral name aliases + solvus discrimination.
"""

# TOML: prefer tomllib (Py>=3.11), fallback to tomli for older Pythons
try:
    import tomllib  # Python 3.11+
except Exception:
    import tomli as tomllib  # type: ignore

from pathlib import Path
import logging
import numpy as np
from typing import Dict, Any

logger = logging.getLogger(__name__)
logging.getLogger("pyIntersecT.mineral_aliases").setLevel(logging.INFO)

# Cache for normalized endmember arrays:
_NORMALIZED_ENDMEMBERS_ARRAY: Dict[Any, Any] = {}


class MineralAliases:
    """
    Load mineral aliases from a TOML and resolve phase names case-insensitively.

    toml_path: path to mineral_symbols.toml which should contain sections
        like ["MAGEMin to Warr2021"] and ["PerpleX to Warr2021"] mapping names.
    """

    def __init__(self, toml_path, strict: bool = True, persist_unknowns: bool = False):
        self.strict = strict
        self.persist_unknowns = persist_unknowns

        # store mappings with normalized (lower().strip()) keys
        self.magemin: Dict[str, str] = {}
        self.perplex: Dict[str, str] = {}

        self._load(toml_path)

    def _load(self, toml_path):
        toml_path = Path(toml_path)
        if not toml_path.exists():
            raise FileNotFoundError(f"mineral_symbols TOML not found: {toml_path}")

        with open(toml_path, "rb") as f:
            data = tomllib.load(f)

        for section_name, mapping in data.items():
            section_lower = section_name.lower()
            if "magemin" in section_lower:
                for k, v in mapping.items():
                    if isinstance(k, str):
                        self.magemin[k.strip().lower()] = v
            elif "perple" in section_lower or "perplex" in section_lower:
                for k, v in mapping.items():
                    if isinstance(k, str):
                        self.perplex[k.strip().lower()] = v

        logger.debug("Loaded %d magemin aliases, %d perplex aliases",
                     len(self.magemin), len(self.perplex))

    def resolve(self, phase_name: str, system: str) -> str:
        """
        Resolve phase_name to canonical name using the chosen system ('magemin'|'perplex').
        Case-insensitive. In non-strict mode unknowns are returned as-is (and optionally persisted).
        """
        if not isinstance(phase_name, str):
            raise TypeError("phase_name must be a string")

        key = phase_name.strip().lower()
        system = system.strip().lower()

        if system == "magemin":
            table = self.magemin
        elif system == "perplex" or system == "perple_x" or system == "perplex":
            table = self.perplex
        else:
            raise ValueError(f"Unknown system: {system!r}")

        if key in table:
            return table[key]

        if self.strict:
            raise KeyError(f"Phase '{phase_name}' not found in {system} dictionary")

        # non-strict: return original (or persist a mapping to itself)
        if self.persist_unknowns:
            table[key] = phase_name
        return phase_name


# -------------------------
# Endmember compositions and solvus discrimination
# -------------------------
ENDMEMBER_COMPOSITIONS = {
    "Mca": {
        "Pg": {"Na": 1, "K": 0, "Al": 3, "Si": 3},  # Paragonite (Na-rich)
        "Ms": {"Na": 0, "K": 1, "Al": 3, "Si": 3},  # Muscovite (K-rich)
        "Mrg": {"Ca": 1, "Al": 4, "Si": 2},         # Margarite (Ca-rich)
    },
    "Fsp": {
        "Ab": {"Na": 1, "Al": 1, "Si": 3},      # Albite
        "An": {"Ca": 1, "Al": 2, "Si": 2},      # Anorthite  
        "Or": {"K": 1, "Al": 1, "Si": 3},       # Orthoclase
    },
    "Cpx": {
        "Di": {"Ca": 1, "Mg": 1, "Si": 2},           # Diopside (high-Ca, Mg-rich)
        "Hd": {"Ca": 1, "Fe": 1, "Si": 2},           # Hedenbergite (high-Ca, Fe-rich)
        "Aug": {"Ca": 0.7, "Mg": 0.7, "Fe": 0.3, "Al": 0.3, "Si": 2},  # Augite (intermediate-Ca)
        "Pig": {"Ca": 0.1, "Mg": 0.9, "Fe": 0.1, "Si": 2},  # Pigeonite (low-Ca, Mg-rich)
        "Jd": {"Na": 1, "Al": 1, "Si": 2},           # Jadeite (Na-Al endmember)
        "Aeg": {"Na": 1, "Fe": 1, "Si": 2},          # Aegirine (Na-Fe endmember)
    },
    "Amp": {
        "Tr": {"Ca": 2, "Mg": 5, "Si": 8},           # Tremolite (calcic, Mg-rich)
        "Act": {"Ca": 2, "Mg": 3, "Fe": 2, "Si": 8}, # Actinolite (calcic, intermediate Mg-Fe)
        "Gln": {"Na": 2, "Mg": 3, "Al": 2, "Si": 8}, # Glaucophane (sodic, Al-bearing)
        "Cum": {"Mg": 5, "Fe": 2, "Si": 8},          # Cummingtonite (Ca-poor, Mg-Fe amphibole)
    },
    "Spl": {
        "Chr": {"Fe": 1, "Cr": 2},          # Chromite (Cr-rich spinel)
        "Usp": {"Fe": 2, "Ti": 1},          # Ulvöspinel-like (Ti-bearing spinel)
        "Mag": {"Fe": 3},                   # Magnetite-like (Fe-dominant)
    },
    "Ilm": {
        "Ilm": {"Fe": 1, "Ti": 1},      # Ilmenite (Fe-Ti oxide)
        "Hem": {"Fe": 2},               # Hematite (Fe2O3)
    },
    "Carb": {
        "Cal": {"Ca": 1},               # Calcite (Ca-only carbonate)
        "Dol": {"Ca": 1, "Mg": 1},      # Dolomite (Ca-Mg carbonate)
        "Ank": {"Ca": 1, "Fe": 1},      # Ankerite (Ca-Fe carbonate)
        "Sid": {"Fe": 1},               # Siderite (Fe carbonate)
    },
}

ENDMEMBER_DISTANCE_THRESHOLD = {
    k: 0.3 for k in ENDMEMBER_COMPOSITIONS.keys()
}
# Stricter threshold for clinopyroxenes due to more endmembers
ENDMEMBER_DISTANCE_THRESHOLD["Cpx"] = 0.35

def _get_element_value(phase: dict, element: str) -> float:
    """Get element value from phase dict, trying multiple key patterns."""
    # Try simple element name first
    if element in phase:
        val = phase[element]
        if val is not None:
            try:
                return float(val)
            except (ValueError, TypeError):
                pass
    
    # Try with _apfu suffix
    key_apfu = f"{element}_apfu"
    if key_apfu in phase:
        val = phase[key_apfu]
        if val is not None:
            try:
                return float(val)
            except (ValueError, TypeError):
                pass
    
    return 0.0

def discriminate_solvus(phases):
    """
    Group phases by parent phase name and apply endmember classification.
    For phases with defined endmember compositions, apply compositional 
    discrimination regardless of whether multiple phases exist.
    
    Returns:
        Tuple of (processed_phases, transformation_stats)
        where transformation_stats is dict: {(parent, assigned): count}
    """
    grouped = {}
    for p in phases:
        name = p.get("phase", "")
        grouped.setdefault(name, []).append(p)

    output = []
    transformations = {}  # Track (original_name, assigned_name) -> count
    
    for parent_name, group in grouped.items():
        parent_key = None
        for k in ENDMEMBER_COMPOSITIONS:
            if k.lower() == parent_name.lower():
                parent_key = k
                break

        if parent_key is None:
            logger.debug("No endmember rules for phase '%s', keeping as-is", parent_name)
            output.extend(group)
            # Count unchanged phases
            for phase in group:
                key = (parent_name, parent_name)
                transformations[key] = transformations.get(key, 0) + 1
            continue

        logger.debug("Applying endmember classification to %d phase(s) of '%s'", len(group), parent_name)
        
        if parent_key == "Cpx":
            logger.debug("  Using omphacite heuristic for Cpx phases")
            discriminated = _discriminate_cpx(group)
            output.extend(discriminated)
            # Count transformations
            for phase in discriminated:
                assigned = phase.get("phase")
                key = (parent_name, assigned)
                transformations[key] = transformations.get(key, 0) + 1
            continue
        
        # Special handling for feldspars
        if parent_key == "Fsp":
            logger.debug("  Using Ca-alkali criteria for Fsp discrimination")
            for phase in group:
                phase["orig_phase"] = phase.get("phase")
                
                na = _get_element_value(phase, "Na")
                ca = _get_element_value(phase, "Ca")
                k = _get_element_value(phase, "K")
                
                total_alkali = na + ca + k
                if total_alkali <= 0:
                    assigned_name = "Fsp"
                elif ca / total_alkali > 0.10:
                    assigned_name = "Pl"
                elif na > k:
                    assigned_name = "Pl"
                else:
                    assigned_name = "Afs"
                
                phase["phase"] = assigned_name
                key = (parent_name, assigned_name)
                transformations[key] = transformations.get(key, 0) + 1
                output.append(phase)
            continue

        # Standard endmember distance-based discrimination for ALL phases
        threshold = ENDMEMBER_DISTANCE_THRESHOLD.get(parent_key, 0.3)
        for idx, phase in enumerate(group):
            phase["orig_phase"] = phase.get("phase")
            best_em, best_dist = _find_nearest_endmember(phase, parent_key)
            
            logger.debug("  Phase %d/%d: nearest endmember='%s', distance=%.4f (threshold=%.2f)", 
                        idx+1, len(group), best_em, best_dist, threshold)
            
            if best_dist <= threshold:
                old_name = phase["phase"]
                phase["phase"] = best_em
                logger.debug("    Assigned '%s' -> '%s' (distance %.4f <= threshold %.2f)", 
                           old_name, best_em, best_dist, threshold)
                # Count transformation
                key = (parent_name, best_em)
                transformations[key] = transformations.get(key, 0) + 1
            else:
                logger.debug("    Keeping original name '%s' (distance %.4f > threshold %.2f)", 
                           phase["phase"], best_dist, threshold)
                # Count unchanged
                key = (parent_name, parent_name)
                transformations[key] = transformations.get(key, 0) + 1
            output.append(phase)

    return output, transformations


def _discriminate_cpx(cpx_group):
    """
    Discriminate clinopyroxenes using compositional heuristics.
    
    Classifies clinopyroxenes based on the standard nomenclature following
    Morimoto et al. (1988) pyroxene classification scheme. Recognizes six
    endmember types plus omphacite as an intermediate composition:
    
    - Jd (jadeite): Na-Al pyroxene, very high jadeite component (>0.8), very low Ca (<0.15)
    - Omp (omphacite): Intermediate Di-Jd solid solution, significant jadeite component (>0.25)
    - Acm (aegirine): Na-Fe pyroxene, high Na, low Ca, Fe-dominant
    - Di (diopside): High-Ca Mg-rich pyroxene (Ca>0.5, Mg#>0.6)
    - Hd (hedenbergite): High-Ca Fe-rich pyroxene (Ca>0.5, Mg#<0.4)
    - Aug (augite): Intermediate-Ca pyroxene (0.2<Ca<0.5)
    - Pig (pigeonite): Low-Ca pyroxene (0.05<Ca<0.2)
    - Cpx (generic): Compositions not matching specific endmember criteria
    
    Returns list of phases with updated names.
    """
    output = []
    
    for phase in cpx_group:
        phase["orig_phase"] = phase.get("phase")
        
        # Get cation abundances
        na = _get_element_value(phase, "Na")
        al = _get_element_value(phase, "Al")
        ca = _get_element_value(phase, "Ca")
        mg = _get_element_value(phase, "Mg")
        fe = _get_element_value(phase, "Fe")
        
        logger.debug("    Cpx composition: Ca=%.3f, Mg=%.3f, Fe=%.3f, Na=%.3f, Al=%.3f", 
                    ca, mg, fe, na, al)
        
        # Calculate diagnostic ratios
        total_cations = ca + mg + fe + na + al
        if total_cations <= 0:
            phase["phase"] = "Cpx"
            logger.debug("  No cation data available, keeping as 'Cpx'")
            output.append(phase)
            continue
        
        # Normalized fractions
        jd_component = (na + al) / total_cations  # Jadeite component
        ca_frac = ca / total_cations              # Calcium content
        na_frac = na / total_cations              # Sodium content
        mg_number = mg / (mg + fe) if (mg + fe) > 0 else 0.0  # Mg/(Mg+Fe)
        
        # Hierarchical classification based on standard pyroxene nomenclature
        
        # First branch: Sodic pyroxenes (high Na+Al content)
        
        # 1) Jadeite: Pure sodic pyroxene with very high jadeite component and minimal Ca
        # Corresponds to compositions near the Jd apex in the Jd-Ae-Di+Hd diagram
        if jd_component > 0.8 and ca_frac < 0.15:
            assigned_name = "Jd"
            logger.debug("  Assigned 'Jd' (jadeite): jd_component=%.3f, ca_frac=%.3f", 
                        jd_component, ca_frac)
        
        # 2) Omphacite: Jadeite with significant Ca, representing Di-Jd solid solution
        # Occupies the region between pure jadeite and the diopside-hedenbergite join
        # Critical for high-pressure metamorphic assemblages (eclogites)
        elif jd_component > 0.25:
            assigned_name = "Omp"
            logger.debug("  Assigned 'Omp' (omphacite): jd_component=%.3f, ca_frac=%.3f", 
                        jd_component, ca_frac)
        
        # 3) Aegirine: Sodic Fe-rich pyroxene with high Na, low Ca, Fe-dominant
        # Forms solid solution toward aegirine-augite with increasing Ca
        elif na_frac > 0.4 and ca_frac < 0.25 and mg_number < 0.3:
            assigned_name = "Acm"
            logger.debug("  Assigned 'Acm' (aegirine): na_frac=%.3f, ca_frac=%.3f, mg_number=%.3f", 
                        na_frac, ca_frac, mg_number)
        
        # Second branch: Non-sodic pyroxenes (low Na+Al, classified by Ca content)
        
        # 4) Diopside: High-Ca, Mg-rich endmember
        # Requires Ca > 0.5 (>50% of M2 site) and high Mg# characteristic of diopside
        elif ca_frac > 0.5 and mg_number > 0.6:
            assigned_name = "Di"
            logger.debug("  Assigned 'Di' (diopside): ca_frac=%.3f, mg_number=%.3f", 
                        ca_frac, mg_number)
        
        # 5) Hedenbergite: High-Ca, Fe-rich endmember
        # Requires Ca > 0.5 and Fe-dominant composition
        elif ca_frac > 0.5 and mg_number < 0.4:
            assigned_name = "Hd"
            logger.debug("  Assigned 'Hd' (hedenbergite): ca_frac=%.3f, mg_number=%.3f", 
                        ca_frac, mg_number)
        
        # 6) Augite: Intermediate-Ca pyroxene
        # Standard definition: 20-45% Ca in M2 site (Wo20-45 in Wo-En-Fs diagram)
        # Characteristic of many igneous and metamorphic rocks
        elif 0.2 <= ca_frac <= 0.5:
            assigned_name = "Aug"
            logger.debug("  Assigned 'Aug' (augite): ca_frac=%.3f, mg_number=%.3f", 
                        ca_frac, mg_number)
        
        # 7) Pigeonite: Low-Ca pyroxene
        # Standard definition: 5-20% Ca in M2 site (Wo5-20)
        # Monoclinic structure distinguishes it from orthopyroxene (Wo < 5)
        # Common in rapidly cooled igneous rocks and some metamorphic assemblages
        elif 0.05 <= ca_frac < 0.2:
            assigned_name = "Pig"
            logger.debug("  Assigned 'Pig' (pigeonite): ca_frac=%.3f, mg_number=%.3f", 
                        ca_frac, mg_number)
        
        # 8) Generic Cpx: Very low Ca (<5%) or compositions not matching criteria
        # Very low Ca compositions (Ca < 0.05) approach orthopyroxene field
        # but retain monoclinic structure, hence generic Cpx designation
        else:
            assigned_name = "Cpx"
            logger.debug("  Kept as 'Cpx' (intermediate/low-Ca composition): jd_comp=%.3f, ca_frac=%.3f, mg_number=%.3f", 
                        jd_component, ca_frac, mg_number)
        
        phase["phase"] = assigned_name
        output.append(phase)
    
    return output


def _find_nearest_endmember(phase: dict, parent_name: str):
    endmembers = ENDMEMBER_COMPOSITIONS[parent_name]

    phase_elems = set()
    for sto in endmembers.values():
        phase_elems |= set(sto.keys())
    elems = tuple(sorted(phase_elems))

    key = (parent_name, elems)
    cached = _NORMALIZED_ENDMEMBERS_ARRAY.get(key)
    if cached is not None:
        em_names, em_array = cached
    else:
        em_names = list(endmembers.keys())
        em_vecs = []
        for nm in em_names:
            sto = endmembers[nm]
            vec = np.array([sto.get(e, 0.0) for e in elems], dtype=float)
            total = vec.sum()
            if total <= 0:
                total = 1.0
            vec = vec / total
            em_vecs.append(vec)
        em_array = np.vstack(em_vecs)
        _NORMALIZED_ENDMEMBERS_ARRAY[key] = (em_names, em_array)

    p_vec = np.array([_get_element_value(phase, e) for e in elems], dtype=float)
    p_total = p_vec.sum()
    if p_total <= 0:
        p_total = 1.0
    p_vec = p_vec / p_total

    diffs = em_array - p_vec
    dists = np.linalg.norm(diffs, axis=1)
    idx = int(dists.argmin())
    
    phase_comp_str = ", ".join([f"{e}={_get_element_value(phase, e):.3f}" for e in elems if _get_element_value(phase, e) > 0.001])
    logger.debug("      Phase composition (raw): %s", phase_comp_str)
    logger.debug("      Distances to endmembers: %s", 
                ", ".join([f"{em_names[i]}={dists[i]:.4f}" for i in range(len(em_names))]))
    
    return em_names[idx], float(dists[idx])