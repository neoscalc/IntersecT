# IntersecT Technical Documentation

This document provides technical information about IntersecT's implementation for developers extending the system, advanced users debugging issues, and researchers seeking implementation details. The package structure and algorithms follow the methodology described in Nerone et al. (2025, https://doi.org/10.1016/j.cageo.2025.105949).

## Architecture Overview

IntersecT consist of data parsing, phase identification and quality factor calculations. The parsing module handles file formats from Perple_X and MAGEMin and extracts compositional data with appropriate stoichiometric coefficients. The mineral aliases module resolves phase names from software-specific abbreviations to Warr (2021) standard nomenclature and applies compositional discrimination to solid solution series. The IntersecT module integrates these components to construct data tables matching observed and predicted compositions, then calculates quality factors with statistical weighting based on reduced χ² statistics.

New thermodynamic software packages can be supported by extending the parser without modifying phase resolution or calculation logic. Additional solid solution series can be added to the discrimination system by updating endmember definitions without touching parsing or quality factor code.

## Parsing System (parsers.py)

The parsing system handles thermodynamic model outputs from multiple software packages with different file formats and organizational structures. The parser detects file format, identifies coordinate columns representing independent variables of the calculation, extracts phase-specific compositional data, and handles stoichiometric coefficients in oxide notation.

### File Format Detection

The parser examines the first non-empty line to determine file format. Files containing two or more commas are treated as CSV format characteristic of MAGEMin output and processed with comma separation. Files lacking this comma density are assumed to follow Perple_X column format and processed with whitespace separation.

For Perple_X files, the parser searches for a header line containing "Name" with coordinate identifiers like "T(K)" or "P(bar)". This header defines the column structure for subsequent data rows. If no header is found, the first non-comment line is treated as the header to accommodate variant output formats.

### Coordinate Column Identification

Coordinate columns define the independent variables of the thermodynamic calculation. The system recognizes temperature columns through variants including "T(K)", "T", "T[°C]", and "temperature", pressure columns through "P(bar)", "P[kbar]", "P", and "pressure", and compositional variables through "X", "x", "X[0.0-1.0]", and similar patterns. All matching operates case-insensitively to handle variations across software versions.

When all three coordinate types are present, as occurs in MAGEMin calculations where bulk composition varies as an independent variable alongside pressure and temperature, the parser returns all three to enable subsequent structural analysis. For files with two varying coordinates, the parser returns the identified pair. If automatic detection fails, the parser uses the first two non-excluded columns where exclusions include "phase", "name", "point", and "counter".

### Element Information Extraction

Compositional data in thermodynamic model outputs appears in oxide notation with stoichiometric coefficients embedded in column names. The system identifies these patterns to separate the base element symbol from its coefficient. The regular expression pattern `^([A-Z][a-z]?)(\d*)O\d*` matches oxide notation such as Al₂O₃, Na₂O, K₂O, FeO, MgO, and CaO, capturing the element symbol and coefficient. When no coefficient appears before oxygen, the value defaults to one as in FeO, MgO, and CaO.

Alternative patterns include bracket notation like Mg[apfu] and parenthetical notation such as Si(apfu). These typically indicate direct elemental abundances rather than oxide compositions and receive a coefficient of one. For ambiguous cases, the system extracts the first one or two characters as the element symbol with default coefficient of one.

Element names are normalized through a lookup table handling common variations. This ensures sodium appears as Na rather than na, aluminum as Al rather than al, and similarly for other elements. Normalization operates independently of coefficient extraction, maintaining separation between identification of element presence and determination of stoichiometric contribution.

### Compositional Data Filtering

Thermodynamic model outputs contain numerous columns representing thermodynamic properties, elastic moduli, seismic velocities, and other derived quantities that should not be treated as compositional variables. The filtering system identifies and excludes these columns. Explicit exclusion patterns include Greek letter identifiers for thermodynamic derivatives, heat capacity indicators like cp and cv, entropy columns, elastic properties including "Gruneisen", "Ks,", and "Gs,", thermodynamic potentials marked by "V,J/", "H,J/", and "G,J/", and seismic properties containing "Vp,", "Vs,", or "rho,".

The system requires columns to contain either oxide notation with formula unit indicators or explicit apfu labeling. Columns must contain the letter O for oxygen with indicators like "mol", "pfu", or "apfu", or they must contain explicit markers like "[apfu]" or "(apfu)". Weight percent columns and mode percentages are excluded. Chemical potential columns identified by "mu[" are excluded as they represent thermodynamic derivatives rather than compositional abundances.

### Data Block Construction

The parser groups rows sharing identical coordinate values into calculation blocks representing single points in pressure-temperature-composition space. Each block contains a coordinate dictionary mapping variable names to numerical values and a list of phase dictionaries representing the stable mineral assemblage at those conditions.

Phase dictionaries include raw column data preserving original naming conventions and processed element data with normalized names and applied stoichiometric coefficients. The raw data preserves complete information for diagnostic examination. The processed data with normalized element names and multiplied coefficients enables direct comparison with observed compositions without subsequent coefficient handling.

System detection of file origin (MAGEMin or Perple_X) is preserved through a "system" field that downstream processes use to select appropriate name resolution mappings. This system tagging occurs during initial format detection and propagates through all processing stages.

## Phase Name Resolution and Solvus Discrimination (mineral_aliases.py)

The system translates software-specific abbreviations to standardized nomenclature following Warr (2021) and distinguishes chemically distinct phases within solid solution series based on compositional criteria. The implementation separates name resolution, which maps variant abbreviations to standard nomenclature, from endmember discrimination, which assigns specific compositional varieties based on measured cation abundances.

### Name Resolution Mechanism

The MineralAliases class loads translation tables from a TOML configuration file containing mappings from software-specific abbreviations to Warr (2021) standard nomenclature. Separate tables exist for MAGEMin and Perple_X because these packages employ different abbreviation conventions. All mappings use normalized lowercase keys for case-insensitive matching to accommodate variations in capitalization across file formats and software versions.

Resolution queries specify phase name and source system, allowing selection of the appropriate translation table. When a mapping exists, the resolver returns the Warr (2021) standard abbreviation. In strict mode, unmapped names raise exceptions to alert users that the configuration lacks entries for phases in their model output. In non-strict mode, unmapped names pass through unchanged with optional persistence for tracking.

Case-insensitive matching handles variations where different Perple_X versions may capitalize phase names as GARNET, Garnet, or garnet. Normalization to lowercase keys during table loading and query formulation ensures all variants map correctly to the Grt abbreviation.

### Endmember Discrimination Framework

Phases sharing a parent name after initial resolution undergo compositional analysis for endmember assignment. The discrimination process applies to all phases with defined endmember compositions in the database, operating independently on each phase regardless of whether other phases with the same parent name exist at the calculation point. This approach enables consistent classification based on absolute compositional criteria rather than relative comparisons between coexisting phases.

The endmember database defines stoichiometric ratios for reference compositions used in discrimination calculations. Each solid solution series contains a dictionary mapping endmember names to characteristic cation abundances. The definitions focus on major elements that distinguish compositional types while omitting minor components that would complicate calculations without improving discrimination. These compositions represent idealized endmember stoichiometries as they appear in thermodynamic databases rather than natural mineral compositions that typically show solid solution.

For white micas, the database distinguishes paragonite with Na=1, K=0, Al=3, Si=3, muscovite with Na=0, K=1, Al=3, Si=3, and margarite with Ca=1, Al=4, Si=2. Feldspars separate plagioclase with Ca=1, Al=2, Si=2 from alkali feldspar with K=1, Al=1, Si=3. Amphiboles differentiate tremolite with Ca=2, Mg=5, Si=8, actinolite with Ca=2, Mg=3, Fe=2, Si=8, and glaucophane with Na=2, Mg=3, Al=2, Si=8. Spinels include chromite at Fe=1, Cr=2, ulvöspinel at Fe=2, Ti=1, and magnetite at Fe=3. The ilmenite-hematite series distinguishes ilmenite at Fe=1, Ti=1 from hematite at Fe=2. Carbonates differentiate calcite at Ca=1, dolomite at Ca=1, Mg=1, ankerite at Ca=1, Fe=1, and siderite at Fe=1.

Clinopyroxenes require special consideration in both database structure and discrimination methodology. The endmember compositions define six endmembers representing the vertices of the clinopyroxene compositional space: diopside at Ca=1, Mg=1, Si=2, hedenbergite at Ca=1, Fe=1, Si=2, augite at Ca=0.7, Mg=0.7, Fe=0.3, Al=0.3, Si=2, pigeonite at Ca=0.1, Mg=0.9, Fe=0.1, Si=2, jadeite at Na=1, Al=1, Si=2, and aegirine at Na=1, Fe=1, Si=2. These compositions establish the framework within which clinopyroxene classification operates. Notably absent from this database is omphacite, which represents a compositionally intermediate series along the join between diopside and jadeite rather than an independent endmember. Natural omphacites contain variable proportions of both diopside and jadeite components, typically expressed as compositions with significant contents of both calcium-magnesium and sodium-aluminum cations.

The inclusion of clinopyroxene endmembers in the database serves primarily to maintain consistency with the general database structure and to document the compositional framework that guides discrimination. However, clinopyroxenes undergo compositional heuristic discrimination as described in the following section.

### Distance-Based Discrimination for Most Phases

For phases other than clinopyroxenes, the system applies distance-based discrimination using Euclidean distances in normalized compositional space. The algorithm computes normalized compositional vectors for each phase and compares them to endmember stoichiometries defined in the database. Normalization divides each element abundance by the sum of relevant cations, creating unit-normalized vectors that compare stoichiometric ratios rather than absolute concentrations. This normalization eliminates effects of analytical totals and formula unit conventions, ensuring that compositional comparisons reflect actual cation proportions.

Distance thresholds determine whether composition matches an endmember sufficiently for name assignment. Most systems employ a threshold of 0.3 in normalized compositional space, corresponding to compositional differences that remain within typical solid solution ranges observed in natural minerals and thermodynamic models. When multiple endmembers exist for a parent phase, the system identifies the nearest endmember by Euclidean distance and assigns that name if the distance falls below the threshold. Phases with distances exceeding the threshold for all defined endmembers retain the parent phase name, indicating compositional characteristics that place them between idealized endmember compositions.

The distance calculation incorporates all elements present in either the phase composition or the endmember definitions, constructing vectors in a multidimensional space where each dimension represents a different cationic species. This comprehensive approach ensures that discrimination considers the complete compositional signature rather than focusing on a subset of diagnostic elements. The system caches normalized endmember arrays to minimize redundant calculations when processing datasets with many calculation points containing the same phase assemblages, improving computational efficiency without affecting results.

White micas exemplify this discrimination approach. When the system encounters a phase resolved as Mca (white mica) from the software output, it extracts the sodium, potassium, calcium, aluminum, and silicon abundances, constructs a normalized vector, and calculates distances to the paragonite, muscovite, and margarite reference compositions. A mica with high potassium relative to sodium will show small distance to muscovite and large distance to paragonite, receiving the muscovite endmember name if the distance falls below 0.3. A mica with intermediate sodium-potassium ratios that produces distances exceeding 0.3 to all three endmembers retains the generic Mca name, indicating a composition between idealized endmembers. Amphiboles follow similar logic with distance calculations to tremolite, actinolite, and glaucophane determining which calcic or sodic amphibole variety best represents the observed composition.

### Clinopyroxene Specialized Discrimination

Clinopyroxenes receive specialized handling through compositional heuristics rather than distance-based matching. The classification distinguishes eight compositional varieties including six true thermodynamic endmembers plus omphacite as a petrologically critical intermediate composition and a generic clinopyroxene category for compositions not matching specific criteria.

The discrimination algorithm applies a hierarchical cascade of compositional criteria to each clinopyroxene individually, evaluating conditions from most to least restrictive. The system calculates diagnostic ratios including jadeite component as the sum of sodium and aluminum divided by total major cations, calcium fraction as calcium divided by total cations, sodium fraction as sodium divided by total cations, and magnesium number as magnesium divided by the sum of magnesium and iron. These ratios capture the compositional variations that distinguish different clinopyroxene varieties in both petrological practice and standard nomenclature schemes.

The classification proceeds through eight hierarchical conditions organized into two primary branches. The first branch addresses sodic pyroxenes characterized by high sodium-aluminum content. Jadeite requires very high jadeite component exceeding 0.8 combined with very low calcium fraction below 0.15, capturing pure sodic pyroxenes characteristic of ultra-high pressure metamorphism where the jadeite endmember dominates. Omphacite requires jadeite component exceeding 0.25 without stringent calcium restrictions, recognizing all sodic-calcic pyroxenes with significant jadeite component typical of eclogite-facies assemblages. This threshold-based recognition of omphacite addresses the fundamental challenge that omphacite occupies compositional space between diopside and jadeite endmembers while requiring explicit identification for petrological interpretation. The classification of omphacite as jadeite with significant calcium content directly reflects the triangular compositional diagram showing omphacite positioned between the jadeite apex and the diopside-hedenbergite join. Aegirine requires high sodium fraction exceeding 0.4, low calcium fraction below 0.25, and iron-dominant character with magnesium number below 0.3, identifying sodic ferric pyroxenes associated with alkaline igneous rocks or metasomatized assemblages.

The second branch addresses non-sodic pyroxenes classified according to calcium content following the standard quadrilateral classification scheme. Diopside requires high calcium fraction exceeding 0.5 combined with magnesium-dominant character indicated by magnesium number exceeding 0.6, recognizing calcic magnesian pyroxenes typical of metabasites and metacarbonates. Hedenbergite requires high calcium fraction exceeding 0.5 with iron-dominant character indicated by magnesium number below 0.4, identifying calcic ferrous pyroxenes that appear in iron-rich metamorphic assemblages. Augite requires intermediate calcium fraction between 0.2 and 0.5, corresponding to the standard nomenclature definition of twenty to forty-five percent wollastonite component. Augite represents one of the most common clinopyroxene varieties in igneous and metamorphic rocks, occupying the central region of the pyroxene quadrilateral diagram with intermediate calcium content and variable magnesium-iron ratios. Pigeonite requires low calcium fraction between 0.05 and 0.2, corresponding to five to twenty percent wollastonite component in the standard classification. Pigeonite forms a monoclinic low-calcium pyroxene distinct from orthopyroxenes that contain less than five percent wollastonite component, commonly appearing in rapidly cooled igneous rocks and high-temperature metamorphic assemblages. Compositions satisfying none of these specific criteria, including those with very low calcium content below 0.05 that approach the orthopyroxene compositional field while retaining monoclinic structure, receive the generic clinopyroxene designation.

This hierarchical structure ensures that the most restrictive conditions evaluate first, preventing misclassification of rare endmember compositions. The jadeite condition must be satisfied before omphacite consideration, guaranteeing that nearly pure jadeite receives its proper name rather than being grouped with jadeite-bearing omphacite. Similarly, omphacite recognition precedes evaluation of other sodic varieties like aegirine, ensuring that sodic-calcic compositions receive appropriate classification before considering purely sodic alternatives. Among non-sodic pyroxenes, the progression from high-calcium varieties through intermediate-calcium augite to low-calcium pigeonite follows decreasing calcium content, matching the logical flow of standard pyroxene classification diagrams.

Users can specify any of the eight recognized clinopyroxene varieties in their measured composition files and the system will successfully match these specifications with model phases classified automatically through the compositional heuristic. This capability enables quality factor calculations for assemblages where specific clinopyroxene varieties provide key constraints on pressure-temperature-composition conditions. For instance, users analyzing eclogites can specify omphacite compositions, those working with basaltic systems can specify augite, and those studying rapidly cooled igneous rocks can specify pigeonite, with the system appropriately identifying and matching the corresponding phases predicted by the thermodynamic model.

### Computational Efficiency

The system caches normalized endmember arrays to minimize redundant calculations when processing datasets with many calculation points containing the same phase assemblages. The first encounter with a parent phase and element set computes normalized vectors for all endmembers and stores them in a global cache keyed by parent name and element tuple. Subsequent calculations retrieve cached vectors rather than recomputing them. This caching proves effective when processing regular grids where the same assemblage appears at many points, with arrays computed once during processing of the first point and reused hundreds or thousands of times. The cache persists within a single analysis run but not between executions, requiring recomputation when the program restarts.

The discrimination system operates on phases that have already undergone name resolution, receiving as input phases with standardized Warr (2021) nomenclature rather than software-specific abbreviations. This separation of concerns enables independent extension of translation tables without modifying discrimination logic, and conversely allows refinement of discrimination criteria without touching name resolution. The modular architecture facilitates maintenance and future development while ensuring that each component performs a well-defined function in the overall phase identification pipeline.

## Quality Factor Analysis Framework (IntersecT.py)

The QualityFactorAnalysis class integrates parsing, phase resolution, and discrimination to construct data tables matching observed and predicted compositions, then calculates quality factors with statistical weighting following the methodology described in Nerone et al. (2025).

### Data Structure Integration

The class maintains dual structures representing modern and legacy formats to enable incremental migration of calculation methods while preserving compatibility with original algorithms. The modern structure stores parsed blocks as a list of dictionaries with coordinate values, phase compositions, and system identifiers. The legacy structure flattens this into a pandas DataFrame matching original organization.

Construction of the legacy structure occurs through build_intersect_table, which iterates through blocks extracting coordinates and compositions matching observed phase-element combinations. For each block, the method resolves phase names, applies solvus discrimination, searches for matching phases, and extracts element values, recording NaN for absent phases.

The DataFrame ensures coordinate columns appear first followed by composition columns in user-specified order. The labels array stores column names as a NumPy array for advanced indexing. The phase_id array assigns integer identifiers grouping columns by parent phase, enabling iteration over phases to compute phase-specific quality factors.

### Coordinate System Configuration

The system accommodates pressure-temperature, composition-pressure, and temperature-composition diagrams. The suggest_plot_coordinates method analyzes block structure to identify varying coordinates and organizational patterns, then returns appropriate axis assignments.

For files with two varying coordinates, the method applies priority logic preferring temperature-pressure combinations. Files with three varying coordinates require structural analysis through is_independent_variable, which examines whether a coordinate exhibits the organizational pattern of an independent variable characterized by multiple combinations of other coordinates for each unique value.

When analysis identifies composition as independent, the system recommends composition horizontally with pressure or temperature vertically. When analysis indicates pressure-temperature organization, the system recommends temperature horizontally and pressure vertically. Users can override recommendations through set_plot_coordinates.

Unit conversions transform temperatures from Kelvin to Celsius by subtracting 273 and pressures from bars to gigapascals by dividing by 10000.

### Quality Factor Calculation Pipeline

Quality factors quantify agreement between observed and predicted compositions while accounting for analytical uncertainties. Calculation proceeds through element-specific, phase-specific, and total quality factors. Element-specific values evaluate individual compositional variables. Phase-specific values average within each mineral. Total values integrate across all phases.

The formula implements Duesterhoeft and Lanari (2020) methodology incorporating an uncertainty-normalized difference raised to a power dependent on predicted value. Power dependence ensures small values near detection limits do not dominate. Uncertainty normalization accounts for varying analytical precision.

Uncertainties enter from user-specified values or automatic calculation based on technique. The calc_obs_err method implements empirical relationships between compositional magnitude and precision for electron microprobe analyses derived from Lanari and Duesterhoeft (2019) and Lanari and Hermann (2021). Separate parameterizations exist for EDS, WDS mapping, and WDS spot analyses.

### Statistical Weighting System

Reduced χ² statistics provide independent assessment of model-data agreement complementing quality factors. Values near unity suggest differences consistent with analytical uncertainty alone. Values substantially exceeding unity indicate systematic misfit potentially arising from disequilibrium, analytical bias, or thermodynamic database limitations.

Weighting assigns each phase a weight inversely proportional to minimum reduced χ², giving greater influence to phases with good agreement and down-weighting those with poor fits. Minimum values below unity are reset to one before inversion to prevent artificial over-weighting. Weights are normalized to sum to unity.

Application of weights to phase-specific quality factors yields the weighted total quality factor representing the final statistical estimate of best fit conditions. Maximum locations identify conditions where the model best reproduces observations while accounting for phase reliability. Comparison between weighted and unweighted results reveals whether phases indicate consistent conditions or whether systematic disagreements suggest disequilibrium.

## Extension Guidelines

Supporting additional thermodynamic software requires implementing parsing logic for the specific file format. The parse_table function provides the extension point for format detection and parsing branches. New parsers should return data in the standard block format with coordinate dictionaries and phase lists for direct integration with existing phase resolution and calculation code.

Adding solid solution series requires defining endmember compositions in ENDMEMBER_COMPOSITIONS and setting thresholds in ENDMEMBER_DISTANCE_THRESHOLD. Definitions should specify characteristic cation ratios for each endmember, focusing on major elements distinguishing compositional types while omitting minor components. Thresholds should be calibrated through examination of model outputs to find values correctly discriminating known compositions without excessive sensitivity to minor variations.

For solid solutions requiring specialized discrimination beyond distance calculations, the discriminate_solvus function provides the extension point for custom algorithms. The clinopyroxene omphacite discrimination demonstrates the pattern of checking parent names and delegating to specialized functions implementing appropriate compositional heuristics for specific mineral systems.