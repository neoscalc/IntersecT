from pathlib import Path
import os
import sys
import logging
from tkinter import filedialog
from typing import Sequence, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from .parsers import parse_table, extract_element_info, normalize_element_name
from .mineral_aliases import MineralAliases, discriminate_solvus


class QualityFactorAnalysis:
    
    def __init__(self, mineral_symbols_toml: Optional[Path] = None):
        # Core data structures
        self.data = pd.DataFrame()  # Legacy compatibility
        self.labels = []
        self.x = np.array([])
        self.y = np.array([])
        
        # Phase tracking
        self.labels_new = []
        self.phase_id = np.array([])
        self.phase_name = []
        
        # Observations
        self.apfu_name = []
        self.apfu_obs = np.array([])
        self.obs_err = np.array([])
        self.analysis_type = ''
        
        # Output settings
        self.color_scheme = 'viridis'
        self.output_dir = ''
        self.min_redchi2 = np.array([])
        
        # Modern parsing support
        self.model_blocks: list[dict] = []
        self.coord_columns: Optional[list[str]] = None
        self.coords_x: Optional[str] = None
        self.coords_y: Optional[str] = None
        self.phase_data: Optional[pd.DataFrame] = None
        
        # Mineral aliases
        if mineral_symbols_toml is None:
            mineral_symbols_toml = Path(__file__).parent / "mineral_symbols.toml"
        self.aliases = MineralAliases(mineral_symbols_toml)
        

    @classmethod
    def from_default_symbols(cls):
        """Create instance with default mineral symbols file."""
        return cls()

    def _setup_logging(self):
        if not self.output_dir:
            return
        log_path = os.path.join(self.output_dir, "log_IntersecT.txt")
        logging.basicConfig(
            filename=log_path,
            filemode='w',
            encoding='utf-8',
            level=logging.INFO,
            format='%(message)s',
            force=True
        )
    
    def _log_print(self, message):
        print(message)
        logging.info(message)

    # ================================
    # I/O Methods
    # ================================
    
    def set_output_directory(self, path: Optional[str] = None):
        """Prompt user to select an output directory."""
        if path is None:
            self._log_print("Please select an output directory.")
            self.output_dir = filedialog.askdirectory()
        else:
            self.output_dir = str(Path(path).expanduser().resolve())
        
        self._log_print("The output directory is:  " + self.output_dir)
        self._setup_logging()

    def load_model_output(self, coord_columns: Optional[Sequence[str]] = None, 
                         filename: Optional[str] = None):
        """Load model output and apply phase name resolution with solvus discrimination."""
        if filename is None:
            self._log_print("Please select a csv/phm file from MAGEMin/Perple_X.")
            filename = filedialog.askopenfilename()
            if not filename:
                return

        self._log_print("The data file is:  " + filename)

        self.model_blocks = parse_table(filename, coord_columns=coord_columns)

        if self.model_blocks:
            fc = list(self.model_blocks[0].get("coords", {}).keys())
            self.coord_columns = fc if fc else (list(coord_columns) if coord_columns else None)

            # Apply phase name resolution and solvus discrimination to all blocks
            total_transformations = {}

            for block in self.model_blocks:
                system = block.get("system", "perplex")
                phases = []

                for ph in block.get("phases", []):
                    phc = ph.copy()
                    original_phase_name = str(phc.get("phase", "")).strip()

                    try:
                        resolved_name = self.aliases.resolve(original_phase_name, system)
                        phc["phase"] = resolved_name
                        phc["orig_phase"] = original_phase_name
                    except KeyError:
                        phc["orig_phase"] = original_phase_name

                    phases.append(phc)

                # Apply solvus discrimination and collect statistics
                phases, block_transformations = discriminate_solvus(phases)

                # Update block with processed phases
                block["phases"] = phases

                # Aggregate transformations
                for key, count in block_transformations.items():
                    total_transformations[key] = total_transformations.get(key, 0) + count

            # Log transformation statistics
            self._log_print("\n" + "-"*60)
            self._log_print("Phase discrimination summary:")

            # Group by parent phase and sort
            parent_phases = {}
            for (parent, assigned), count in total_transformations.items():
                if parent not in parent_phases:
                    parent_phases[parent] = []
                parent_phases[parent].append((assigned, count))

            for parent in sorted(parent_phases.keys()):
                assignments = sorted(parent_phases[parent], key=lambda x: x[1], reverse=True)
                total_count = sum(count for _, count in assignments)

                for assigned, count in assignments:
                    if parent == assigned:
                        self._log_print(f"  {parent} → {assigned} : {count} (unchanged)")
                    else:
                        percentage = (count / total_count) * 100
                        self._log_print(f"  {parent} → {assigned} : {count} ({percentage:.1f}%)")

            self._log_print("-"*60 + "\n")
        

    def suggest_plot_coordinates(self):
        """Auto-suggest two coordinates for plotting."""
        if not self.model_blocks or not self.coord_columns:
            return None, None

        # Priority: T then P, or X then P
        temp_candidates = ["T(K)", "T", "T[°C]", "temperature"]
        pres_candidates = ["P(bar)", "P[kbar]", "P", "pressure"]
        x_candidates = ["X", "x", "X[0.0-1.0]", "X(C1)", "x(c1)", "Y(CO2)", "y(co2)"]

        def pick(cands, cols):
            for cand in cands:
                for c in cols:
                    if c.lower() == cand.lower():
                        return c
            return None

        def n_unique(col_name):
            """Count unique values for a coordinate across all blocks."""
            vals = set()
            for block in self.model_blocks:
                v = block.get("coords", {}).get(col_name)
                if v is not None:
                    vals.add(v)
            return len(vals)

        def is_independent_variable(var_name, other_vars):
            """
            Check if var_name is an independent variable by examining data structure.
            Returns True if for each unique value of var_name, there are multiple 
            combinations of other variables (indicating var_name is outer loop).
            """
            # Group data by var_name and collect unique combinations of other vars
            var_groups = {}
            for block in self.model_blocks:
                coords = block.get("coords", {})
                if var_name not in coords:
                    continue
                
                var_val = coords[var_name]
                other_vals = tuple(coords.get(ov) for ov in other_vars if ov in coords)

                if var_val not in var_groups:
                    var_groups[var_val] = set()
                var_groups[var_val].add(other_vals)

            # If most var_name values have multiple other_var combinations, 
            # var_name is likely an independent variable
            if not var_groups:
                return False

            avg_combinations = sum(len(combos) for combos in var_groups.values()) / len(var_groups)
            return avg_combinations > 1.5

        t = pick(temp_candidates, self.coord_columns)
        p = pick(pres_candidates, self.coord_columns)
        x = pick(x_candidates, self.coord_columns)

        # Check variance for each coordinate
        nt = n_unique(t) if t else 0
        np_val = n_unique(p) if p else 0
        nx = n_unique(x) if x else 0

        # Filter out non-varying coordinates
        varying = {}
        if nt > 1:
            varying['t'] = t
        if np_val > 1:
            varying['p'] = p
        if nx > 1:
            varying['x'] = x

        # If all three coordinates vary, use structural analysis
        if 't' in varying and 'p' in varying and 'x' in varying:
            # Check if X is structured as independent variable with PT as path
            x_is_independent = is_independent_variable(x, [p, t])

            if x_is_independent:
                print(f"Auto-detected X-path configuration: X is independent variable with varying P-T conditions")
                return x, p  # X on x-axis, P on y-axis
            else:
                print(f"All coordinates vary but X is not independent variable. Using standard T-P configuration.")
                return t, p

        # Standard priority logic for cases with 2 or fewer varying coordinates
        if 't' in varying and 'p' in varying:
            return t, p
        if 'x' in varying and 'p' in varying:
            return x, p
        if 'x' in varying and 't' in varying:
            return x, t

        # Fallback: first two varying coordinates or first two columns
        varying_list = [varying[k] for k in ['t', 'p', 'x'] if k in varying]
        if len(varying_list) >= 2:
            return varying_list[0], varying_list[1]

        if len(self.coord_columns) >= 2:
            return self.coord_columns[0], self.coord_columns[1]

        return None, None

    def set_plot_coordinates(self, xname: str, yname: str):
        """Set coordinate columns for plotting and build x, y arrays."""
        self.coords_x = xname
        self.coords_y = yname

        # Extract x, y values from blocks
        x_vals = []
        y_vals = []

        if self.model_blocks:
            available_coords = list(self.model_blocks[0].get("coords", {}).keys())

            x_found = xname in available_coords
            y_found = yname in available_coords

            if not x_found or not y_found:
                raise ValueError(
                    f"Requested coordinates not found in data. "
                    f"Requested: x='{xname}', y='{yname}'. "
                    f"Available: {available_coords}"
                )

        for block in self.model_blocks:
            coords = block.get("coords", {})
            if xname in coords and yname in coords:
                x_vals.append(coords[xname])
                y_vals.append(coords[yname])

        self.x = np.array(x_vals)
        self.y = np.array(y_vals)

        # Initialize labels with original coordinate names
        self.labels = [xname, yname]

        # Helper function to detect and convert coordinates
        def convert_coordinate(values, name):
            name_lower = name.lower()
            
            # Temperature: T(K) or T[K] variants, convert to Celsius
            if "t(k)" in name_lower or "t[k]" in name_lower or (("temp" in name_lower or name_lower == "t") and "k" in name_lower):
                return values - 273, "T(°C)"
            
            # Temperature: Already in Celsius but with bracket notation
            if "t[°c]" in name_lower or "t[c]" in name_lower:
                return values, "T(°C)"
            
            # Pressure in bar: P(bar) or variants, excluding kbar
            if "p(bar)" in name_lower or (("p" in name_lower or "pressure" in name_lower) and "bar" in name_lower and "kbar" not in name_lower):
                return values / 10000, "P(GPa)"
            
            # Pressure in kbar: P[kbar] or variants containing kbar
            if "p[kbar]" in name_lower or "kbar" in name_lower:
                return values / 10, "P(GPa)"
            
            # No conversion needed
            return values, name

        # Apply conversions
        self.x, self.labels[0] = convert_coordinate(self.x, xname)
        self.y, self.labels[1] = convert_coordinate(self.y, yname)

        self._log_print("The coordinate variables are:  " + str(self.labels))
        self._log_print(f"The {self.labels[0]} range is: {self.x.min():.4g} to {self.x.max():.4g}")
        self._log_print(f"The {self.labels[1]} range is: {self.y.min():.4g} to {self.y.max():.4g}")

    def import_analytical_compo(self, filename: Optional[str] = None):
        """Import composition data from a file."""
        if filename is None:
            self._log_print("Please select a text file containing measured compositions in apfu.")
            filename = filedialog.askopenfilename()
            if not filename:
                return
        
        self._log_print("The composition file is:  " + filename)
        
        input_data = pd.read_csv(filename, sep='\t', header=None, comment='#')
        arrays = input_data.values
        
        self.apfu_name = [str(i) for i in arrays[0] if str(i) != 'nan']
        self.apfu_obs = np.array([float(i) for i in arrays[1] if str(i) != 'nan'])
        
        # Handle observation errors
        obs_err_raw = np.array([float(i) if str(i) != '-' else '-' 
                               for i in arrays[2] if str(i) != 'nan'])
        
        if np.all(obs_err_raw == '-'):
            self.obs_err = np.array([])
        else:
            self.obs_err = np.where((obs_err_raw != '-') & (obs_err_raw.astype(float) < 0.01), 
                                   0.01, obs_err_raw).astype(float)
        
        self.analysis_type = ''.join([str(i) for i in arrays[3] if str(i) != 'nan'])
        
        # Extract phase names from Phase_Element format and maintain order of first appearance
        phase_name_ordered = []
        for apfu_full_name in self.apfu_name:
            if "_" in apfu_full_name:
                phase = apfu_full_name.split("_", 1)[0]
                if phase not in phase_name_ordered:
                    phase_name_ordered.append(phase)
        self.phase_name = phase_name_ordered
        self.color_scheme = ''.join([str(i) for i in arrays[5] if str(i) != 'nan'])
        
        self._log_print("The input variables are: " + str(self.apfu_name))
        self._log_print("The input compositions are: " + str(self.apfu_obs))
        self._log_print("The input uncertainties are: " + str(self.obs_err))
        self._log_print("The selected analysis type is: " + self.analysis_type)
        self._log_print("The phase names are: " + str(self.phase_name))
        self._log_print("The color scheme is: " + self.color_scheme)
        
        if len(self.obs_err) == 0:
            self.obs_err = self.calc_obs_err(self.apfu_obs)

    # ======================
    # Build intersect table 
    # ======================
    
    def build_intersect_table(self):
        """
        Build phase_data DataFrame and legacy self.data structure.
        """
        if not self.model_blocks:
            raise RuntimeError("No model blocks loaded")
        if not (self.coords_x and self.coords_y):
            raise RuntimeError("Plot coordinates not set")
        if len(self.apfu_name) == 0:
            raise RuntimeError("No observed elements")

        # Parse observed phase-element combinations
        observed_phases_elements = {}
        for apfu_full_name in self.apfu_name:
            if "_" not in apfu_full_name:
                continue
            parts = apfu_full_name.split("_", 1)
            phase_warr = parts[0]
            element = parts[1]

            if phase_warr not in observed_phases_elements:
                observed_phases_elements[phase_warr] = []
            observed_phases_elements[phase_warr].append(element)

        records = []

        for block in self.model_blocks:
            coords = block.get("coords", {})

            try:
                x_val = coords[self.coords_x]
                y_val = coords[self.coords_y]
            except KeyError:
                continue
            
            # Apply unit conversions
            if self.labels[0] == "T(°C)" and self.coords_x == "T(K)":
                x_val -= 273
            elif self.labels[0] == "P(GPa)" and self.coords_x == "P(bar)":
                x_val /= 10000

            if self.labels[1] == "P(GPa)" and self.coords_y == "P(bar)":
                y_val /= 10000
            elif self.labels[1] == "T(°C)" and self.coords_y == "T(K)":
                y_val -= 273

            # Use already processed phases from block
            phases = block.get("phases", [])

            # Build row
            row = {"x": float(x_val), "y": float(y_val)}

            for phase_warr, elements in observed_phases_elements.items():
                matching_phases = [p for p in phases 
                                 if p.get("phase", "").lower() == phase_warr.lower()]

                if not matching_phases:
                    for element in elements:
                        column_name = f"{phase_warr}_{element}"
                        row[column_name] = np.nan
                    continue
                
                model_phase = matching_phases[0]

                for element in elements:
                    value = self._find_element_in_phase(model_phase, element)
                    column_name = f"{phase_warr}_{element}"
                    row[column_name] = float(value) if value is not None else np.nan

            records.append(row)

        if not records:
            raise RuntimeError("No valid records assembled")

        self.phase_data = pd.DataFrame.from_records(records)
        self.phase_data.sort_values(["y", "x"], inplace=True, ignore_index=True)
        self.x = self.phase_data["x"].values
        self.y = self.phase_data["y"].values

        # Build legacy self.data structure for compatibility with original methods
        self._build_legacy_data_structure()
    
    def _find_element_in_phase(self, phase_dict: dict, target_element: str):
        """Find element value in phase dictionary."""
        target_norm = normalize_element_name(target_element)
        
        # Try simple element key first (created by parser for discrimination)
        if target_norm in phase_dict:
            return phase_dict[target_norm]
        
        # Try _apfu key (normalized by parser)
        key_apfu = f"{target_norm}_apfu"
        if key_apfu in phase_dict:
            return phase_dict[key_apfu]
        
        # Fallback: search through all keys and extract with coefficient
        for key, value in phase_dict.items():
            if not isinstance(key, str) or key in ("phase", "orig_phase"):
                continue
            
            elem, coef = extract_element_info(key)
            elem_norm = normalize_element_name(elem)
            
            if elem_norm == target_norm:
                try:
                    num_val = float(value)
                    # Only apply coefficient if this is a raw oxide column
                    # (not already an _apfu column)
                    if not key.endswith("_apfu"):
                        return num_val * coef
                    else:
                        return num_val
                except (ValueError, TypeError):
                    continue
                
        return None
    
    def _build_legacy_data_structure(self):
        """Build legacy self.data structure for compatibility with calculation methods."""
        # Create column labels matching original format
        data_cols = [c for c in self.phase_data.columns if c not in ("x", "y")]
        self.labels = np.array([self.labels[0], self.labels[1]] + data_cols)
        
        # Build data DataFrame
        data_dict = {
            self.labels[0]: self.x,
            self.labels[1]: self.y
        }
        for col in data_cols:
            data_dict[col] = self.phase_data[col].values
        
        self.data = pd.DataFrame(data_dict)
        
        # Build phase_id array (group by phase name)
        phase_ids = {}
        current_id = 1
        self.phase_id = []
        
        for col in data_cols:
            if "_" in col:
                phase = col.split("_")[0]
                if phase not in phase_ids:
                    phase_ids[phase] = current_id
                    current_id += 1
                self.phase_id.append(phase_ids[phase])
            else:
                self.phase_id.append(0)
        
        self.phase_id = np.array(self.phase_id)

    # ====================
    # Error calculation
    # ====================
    
    def calc_obs_err(self, apfu_obs):
        """Calculate observation errors based on analysis type."""
        if self.analysis_type == 'EDS':
            calc_err = 0.0703 * (apfu_obs**0.3574)
            min_err, max_err = 0.01, 0.1
        elif self.analysis_type == 'WDS map':
            calc_err = 0.0434 * (apfu_obs**0.3451)
            min_err, max_err = 0.005, 0.05
        elif self.analysis_type == 'WDS spot':
            calc_err = 0.023 * (apfu_obs**0.2772)
            min_err, max_err = 0.005, 0.05
        else:
            self._log_print('Please enter a valid analysis type (EDS, WDS map, WDS spot)')
            sys.exit()
        
        calc_err = np.clip(calc_err, min_err, max_err)
        self._log_print('The calculated uncertainties are:  ' + str(calc_err))
        return calc_err

    # ================================
    # Quality factor calculations
    # ================================
    
    def Q_elem(self, apfu_obs, model, obs_err):
        """Calculate quality factor for each element."""
        factor_min = 1
        factor_max = 6
        obs_err = np.where(obs_err < 0.01, 0.01, obs_err)
        diff = np.abs(apfu_obs - model)
        num = np.clip(diff - obs_err / factor_min, 0, factor_max * obs_err)
        Qcmp_elem = 100 * np.abs(1 - num / (factor_max * obs_err))**(model + 1)
        return Qcmp_elem
    
    def Q_phase(self, apfu_obs, model, obs_err):
        """Calculate quality factor for each phase."""
        factor_min = 1
        factor_max = 6
        obs_err = np.where(obs_err < 0.01, 0.01, obs_err)
        diff = np.abs(apfu_obs - model)
        num = np.clip(diff - obs_err / factor_min, 0, factor_max * obs_err)
        Qcmp_elem = np.abs(1 - num / (factor_max * obs_err))**(model + 1)
        Qcmp_phase = np.sum(Qcmp_elem) / np.size(Qcmp_elem) * 100
        return Qcmp_phase
    
    def chi2(self, apfu_obs, model, obs_err):
        """Calculate chi-squared statistics."""
        return np.sum((apfu_obs - model)**2 / obs_err**2)
    
    def red_chi2(self, apfu_obs, model, obs_err, f):
        """Calculate reduced chi-squared statistics."""
        return self.chi2(apfu_obs, model, obs_err) / (f - 1)
    
    def norm_weight(self, weight):
        """Normalize the weight array."""
        weight_norm = weight / np.sum(weight)
        return weight_norm
    
    def Q_tot(self, Qcmp_tot, weight_norm):
        """Calculate total quality factor."""
        return np.sum(Qcmp_tot * weight_norm)

    # =====================
    # Plotting methods
    # =====================
    
    def plot_elem(self, Qcmp, i):
        """Plot quality factor for each element."""
        Qcmp_2D = np.reshape(Qcmp, (len(np.unique(self.y)), len(np.unique(self.x))))
        plt.imshow(Qcmp_2D, cmap=self.color_scheme, aspect='auto', origin='lower', 
                  extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.colorbar()
        plt.title('Quality factor for ' + self.apfu_name[i])
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        plt.clim(0, 100)
        
        self._log_print('The maximum value of the quality factor for ' + self.apfu_name[i] + ' is: ' + str(np.nanmax(Qcmp_2D)))
        
        contoured = plt.contour(Qcmp_2D, levels=np.arange(0, 100, 10), colors="white", 
                              linewidths=0.5, origin="lower", 
                              extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.0f')
        
        output_path = os.path.join(self.output_dir)
        os.makedirs(output_path, exist_ok=True)
        plt.savefig(os.path.join(output_path, "Qcmp_" + self.apfu_name[i] + '.pdf'), format='pdf')
        plt.show()
        plt.close()
    
    def plot_phase(self, Qcmp, i):
        """Plot quality factor for each phase."""
        Qcmp_2D = np.reshape(Qcmp, (len(np.unique(self.y)), len(np.unique(self.x))))
        plt.imshow(Qcmp_2D, cmap=self.color_scheme, aspect='auto', origin='lower', 
                  extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.colorbar()
        plt.title('Quality factor for ' + self.phase_name[i-1])
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        plt.clim(0, 100)
        
        self._log_print('The maximum value of the quality factor for ' + self.phase_name[i-1] + ' is: ' + str(np.nanmax(Qcmp_2D)))
        
        contoured = plt.contour(Qcmp_2D, levels=np.arange(0, 100, 10), colors="white", 
                              linewidths=0.5, origin="lower", 
                              extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.0f')
        
        if np.nanmax(Qcmp_2D) < 100:
            n_p = np.count_nonzero(self.y == self.y[0])
            max_Qcmp = np.where(Qcmp_2D == np.nanmax(Qcmp_2D))
            max_Qcmp_y = self.y[max_Qcmp[0]*n_p]
            max_Qcmp_x = self.x[max_Qcmp[1]]
            self._log_print('The ' + self.labels[0] + ' and ' + self.labels[1] + ' position of the maximum Qcmp of ' + self.phase_name[i-1] + ' is: ' + str(max_Qcmp_x) + ' , ' + str(max_Qcmp_y))
        
        output_path = os.path.join(self.output_dir)
        os.makedirs(output_path, exist_ok=True)
        plt.savefig(os.path.join(output_path, "Qcmp_" + self.phase_name[i-1] + '.pdf'), format='pdf')
        plt.show()
        plt.close()
    
    def plot_tot(self, Qcmp, title):
        """Plot total quality factor."""
        Qcmp_2D = np.reshape(Qcmp, (len(np.unique(self.y)), len(np.unique(self.x))))
        plt.imshow(Qcmp_2D, cmap=self.color_scheme, aspect='auto', origin='lower', 
                  extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.colorbar()
        plt.title(title + ' Q*cmp')
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        plt.clim(0, 100)
        
        self._log_print('The maximum value of the Q*cmp is: ' + str(np.nanmax(Qcmp_2D)))
        
        contoured = plt.contour(Qcmp_2D, levels=np.arange(0, 110, 10), colors="white", 
                              linewidths=0.5, origin="lower", 
                              extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.0f')
        
        max_Qcmp = np.where(Qcmp_2D == np.nanmax(Qcmp_2D))
        n_p = np.count_nonzero(self.y == self.y[0])
        max_Qcmp_y = self.y[max_Qcmp[0]*n_p]
        max_Qcmp_x = self.x[max_Qcmp[1]]
        max_Qcmp_x = np.mean(max_Qcmp_x)
        max_Qcmp_y = np.mean(max_Qcmp_y)
        
        self._log_print('The ' + self.labels[0] + ' and ' + self.labels[1] + ' position of the maximum Q*cmp is: ' + str(max_Qcmp_x) + ' , ' + str(max_Qcmp_y))
        
        plt.plot(max_Qcmp_x, max_Qcmp_y, "ro", markersize=2)
        os.makedirs(self.output_dir, exist_ok=True)
        plt.savefig(os.path.join(self.output_dir, title + "_Qcmp_tot.pdf"), format='pdf')
        plt.show()
        plt.close()
    
    def plot_redchi2_phase(self, redchi2, i, f):
        """Plot reduced chi-squared for phase."""
        redchi2_2D = np.reshape(redchi2, (len(np.unique(self.y)), len(np.unique(self.x))))
        
        if f > 2:
            self._log_print("The number of elements in " + self.phase_name[i-1] + " is: " + str(f))
            plt.imshow(redchi2_2D, cmap=self.color_scheme, aspect='auto', origin='lower', 
                      extent=[min(self.x), max(self.x), min(self.y), max(self.y)], 
                      norm=colors.LogNorm())
            plt.colorbar()
            plt.title('Reduced χ2 ' + self.phase_name[i-1])
            plt.xlabel(self.labels[0])
            plt.ylabel(self.labels[1])
            
            min_redchi2 = np.nanmin(redchi2_2D)
            max_redchi2 = np.nanmax(redchi2_2D)
            step = (max_redchi2 - min_redchi2) / 10
            
            self._log_print('The minimum reduced χ2 value for ' + self.phase_name[i-1] + ' is: ' + str(min_redchi2))
            
            if min_redchi2 <= 1:
                contoured = plt.contour(redchi2_2D, levels=np.arange(1, max_redchi2, step), 
                                      colors="white", linewidths=0.5, origin="lower", 
                                      extent=[min(self.x), max(self.x), min(self.y), max(self.y)], 
                                      norm=colors.LogNorm())
            else:
                contoured = plt.contour(redchi2_2D, 
                                      levels=np.arange(np.round(min_redchi2, decimals=1)+0.1, max_redchi2, step), 
                                      colors="white", linewidths=0.5, origin="lower", 
                                      extent=[min(self.x), max(self.x), min(self.y), max(self.y)], 
                                      norm=colors.LogNorm())
            plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.1f')
            
            os.makedirs(self.output_dir, exist_ok=True)
            plt.savefig(os.path.join(self.output_dir, "redχ2_" + self.phase_name[i-1] + '.pdf'), format='pdf')
            plt.show()
            plt.close()
        else:
            self._log_print("The number of elements in " + self.phase_name[i-1] + " is: " + str(f))
            plt.imshow(redchi2_2D, cmap=self.color_scheme, aspect='auto', origin='lower', 
                      extent=[min(self.x), max(self.x), min(self.y), max(self.y)], 
                      norm=colors.LogNorm())
            plt.colorbar()
            plt.title('χ2 ' + self.phase_name[i-1])
            plt.xlabel(self.labels[0])
            plt.ylabel(self.labels[1])
            
            min_redchi2 = np.nanmin(redchi2_2D)
            max_redchi2 = np.nanmax(redchi2_2D)
            step = (max_redchi2 - min_redchi2) / 10
            
            self._log_print('The minimum χ2 value for ' + self.phase_name[i-1] + ' is: ' + str(min_redchi2))
            
            if min_redchi2 <= 1 and min_redchi2 != max_redchi2:
                contoured = plt.contour(redchi2_2D, levels=np.arange(1, max_redchi2, step), 
                                      colors="white", linewidths=0.5, origin="lower", 
                                      extent=[min(self.x), max(self.x), min(self.y), max(self.y)], 
                                      norm=colors.LogNorm())
                plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.1f')
            elif min_redchi2 > 1 and min_redchi2 != max_redchi2:
                contoured = plt.contour(redchi2_2D, 
                                      levels=np.arange(np.round(min_redchi2, decimals=1)+0.1, max_redchi2, step), 
                                      colors="white", linewidths=0.5, origin="lower", 
                                      extent=[min(self.x), max(self.x), min(self.y), max(self.y)], 
                                      norm=colors.LogNorm())
                plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.1f')
            
            os.makedirs(self.output_dir, exist_ok=True)
            plt.savefig(os.path.join(self.output_dir, "χ2_" + self.phase_name[i-1] + '.pdf'), format='pdf')
            plt.show()
            plt.close()
        
        return min_redchi2
    
    def plot_redchi2_tot(self, redchi2):
        """Plot total reduced chi-squared."""
        redchi2_2D = np.reshape(redchi2, (len(np.unique(self.y)), len(np.unique(self.x))))
        plt.imshow(redchi2_2D, cmap=self.color_scheme, aspect='auto', origin='lower', 
                  extent=[min(self.x), max(self.x), min(self.y), max(self.y)], 
                  norm=colors.LogNorm())
        plt.colorbar()
        plt.title('Total reduced χ2')
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        
        min_redchi2 = np.nanmin(redchi2_2D)
        max_redchi2 = np.nanmax(redchi2_2D)
        step = (max_redchi2 - min_redchi2) / 10
        
        if min_redchi2 <= 1:
            contoured = plt.contour(redchi2_2D, levels=np.arange(1, max_redchi2, step), 
                                  colors="white", linewidths=0.5, origin="lower", 
                                  extent=[min(self.x), max(self.x), min(self.y), max(self.y)], 
                                  norm=colors.LogNorm())
        else:
            contoured = plt.contour(redchi2_2D, 
                                  levels=np.arange(np.round(min_redchi2, decimals=1)+0.1, max_redchi2, step), 
                                  colors="white", linewidths=0.5, origin="lower", 
                                  extent=[min(self.x), max(self.x), min(self.y), max(self.y)], 
                                  norm=colors.LogNorm())
        plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.1f')
        
        self._log_print('The minimum value of the total reduced χ2 is: ' + str(min_redchi2))
        
        n_p = np.count_nonzero(self.y == self.y[0])
        min_redchi2_pos = np.where(redchi2_2D == np.nanmin(redchi2_2D))
        min_redchi2_y = self.y[min_redchi2_pos[0]*n_p]
        min_redchi2_x = self.x[min_redchi2_pos[1]]
        
        self._log_print('The ' + self.labels[0] + ' and ' + self.labels[1] + ' position of the minimum total reduced χ2 is: ' + str(min_redchi2_x) + ' , ' + str(min_redchi2_y))
        
        os.makedirs(self.output_dir, exist_ok=True)
        plt.savefig(os.path.join(self.output_dir, "redχ2_tot.pdf"), format='pdf')
        plt.show()
        plt.close()
        
        return min_redchi2

    # =========================
    # Main calculation methods 
    # =========================
    
    def Qcmp_elem(self):
        """Calculate and plot quality factor for each element."""
        Qcmp_elem = np.empty(len(self.data[self.labels[2]]))
        
        for i in range(len(self.labels[2:])):
            Model_elem = np.array(self.data[self.labels[i+2]])
            for j in range(len(Model_elem)):
                Qcmp_elem[j] = self.Q_elem(self.apfu_obs[i], Model_elem[j], self.obs_err[i])
            self.plot_elem(Qcmp_elem, i)
    
    def Qcmp_phase(self):
        """Calculate and plot quality factor for each phase."""
        phase_idx = []
        Qcmp_phase_tot = np.empty((len(self.y), len(np.unique(self.phase_id))))
        
        for i in np.unique(self.phase_id):
            idx = np.where(self.phase_id == i)[0]+2
            phase_idx.append(idx)
            
            Model_phase = np.array(self.data[self.labels[idx]])
            apfu_obs_idx = self.apfu_obs[idx-2]
            obs_err_idx = self.obs_err[idx-2]
            idx_phase = np.arange(len(idx))
            Qcmp_phase = np.empty(len(Model_phase))
            
            for j in range(len(Model_phase)):
                Qcmp_phase[j] = self.Q_phase(apfu_obs_idx, Model_phase[j, idx_phase], obs_err_idx)
            
            self.plot_phase(Qcmp_phase, i)
            Qcmp_phase_tot[:, i-1] = Qcmp_phase
        
        return Qcmp_phase_tot
    
    def redchi2_phase(self):
        """Calculate and plot reduced chi-squared for each phase."""
        phase_idx = []
        redchi2_phase_tot = np.empty((len(self.y), len(np.unique(self.phase_id))))
        
        for i in np.unique(self.phase_id):
            idx = np.where(self.phase_id == i)[0]+2
            phase_idx.append(idx)
            f = len(idx)
            
            Model_phase = np.array(self.data[self.labels[idx]])
            apfu_obs_idx = self.apfu_obs[idx-2]
            obs_err_idx = self.obs_err[idx-2]
            idx_phase = np.arange(len(idx))
            redchi2_phase = np.empty(len(Model_phase))
            
            if f > 2:
                for j in range(len(Model_phase)):
                    redchi2_phase[j] = self.red_chi2(apfu_obs_idx, Model_phase[j, idx_phase], obs_err_idx, f)
                min_redchi2_phase = self.plot_redchi2_phase(redchi2_phase, i, f)
            else:
                for j in range(len(Model_phase)):
                    redchi2_phase[j] = self.chi2(apfu_obs_idx, Model_phase[j, idx_phase], obs_err_idx)
                min_redchi2_phase = self.plot_redchi2_phase(redchi2_phase, i, f)
                min_redchi2_phase = 1 + min_redchi2_phase
                if f == 1:
                    self._log_print('WARNING: The number of elements in ' + self.phase_name[i-1] + ' is 1, the χ2 may not be accurate')
            
            redchi2_phase_tot[:, i-1] = redchi2_phase
            self.min_redchi2 = np.append(self.min_redchi2, min_redchi2_phase)
        
        return redchi2_phase_tot
    
    def redchi2_tot(self):
        """Calculate and plot total reduced chi-squared."""
        f = len(self.labels) - 2
        redchi2_tot = np.empty(len(self.data[self.labels[2]]))
        
        for i in range(len(self.data[self.labels[2]])):
            Model_tot = np.array(self.data[self.labels[2:]])
            apfu_obs_tot = self.apfu_obs
            obs_err_tot = self.obs_err
            redchi2_tot[i] = self.red_chi2(apfu_obs_tot, Model_tot[i], obs_err_tot, f)
        
        min_redchi2_tot = self.plot_redchi2_tot(redchi2_tot)
        return min_redchi2_tot
    
    def Qcmp_tot(self, Qcmp_phase_tot, redchi2_phase_tot):
        """Calculate total quality factor (unweighted)."""
        if self.phase_id is None:
            raise ValueError("phase_id is not set.")
        
        weight = np.ones(len(np.unique(self.phase_id)))
        weight_norm = self.norm_weight(weight)
        
        Qcmp_allphases = np.empty(len(Qcmp_phase_tot[:, 0]))
        for i in range(len(Qcmp_phase_tot[:, 0])):
            Qcmp_allphases[i] = self.Q_tot(Qcmp_phase_tot[i, :], weight_norm)
        
        max_Qcmp = np.where(Qcmp_allphases == np.nanmax(Qcmp_allphases))
        Qcmpmax_redchi2_value = redchi2_phase_tot[max_Qcmp[0][0]]
        
        self._log_print("The reduced χ2 values for the phases at the maximum Q*cmp is: " + str(Qcmpmax_redchi2_value))
        
        self.plot_tot(Qcmp_allphases, "Unweighted")
    
    def Qcmp_tot_weight(self, Qcmp_phase_tot, redchi2_phase_tot):
        """Calculate total quality factor (weighted)."""
        self.min_redchi2[self.min_redchi2 < 1] = 1
        weight = 1 / self.min_redchi2
        
        self._log_print('The weight is:  ' + str(weight))
        
        weight_norm = self.norm_weight(weight)
        self._log_print('The normalized weight fraction is:  ' + str(weight_norm))
        
        Qcmp_allphases_weight = np.empty(len(Qcmp_phase_tot))
        for i in range(len(Qcmp_phase_tot)):
            Qcmp_allphases_weight[i] = self.Q_tot(Qcmp_phase_tot[i, :], weight_norm)
        
        max_Qcmp = np.where(Qcmp_allphases_weight == np.nanmax(Qcmp_allphases_weight))
        Qcmpmax_redchi2_value = redchi2_phase_tot[max_Qcmp[0][0]]
        
        self._log_print("The reduced χ2 values for the phases at the maximum Q*cmp are: " + str(Qcmpmax_redchi2_value))
        
        self.plot_tot(Qcmp_allphases_weight, "Weighted")
        return Qcmp_allphases_weight


if __name__ == "__main__":
    QFA = QualityFactorAnalysis()
    QFA.run_analysis(0)