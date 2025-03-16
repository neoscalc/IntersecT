#%%
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from tkinter import filedialog

class QualityFactorAnalysis:
    def __init__(self):
        self.data = pd.DataFrame()
        self.labels = []
        self.x = []
        self.y = []
        self.labels_new = []
        self.phase_id = []
        self.apfu_name = []
        self.apfu_obs = []
        self.obs_err = []
        self.analysis_type = ''
        self.color_scheme = 'viridis'
        self.phase_name = []
        self.output_dir = ''
        self.min_redchi2 = []

    def set_output_directory(self):
        """Prompt user to select an output directory."""
        print("Please select an output directory.")

        self.output_dir = filedialog.askdirectory()
        print("The output directory is: ", self.output_dir)

    def import_output_perplex(self):
        """Import data from a perplex tab file and process it."""
        print("Please select a tab file from perplex.")

        filename = filedialog.askopenfilename()
        self.data = pd.read_csv(filename, sep=r'\s+', header=0, skiprows=12)
        print("The data file is: ", filename)

        self.labels = np.array(self.data.columns)  # Convert Index to a list
        self.x, self.y = np.array(self.data[self.labels[0]]), np.array(self.data[self.labels[1]])

        # Convert units if necessary
        if self.labels[0] == "T(K)":
            self.x -= 273
            self.labels[0] = "T(°C)"
        elif self.labels[0] == "P(bar)":
            self.x /= 10000
            self.labels[0] = "P(GPa)"

        if self.labels[1] == "P(bar)":
            self.y /= 10000
            self.labels[1] = "P(GPa)"
        elif self.labels[1] == "T(K)":
            self.y -= 273
            self.labels[1] = "T(°C)"

        print("The werami variables are: ", self.labels)
        print("The", self.labels[0], "range is:", self.x)
        print("The", self.labels[1], "range is:", self.y)

        self.labels_new = np.array([''.join([j for j in label if not j.isdigit()]) for label in self.labels[2:]])

        id = 0
        self.phase_id = np.array([id if label == self.labels_new[i-1] else (id := id + 1) for i, label in enumerate(self.labels_new)])


    def import_analytical_compo(self):
        """Import composition data from a file."""

        print("Please select a text file containing analyses.")

        filename = filedialog.askopenfilename()
        input_data = pd.read_csv(filename, sep='\t', header=None, comment='#')
        print("The composition file is: ", filename)

        arrays = input_data.values

        self.apfu_name = np.array([str(i) for i in arrays[0] if str(i) != 'nan'])
        self.apfu_obs = np.array([float(i) for i in arrays[1] if str(i) != 'nan'])
        self.obs_err = np.array([float(i) if str(i) != '-' else '-' for i in arrays[2] if str(i) != 'nan'])

        # If self.obs_err contains only '-', set it to an empty array
        if np.all(self.obs_err == '-'):
            self.obs_err = np.array([])  # Set to empty array if all values are '-'
        else:
            # Replace valid numeric entries below 0.01 with 0.01
            self.obs_err = np.where((self.obs_err != '-') & (self.obs_err < 0.01), 0.01, self.obs_err)

        self.analysis_type = ''.join([str(i) for i in arrays[3] if str(i) != 'nan'])
        self.phase_name = np.array([str(i) for i in arrays[5] if str(i) != 'nan'])
        self.color_scheme = ''.join([str(i) for i in arrays[7] if str(i) != 'nan'])

        print("The input variables are:", self.apfu_name)
        print("The input compositions are:", self.apfu_obs)
        print("The input uncertainties are:", self.obs_err)
        print("The selected analysis type is:", self.analysis_type)
        print("The phase names are:", self.phase_name)
        print("The color scheme is:", self.color_scheme)

        if len(self.obs_err) == 0:
            self.obs_err = self.calc_obs_err(self.apfu_obs)

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
            print('Please enter a valid analysis type (EDS, WDS map, WDS spot)')
            sys.exit()

        calc_err = np.clip(calc_err, min_err, max_err)
        print('The calculated uncertainties are: ', calc_err)
        return calc_err

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
        """Calculate chi-squared value."""
        return np.sum((apfu_obs - model)**2 / obs_err**2)

    def red_chi2(self, apfu_obs, model, obs_err, f):
        """Calculate reduced chi-squared value."""
        return self.chi2(apfu_obs, model, obs_err) / (f - 1)

    def norm_weight(self, weight):
        """Normalize the weight array by dividing each element by the sum of the weights."""
        weight_norm = weight / np.sum(weight)
        return weight_norm

    def Q_tot(self, Qcmp_tot, weight_norm):
        """Calculate total quality factor."""
        return np.sum(Qcmp_tot * weight_norm)

    def plot_elem(self, Qcmp, i):
        """Plot quality factor for each element."""
        Qcmp_2D = np.reshape(Qcmp, (len(np.unique(self.y)), len(np.unique(self.x))))
        plt.imshow(Qcmp_2D, cmap=self.color_scheme, aspect='auto', origin='lower', extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.colorbar()
        plt.title('Quality factor for ' + self.apfu_name[i])
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        plt.clim(0, 100)

        print('The maximum value of the quality factor for', self.apfu_name[i], 'is:', np.nanmax(Qcmp_2D))
        contoured = plt.contour(Qcmp_2D, levels=np.arange(0, 100, 10), colors="white", linewidths=0.5, origin="lower", extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.0f')

        # Save plot as PDF
        output_path = os.path.join(self.output_dir)
        os.makedirs(output_path, exist_ok=True)
        plt.savefig(os.path.join(output_path, "Qcmp_" + self.apfu_name[i] + '.pdf'), format='pdf')

        plt.show()
        plt.close()

    def plot_phase(self, Qcmp, i):
        """Plot quality factor for each phase."""
        Qcmp_2D = np.reshape(Qcmp, (len(np.unique(self.y)), len(np.unique(self.x))))
        plt.imshow(Qcmp_2D, cmap=self.color_scheme, aspect='auto', origin='lower', extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.colorbar()
        plt.title('Quality factor for ' + self.phase_name[i-1])
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        plt.clim(0, 100)

        print('The maximum value of the quality factor for', self.phase_name[i-1], 'is:', np.nanmax(Qcmp_2D))
        contoured = plt.contour(Qcmp_2D, levels=np.arange(0, 100, 10), colors="white", linewidths=0.5, origin="lower", extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.0f')

        if np.nanmax(Qcmp_2D) < 100:
            n_p = np.count_nonzero(self.y == self.y[0])
            max_Qcmp = np.where(Qcmp_2D == np.nanmax(Qcmp_2D))
            max_Qcmp_y = self.y[max_Qcmp[0]*n_p]
            max_Qcmp_x = self.x[max_Qcmp[1]]
            print('The ' + self.labels[0] + ' and ' + self.labels[1] + ' position of the maximum Qcmp of', self.phase_name[i-1], 'is:', max_Qcmp_x, ',', max_Qcmp_y)

        # Save plot as PDF
        output_path = os.path.join(self.output_dir)
        os.makedirs(output_path, exist_ok=True)
        plt.savefig(os.path.join(output_path, "Qcmp_" + self.phase_name[i-1] + '.pdf'), format='pdf')

        plt.show()
        plt.close()

    def plot_tot(self, Qcmp, title):
        Qcmp_2D = np.reshape(Qcmp, (len(np.unique(self.y)), len(np.unique(self.x))))
        plt.imshow(Qcmp_2D, cmap=self.color_scheme, aspect='auto', origin='lower', extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.colorbar()
        plt.title(title + ' Q*cmp')
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        plt.clim(0, 100)
        print('The maximum value of the Q*cmp is:', np.nanmax(Qcmp_2D))

        contoured = plt.contour(Qcmp_2D, levels=np.arange(0, 110, 10), colors="white", linewidths=0.5, origin="lower", extent=[min(self.x), max(self.x), min(self.y), max(self.y)])
        plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.0f')

        max_Qcmp = np.where(Qcmp_2D == np.nanmax(Qcmp_2D))
        n_p = np.count_nonzero(self.y == self.y[0])
        max_Qcmp_y = self.y[max_Qcmp[0]*n_p]
        max_Qcmp_x = self.x[max_Qcmp[1]]
        max_Qcmp_x = np.mean(max_Qcmp_x)
        max_Qcmp_y = np.mean(max_Qcmp_y)
        print('The ' + self.labels[0] + ' and ' + self.labels[1] + ' position of the maximum Q*cmp is:', max_Qcmp_x, ',', max_Qcmp_y)

        plt.plot(max_Qcmp_x, max_Qcmp_y, "ro", markersize=2)
        os.makedirs(self.output_dir, exist_ok=True)
        plt.savefig(os.path.join(self.output_dir, title + "_Qcmp_tot.pdf"), format='pdf')
        plt.show()
        plt.close()

    def plot_redchi2_phase(self, redchi2, i, f):
        redchi2_2D = np.reshape(redchi2, (len(np.unique(self.y)), len(np.unique(self.x))))
        if f > 2:
            print("The number of elements in", self.phase_name[i-1], "is:", f)
            plt.imshow(redchi2_2D, cmap=self.color_scheme, aspect='auto', origin='lower', extent=[min(self.x), max(self.x), min(self.y), max(self.y)], norm=colors.LogNorm())
            plt.colorbar()
            plt.title('Reduced χ2 ' + self.phase_name[i-1])
            plt.xlabel(self.labels[0])
            plt.ylabel(self.labels[1])

            min_redchi2 = np.nanmin(redchi2_2D)
            max_redchi2 = np.nanmax(redchi2_2D)
            step = (max_redchi2 - min_redchi2) / 10
            print('The minimum reduced χ2 value for', self.phase_name[i-1], 'is:', min_redchi2)

            if min_redchi2 <= 1:
                contoured = plt.contour(redchi2_2D, levels=np.arange(1, max_redchi2, step), colors="white", linewidths=0.5, origin="lower", extent=[min(self.x), max(self.x), min(self.y), max(self.y)], norm=colors.LogNorm())
            else:
                contoured = plt.contour(redchi2_2D, levels=np.arange(np.round(min_redchi2, decimals=1)+0.1, max_redchi2, step), colors="white", linewidths=0.5, origin="lower", extent=[min(self.x), max(self.x), min(self.y), max(self.y)], norm=colors.LogNorm())

            plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.1f')
            os.makedirs(self.output_dir, exist_ok=True)
            plt.savefig(os.path.join(self.output_dir, "redχ2_" + self.phase_name[i-1] + '.pdf'), format='pdf')
            plt.show()
            plt.close()
        
        else:
            print("The number of elements in", self.phase_name[i-1], "is:", f)
            plt.imshow(redchi2_2D, cmap=self.color_scheme, aspect='auto', origin='lower', extent=[min(self.x), max(self.x), min(self.y), max(self.y)], norm=colors.LogNorm())
            plt.colorbar()
            plt.title('χ2 ' + self.phase_name[i-1])
            plt.xlabel(self.labels[0])
            plt.ylabel(self.labels[1])

            min_redchi2 = np.nanmin(redchi2_2D)
            max_redchi2 = np.nanmax(redchi2_2D)
            step = (max_redchi2 - min_redchi2) / 10
            print('The minimum χ2 value for', self.phase_name[i-1], 'is:', min_redchi2)

            if min_redchi2 <= 1 and min_redchi2 != max_redchi2:
                contoured = plt.contour(redchi2_2D, levels=np.arange(1, max_redchi2, step), colors="white", linewidths=0.5, origin="lower", extent=[min(self.x), max(self.x), min(self.y), max(self.y)], norm=colors.LogNorm())
                plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.1f')
            elif min_redchi2 > 1 and min_redchi2 != max_redchi2:
                contoured = plt.contour(redchi2_2D, levels=np.arange(np.round(min_redchi2, decimals=1)+0.1, max_redchi2, step), colors="white", linewidths=0.5, origin="lower", extent=[min(self.x), max(self.x), min(self.y), max(self.y)], norm=colors.LogNorm())
                plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.1f')

            os.makedirs(self.output_dir, exist_ok=True)
            plt.savefig(os.path.join(self.output_dir, "χ2_" + self.phase_name[i-1] + '.pdf'), format='pdf')
            plt.show()
            plt.close()

        return min_redchi2

    def plot_redchi2_tot(self, redchi2):
        redchi2_2D = np.reshape(redchi2, (len(np.unique(self.y)), len(np.unique(self.x))))
        plt.imshow(redchi2_2D, cmap=self.color_scheme, aspect='auto', origin='lower', extent=[min(self.x), max(self.x), min(self.y), max(self.y)], norm=colors.LogNorm())
        plt.colorbar()
        plt.title('Total reduced χ2')
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])

        min_redchi2 = np.nanmin(redchi2_2D)
        max_redchi2 = np.nanmax(redchi2_2D)
        step = (max_redchi2 - min_redchi2) / 10

        if min_redchi2 <= 1:
            contoured = plt.contour(redchi2_2D, levels=np.arange(1, max_redchi2, step), colors="white", linewidths=0.5, origin="lower", extent=[min(self.x), max(self.x), min(self.y), max(self.y)], norm=colors.LogNorm())
        else:
            contoured = plt.contour(redchi2_2D, levels=np.arange(np.round(min_redchi2, decimals=1)+0.1, max_redchi2, step), colors="white", linewidths=0.5, origin="lower", extent=[min(self.x), max(self.x), min(self.y), max(self.y)], norm=colors.LogNorm())

        plt.clabel(contoured, inline=True, fontsize=10, fmt='%1.1f')
        print('The minimum value of the total reduced χ2 is:', min_redchi2)

        n_p = np.count_nonzero(self.y == self.y[0])
        min_redchi2_pos = np.where(redchi2_2D == np.nanmin(redchi2_2D))
        min_redchi2_y = self.y[min_redchi2_pos[0]*n_p]
        min_redchi2_x = self.x[min_redchi2_pos[1]]
        print('The ' + self.labels[0] + ' and ' + self.labels[1] + ' position of the minimum total reduced χ2 is:', min_redchi2_x , ',', min_redchi2_y)

        os.makedirs(self.output_dir, exist_ok=True)
        plt.savefig(os.path.join(self.output_dir, "redχ2_tot.pdf"), format='pdf')
        plt.show()
        plt.close()

        return min_redchi2

    def Qcmp_elem(self):
        Qcmp_elem = np.empty(len(self.data[self.labels[2]]))
        for i in range(len(self.labels[2:])):
            Model_elem = np.array(self.data[self.labels[i+2]])
            for j in range(len(Model_elem)):
                Qcmp_elem[j] = self.Q_elem(self.apfu_obs[i], Model_elem[j], self.obs_err[i])
            self.plot_elem(Qcmp_elem, i)

    def Qcmp_phase(self):
        phase_idx = []
        Qcmp_phase_tot = np.empty((len(self.y), len(np.unique(self.phase_id))))

        for i in np.unique(self.phase_id):
            idx = np.where(self.phase_id == i)[0]+2
            phase_idx.append(idx)
            #print(idx)
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
                    print('WARNING: The number of elements in', self.phase_name[i-1], 'is 1, the χ2 may not be accurate')

            redchi2_phase_tot[:, i-1] = redchi2_phase
            self.min_redchi2 = np.append(self.min_redchi2, min_redchi2_phase)

        return redchi2_phase_tot

    def redchi2_tot(self):
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
        """Calculate the total quality factor (Q*cmp) for the modelled composition."""

        # Ensure phase_id exists and is part of self
        if self.phase_id is None:
            raise ValueError("phase_id is not set. Ensure it is initialized.")

        # All phases have the same weight
        weight = np.ones(len(np.unique(self.phase_id)))

        # Normalize the weight
        weight_norm = self.norm_weight(weight)

        # Calculate the quality factor for each row of mod_elem
        Qcmp_allphases = np.empty(len(Qcmp_phase_tot[:, 0]))

        # Calculate the total quality factor for each row of the modelled composition
        for i in range(len(Qcmp_phase_tot[:, 0])):
            Qcmp_allphases[i] = self.Q_tot(Qcmp_phase_tot[i, :], weight_norm)

        # Find the position of the maximum value of the quality factor
        max_Qcmp = np.where(Qcmp_allphases == np.nanmax(Qcmp_allphases))

        # Accessing the redchi2_phase at the location of the maximum Qcmp_allphases_weight
        Qcmpmax_redchi2_value = redchi2_phase_tot[max_Qcmp[0][0]]

        print("The reduced χ2 values for the phases at the maximum Q*cmp is:", Qcmpmax_redchi2_value)

        # Plot the results
        self.plot_tot(Qcmp_allphases, "Unweighted")


    def Qcmp_tot_weight(self, Qcmp_phase_tot, redchi2_phase_tot):
        """Calculate the weighted total quality factor (Q*cmp_weight) for the modelled composition."""

        # Calculate weight based on the minimum reduced χ2 value of each phase
        self.min_redchi2[self.min_redchi2 < 1] = 1
        weight = 1 / self.min_redchi2
        print('The weight is: ', weight)

        # Normalize the weight
        weight_norm = self.norm_weight(weight)
        print('The normalized weight fraction is: ', weight_norm)

        # Calculate the quality factor for each row of mod_elem
        Qcmp_allphases_weight = np.empty(len(Qcmp_phase_tot))

        # Calculate the quality factor total for each row of the modelled composition
        for i in range(len(Qcmp_phase_tot)):
            Qcmp_allphases_weight[i] = self.Q_tot(Qcmp_phase_tot[i, :], weight_norm)

        # Find the position of the maximum value of the quality factor
        max_Qcmp = np.where(Qcmp_allphases_weight == np.nanmax(Qcmp_allphases_weight))

        # Accessing the redchi2_phase at the location of the maximum Qcmp_allphases_weight
        Qcmpmax_redchi2_value = redchi2_phase_tot[max_Qcmp[0][0]]
        print("The reduced χ2 values for the phases at the maximum Q*cmp are:", Qcmpmax_redchi2_value)

        # Plot the results
        self.plot_tot(Qcmp_allphases_weight, "Weighted")

        return Qcmp_allphases_weight


#! this allows to test the file only if it is run as the main file
if __name__ == "__main__":
    QFA = QualityFactorAnalysis()
    QFA.run_analysis(0)
