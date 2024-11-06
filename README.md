# IntersecT
<img width="150" alt="IntersecT_logo" src="https://github.com/neoscalc/IntersecT/blob/main/src/IntersecT/Logo-IntersecT.png">

IntersecT is a Python script that allows quantitative isopleth thermobarometry on Perple_X pressure (P)-temperature (T)-composition (X) phase diagrams.
The quality factor of the composition (***Q*cmp**, as reported in [Duesterhoeft and Lanari, 2020](https://doi.org/10.1111/jmg.12538)) is calculated for the extracted compositions using the WERAMI routine of Perple_X. The **WERAMI output** file must contain only the phases considered in equilibrium, and their composition must be in **a.p.f.u.** to allow the propagation of uncertainties.

For a step-by-step tutorial on calculating phase diagrams with Perple_X, refer to the [Groppo & Castelli tutorial](https://www.perplex.ethz.ch/perplex/tutorial/Castelli_Groppo_Torino_Tutorial/previous_versions/Castelli_Groppo_Torino_Perple_X_691_Tutorial.pdf)

## How to run IntersecT
IntersecT requires Python v3.10 or more recent versions

Install Jupyter Notebook.
Run the following to install:

```python
pip install jupyterlab
```
```python
pip install notebook
```

Download the [Intersect Jupyter Notebook](https://github.com/neoscalc/IntersecT/blob/main/IntersecT.ipynb) file and save it in the Python\Scripts directory.

In the Python\Script directory, double-click on the **jupyter-notebook.exe** file - it will open a browser tab which is a local folder - and run the **IntersecT.ipynb** file

Read the instructions and run the code cells (**Shift+Enter** when positioned in the cell, or **Run all**)

## Input files
Also refer to the [example folder](https://github.com/neoscalc/IntersecT/tree/main/example)

1. **WERAMI .tab file**

Using the WERAMI routine of Perple_X, export the a.p.f.u. of the elements for each phase that are considered at equilibrium. Columns can be added or deleted manually (e.g., using Excel), but do not modify the rows.

2. **User's .txt file**

The .txt input file must contain:
i. As many rows of comment text as wanted, just make sure they start with ‘#’; these will not be read by the script, and can contain information about the sample or the calculation
        
ii. Any given name for the considered element; these will serve as titles for the plots and the saved .pdf files

iii. The measured composition (in a.p.f.u.) for each element corresponding to the WERAMI calculation
        
iv.	The observed uncertainty (1σ) for each element; alternatively, “-“ must be inserted to let the script calculate an uncertainty based on the type of acquired analysis
        
v. Type of analysis: EDS, WDS map, WDS spot
        
vi. The name of the phases; these will serve as titles for the plots and the saved .pdf files
        
vii. A chosen colour scheme for the plots from the options available in Python

## Explanation of the code cells
```python
%pip install pyIntersecT
```
installs the required packages automatically if not already installed in the local environment.

```python
from pyIntersecT import IntersecT 
InT = IntersecT.QualityFactorAnalysis()
```
imports the **IntersecT library**, which is necessary to run the following calculations and defines the object InT to call the subroutines to calculate the quality factor.

```python
InT.import_output_perplex()
```
imports the **WERAMI .tab file**. It identifies the x and y axes, the labels of the WERAMI columns, and how many columns for each phase are present. 

```python
InT.import_analytical_compo()
```
imports the **.txt file** with the observed compositions and recognizes the input values; if instead of the observed uncertainty (1σ) the user reports “-“, the script will calculate an uncertainty based on the type of measurement (EDS, WDS map, WDS spot).

```python
InT.set_output_directory()
```
lets the user choose a **directory** in which the script will save all the output maps in .pdf format.

```python
InT.Qcmp_elem(i)
```
calculates the ***Q*cmp** for each **element** of each phase.

```python
Qcmp_phase_tot = InT.Qcmp_phase(i)
```
calculates the ***Q*cmp** for each **phase**.

```python
redchi2_phase = InT.redchi2_phase()
```
calculates the **reduced chi-squared statistic** (best fit for values ≤ 1) for each **phase** and finds the lowest value; if only two elements are considered for a phase, the result will be a **chi-squared statistic** (best fit ~ 0).

```python
redchi2_allphases = InT.redchi2_tot()
```
calculates the **total reduced chi-squared statistic** for the **overall dataset**, which can be helpful for inversion.

```python
Qcmp_allphases = InT.Qcmp_tot(Qcmp_phase_tot, redchi2_phase)
```
calculates the **total *Q*cmp** for the **overall dataset**, assuming an equal weight for each phase (i.e., an unweighted total *Q*cmp).

```python
Qcmp_allphases_weight = InT.Qcmp_tot_weight(Qcmp_phase_tot, redchi2_phase)
```
calculates a weighted **total *Q*cmp** for the **overall dataset** based on the best value of the reduced chi-squared statistic for each phase.
