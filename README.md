## Description
The following script perform singular value decomposition (SVD) the chosen and extracts the fist eigenvector. This process mimics the process of 'fitting' a vector to 3D data or a 3D regression. This is particularly useful for measuring angle between two alpha helicies in a simulation dataset. Alternatively, it can also used to angle fluctuations with respect to a particular axis.

Prerequisite to the script are `numpy`, `scipy` and `matplotlib` (optional). The script can be directly modified for use. 

## Usage 

* **Selection :** First an appropriate part(s) of the protein needs to be selected to perform the calculation : 

**Examples**

range of residues : e.g. selecting 100 to 120 residues => `'protein and resid 100:120'`

range of atoms : e.g. selecting 2102 to 2485 residues => `'protein and bynum 2102:2485'`

* **Choosing the files** : 
choose the names of the files `mda.Universe('step7_production.gro','step7_production.xtc')`. Uncommment the plotting it is necessary.

* **Running the script** : `Python angle_estimator.py`
