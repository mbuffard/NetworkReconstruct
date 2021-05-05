Phos2Net V1.1
==================
this version add a log file in order to save the parameters used for the reconstruction


==================

Reconstruct a network associated to a list of target proteins, based on KEGG and Pathway Common databases.
This code is an improvement of the code from:
Naldi A, Larive RM, Czerwinska U, Urbach S, Montcourrier P, Roy C, et al. (2017) Reconstruction and signal propagation analysis of the Syk signaling network in breast cancer cells. PLoS Comput Biol 13(3): e1005432
http://dx.doi.org/10.1371/journal.pcbi.1005432
GUI and options have been added.
It contains all databases needed.

This code, like the publication, is free to use and adapt, as long as the original work is cited (CC BY 4.0).


Requirements
============

* r-base (used for random walk)
* python 3 (may also work with 2.7)
* networkx
* fisher (fast implementation used in rank_pathways.py). Use this for python 3.5 (until next release):
  pip install git+https://github.com/brentp/fishers_exact_test.git

A conda environment with all requirements can be created with the following commands:

```
conda create -n bio python
source activate bio

conda install networkx biopython
pip install git+https://github.com/brentp/fishers_exact_test.git
conda install r
```

Other useful packages:

* pandas (R-like dataframe)
* matplotlib and seaborn (plotting)
* jupyter (notebook)

They can be added to your conda environment with:
```
conda install pandas seaborn jupyter
```

Program will start by calling Phos2Net.py

Source files
============

* converter.py: protein identifier mapping handles Uniprot and HUGO (gene names) data mappings (uniprot covers KEGG IDs, HUGO provides gene names)
* rank_pathways.py: load the lists of members of imported pathways (see import_kegg),
  load the lists of identified proteins and perform enrichment tests (fisher)
  It can also rank the pathways according to pvalue and select the ones to build a network
* nearshortest.py: extraction of shortest paths (needs cleanup)
* parsego.py: load the GO association data and extract associations to specific sets of terms
* GUI.py: Graphical User Interface using tkinter

