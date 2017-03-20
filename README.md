NetworkReconstruct
==================

Reconstruct a network associated to a list of target proteins, based on enriched pathways.
For a complete description and use case, see the associated publication:
Naldi A, Larive RM, Czerwinska U, Urbach S, Montcourrier P, Roy C, et al. (2017) Reconstruction and signal propagation analysis of the Syk signaling network in breast cancer cells. PLoS Comput Biol 13(3): e1005432
http://dx.doi.org/10.1371/journal.pcbi.1005432

This code, like the publication, is free to use and adapt, as long as the original work is cited (CC BY 4.0).


Requirements
============


* python 3 (may also work with 2.7)
* networkx
* biopython (for the KEGG (KGML) parser in import_kegg.py)
* fisher (fast implementation used in rank_pathways.py). Use this for python 3.5 (until next release):
  pip install git+https://github.com/brentp/fishers_exact_test.git

A conda environment with all requirements can be created with the following commands:

conda create -n bio python
source activate bio

conda install networkx biopython
pip install git+https://github.com/brentp/fishers_exact_test.git


Other useful packages:

* pandas (R-like dataframe)
* matplotlib and seaborn (plotting)
* jupyter (notebook)


conda install pandas seaborn jupyter



Source files
============

* converter.py: protein identifier mapping handles Uniprot and HUGO (gene names) data mappings (uniprot covers KEGG IDs, HUGO provides gene names)
* import_kegg.py: parse and convert KEGG pathways (from KGML to TSV)
* rank_pathways.py: load the lists of members of imported pathways (see import_kegg),
  load the lists of identified proteins and perform enrichment tests (fisher)
  It can also rank the pathways according to pvalue and select the ones to build a network
* syk_datasets.py: define input and output and calls code from rank_pathways
* nearshortest.py: extraction of shortest paths (needs cleanup)
* parsego.txt: load the GO association data and extract associations to specific sets of terms

