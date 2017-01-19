NetworkReconstruct
##################


Requirements (using miniconda is encouraged)
============================================


* python 3 (may also work with 2.7)
* networkx
* biopython (for the KEGG (KGML) parser in import_kegg.py)
* fisher (fast implementation used in rank_pathways.py). Use this for python 3.5 (until next release):
  pip install git+https://github.com/brentp/fishers_exact_test.git


conda create -n bio python
source activate bio

conda install networkx biopython
pip install git+https://github.com/brentp/fishers_exact_test.git


Other useful packages:

* pandas (R-like dataframe)
* matplotlib and seaborn (plotting)
* sklearn (machine learning, PCA)
* jupyter (notebook)


conda install pandas seaborn jupyter



Code files
==========

* converter.py: protein identifier mapping handles Uniprot and HUGO (gene names) data mappings (uniprot covers KEGG IDs, HUGO provides gene names)
* import_kegg.py: parse and convert KEGG pathways (from KGML to TSV)
* rank_pathways.py: load the lists of members of imported pathways (see import_kegg),
  load the lists of identified proteins and perform enrichment tests (fisher)
  It can also rank the pathways according to pvalue and select the ones to build a network
* syk_datasets.py: define input and output and calls code from rank_pathways
* nearshortest.py: extraction of shortest paths (needs cleanup)

