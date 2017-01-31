from __future__ import print_function

import os
import os.path
from rank_pathways import *

outfolder = "SYK_output/"
targetfolder = "/home/aurelien/Dropbox/Work/Physica files/results_SYK/dataset-comparison/nodes/"
targetfiles = ('mda231_positif.txt', 'mcf7_positif.txt', 'mda231-mcf7_positif.txt', 'mda231.txt', 'mcf7.txt', 'mda231-mcf7.txt', 'dg75.txt', 'mapped_skov3tr_yu15.tsv')
#targetfiles = ('mda231_positif.txt', 'mcf7_positif.txt', 'dg75.txt')


all_targets = []
for filename in targetfiles:
    targets = load_targets(os.path.join(targetfolder, filename))
    all_targets.append(targets)
    rank_file = os.path.join(outfolder, "%s__pathways.tsv" % filename)
    network_file_name = os.path.join(outfolder, "%s__network" % filename)

    rank(targets, pathways, allmembers, rank_file )
    build_network(rank_file, network_file_name, targets )


merge_ranks(outfolder, targetfiles, "merged_pathways.tsv")


