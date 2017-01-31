from __future__ import print_function

import os
import os.path
from rank_pathways import *

outfolder = "../output/"
targetfolder = ".."
filename = 'syk_published.tsv'


targets = load_targets(os.path.join(targetfolder, filename))
rank_file = os.path.join(outfolder, "%s__pathways.tsv" % filename)
network_file_name = os.path.join(outfolder, "%s__network" % filename)

rank(targets, pathways, allmembers, rank_file )
build_network(rank_file, network_file_name, targets )


