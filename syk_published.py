from __future__ import print_function

import os
import os.path
import networkx as nx
from rank_pathways import *
import nearshortest
import parsego


# Set paths and other parameters
SYK_UID = 'P43405'

work_folder = '../syk_published'
outfolder = os.path.join(work_folder, "output/")
filename = 'syk_published.tsv'

rank_file = os.path.join(outfolder, "%s__pathways.tsv" % filename)
network_file_name = os.path.join(outfolder, "%s__network" % filename)

# The OBO file describes the GO hierarchy, the GOA files associates terms to Uniprot IDs
obo_filename = os.path.join(os.getenv("HOME"), 'pathwaycache', 'GO', "go-basic.obo")
goa_filename = os.path.join(os.getenv("HOME"), 'pathwaycache', 'GO', "goa_human.gaf.gz")


# the groups of GO terms of interest
target_terms = {
    "cell adhesion and motility": ("GO:0048870","GO:0007155", "GO:0034330",  'GO:0022610', 'GO:0060352', 'GO:0030030', ), # localisation: 'GO:0030054'
    "cell growth and death": ("GO:0008283", "GO:0007049", "GO:0008219", "GO:0019835", "GO:0000920", 'GO:0007569', 'GO:0051301', 'GO:0060242'),
    "cell transport and metabolism": ("GO:0044237", "GO:1990748", "GO:0010496", "GO:1990748"),
    "immune system and inflammation": ("GO:0002376", "GO:0001906"),
    "cell differentiation": ("GO:0030154", "GO:0036166"),
    "cancer": ("GO:0046222", "GO:0046223", "GO:0050813", "GO:0044598", "GO:0044597", "GO:0045122", "GO:0050814", "GO:0018983", "GO:0097681"),
    "pTyrMod": ("GO:0004725", "GO:0004713"),
    "yphosphatase": ("GO:0004725",),
    "ykinase": ("GO:0004713",),
}
all_go_targets = [ t for t in target_terms ]


# Rules to assign weights
WEIGHTS = {
    "--":5,
    "-T":8,
    "-t":1,
    "d-":7,
    "dT":9,
    "dt":2,
}

overflow = .2
overflow = 0
target_files_folder = ''
target_files = [ 'ezrin', 'cortactin' ]


def add_weights(filename, indata, targets, ptyr_mod):
    arc_list = []
    with open( filename ) as f:
        header = f.readline().strip('\n').split('\t')
        for line in f:
            row = line.strip('\n').split('\t')
            source_node = row[0]
            target_node = row[1]
            if "None" == source_node or "None" == target_node: continue
            
            if source_node in indata:
                wkey = "d"
            else:
                wkey = "-"
            
            if source_node in ptyr_mod and target_node in targets:
                wkey += "T"
            elif source_node not in ptyr_mod and target_node in targets:
                wkey += "t"
            else:
                wkey += "-"
            
            weight = WEIGHTS[wkey]
            arc_list.append( (source_node, target_node, weight) )
        
        f.close()
    
    return arc_list

def add_directs(arc_list, directs):
    for t in directs:
        w = WEIGHTS['dT']
        arc_list.append( ('P43405',t,w) )


targets = load_targets(os.path.join(work_folder, filename))
rank(targets, pathways, allmembers, rank_file )
build_network(rank_file, network_file_name, targets )


##########################
# Assign edge weights
##########################

# Retrieve GO processes
print('load GO')
uid2go = parsego.load_goa(goa_filename)
uniprot2process = parsego.load_targets(all_go_targets, target_terms, uid2go, obo_filename)
ptyr_mod = set([ uid for uid in uniprot2process if 'pTyrMod' in uniprot2process[uid] ])

# TODO: pick phospho-targets from the full target dataset
phosphotargets = targets
arc_list = add_weights('%s.tsv' % network_file_name, targets, phosphotargets, uniprot2process)


# Load it as a graph and compute the length of the all shortest paths from SYK
print('start shortest')
G = nx.DiGraph()
G.add_weighted_edges_from(arc_list)
bests = nx.single_source_dijkstra_path_length(G, SYK_UID)


# Reconstruct shortest paths to specified sets of targets
for tfile in target_files:
    print( "#############################################################")
    print( "#       ", tfile)
    print( "#############################################################")
    thetargets = []
    with open( os.path.join(work_folder, "%s.txt" % tfile) ) as f:
        for line in f:
            t = line.strip('\n').split('\t')[0]
            if t in bests:
                thetargets.append(t)
    
    print( "Reconstruct shortests")
    se,sn = nearshortest.find_paths(G, bests, thetargets, tfile, overflow, outfolder=outfolder)

