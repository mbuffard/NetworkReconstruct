from __future__ import print_function

import os
import sys
import math
import networkx as nx
import subprocess


############################### Creating the graph ############################
def load_graph( file_name, sep='\t'):
    """Load a graph from a file.
    
    Each line contains an interaction as a (source,target,weight) triplet (tab-separated by default)"""
    with open(file_name) as f:
        arc_list = []
        for row in f:
            row = row.strip('\n').split(sep)
            w = float(row[2])
            arc_list.append( (row[0], row[1], w) )
    f.close()
    return load_interactions(arc_list)


def load_interactions( arc_list ):
    """Build a graph from a list of interactions.
    
    Each interaction is a triplet (source, target,weight)"""
    G = nx.DiGraph()
    G.add_weighted_edges_from(arc_list)
    return G

def load_names( file_name, sep='\t' ):
    names = {}
    with open(file_name) as f:
        for row in f:
            row = row.strip('\n').split(sep)
            names[row[0]] = row[1]
    return names


def rank_paths(G, source):
    "Find the length of the best path for each reachable node"
    return nx.single_source_dijkstra_path_length(G, source)


def find_paths(G, ranks, targets, tfile, overflow=0.1, title="shortests", outfolder='subnetworks'):
    "Use the previously identified ranks to reconstruct the (near-)shortest paths"
    selected_nodes = {}
    selected_edges = {}
    to_explore = {}
    
    for t in targets:
        if t not in ranks: continue
        to_explore[t] = ranks[t] * overflow
    
    while to_explore:
        next_explore = {}
        for t,o in to_explore.items():
            if t not in ranks: continue
            selected_nodes[t] = o
            cur_best = ranks[t]
#            print( "%s (%s)" % ( t, cur_best ))
            for n in G.predecessors(t):
                if n not in ranks: continue # not reachable from the source
                prev_best = ranks[n]
                prev_distance = G.get_edge_data(n,t)['weight']
                cur_o = prev_best +  prev_distance - cur_best
                if cur_o > o: continue # too much longer than the shortest path
#                print( "  %s:   %s (+%s)" % (n,prev_best, cur_o))
                cur_score = o - cur_o
                selected_edges[ (n,t) ] = overflow - cur_score
                best_known = -1
                if n in selected_nodes: best_known = selected_nodes[n]
                if n in to_explore: best_known = max(best_known, to_explore[n])
                if n in next_explore: best_known = max(best_known, next_explore[n])
                if cur_score > best_known:
                    next_explore[n] = cur_score
            
        to_explore = next_explore


    if not os.path.exists(outfolder): os.mkdir(outfolder)
    with open( os.path.join(outfolder, "%s_%s_edges.tsv" % (title,tfile)), 'w' ) as out:
        out.write("Source\tTarget\tType\tSign\tOverflow\tWeight\n")
        for e in selected_edges:
            w = G.get_edge_data(e[0],e[1])['weight']
            out.write("%s\t%s\t%s\t%s\n" % (e[0],e[1],selected_edges[e],w))
    with open( os.path.join(outfolder, "%s_%s_nodes.tsv" % (title,tfile)), 'w' ) as out:
        out.write("UID\tLabel\tBest\tOverflow\n")
        for n in selected_nodes:
            b = ranks[n]
            out.write("%s\t%s\t%s\n" % (n, b, overflow-selected_nodes[n]))

    return selected_edges, selected_nodes



def save_graph(G, filename):
    with open(filename, 'w') as out:
        for source,target,data in G.edges(data=True):
            w = data['weight']
            out.write('%s\t%s\t%s\n' % (source,target,w))

def random_walk(G, helper, source, rpath):
    "Refine the weights of all arcs after running a random walk (external R code)"
    
    filename = os.path.join(rpath,'network')
    f_nodes = '%s_nodes.tsv' % filename
    f_wedges = '%s_wedges.tsv' % filename
    f_wresults = '%s_wresults.tsv' % filename
    f_redges = '%s_redges.tsv' % filename
    f_rresults = '%s_rresults.tsv' % filename
    with open(f_nodes, 'w') as out:
        for node in G.nodes:
            if node == source:
                v = 100
            else: v = 0
            out.write('%s %s\n' % (node, v))

    with open(f_wedges, 'w') as out:
        for source,target,data in G.edges(data=True):
            p = helper(data['weight'])
            out.write('%s %s %s\n' % (source,p,target))

    with open(f_redges, 'w') as out:
        for source,target in G.edges():
            out.write('%s %s %s\n' % (source,1,target))
    
    # Launch the random walk for the weighted and raw networks
    f_rcode = os.path.join(rpath, 'random_walk.r')
    fd = open(f_rcode)
    subprocess.call(['R', '--slave', '--args', f_nodes, f_wedges, f_wresults], stdin=fd)
    fd.close()
    fd = open(f_rcode)
    subprocess.call(['R', '--slave', '--args', f_nodes, f_redges, f_rresults], stdin=fd)
    fd.close()
    
    # TODO: post-processing
    scale = {}
    with open(f_wresults) as f:
        for line in f:
            data = line.strip('\n').split(' ')
            scale[data[0]] = float(data[1])
    with open(f_rresults) as f:
        for line in f:
            data = line.strip('\n').split(' ')
            scale[data[0]] /= float(data[1])

    scaled_arcs = []
    for source,target,data in G.edges(data=True):
        w = scale[source] * data['weight']
        scaled_arcs.append( (source,target,w) )

    Gs = nx.DiGraph()
    Gs.add_weighted_edges_from(scaled_arcs)
    return Gs


if __name__ == '__main__':

    overflow = .2
    overflow = 0
    
    target_files = [ 'ezrin', 'cortactin' ]
    
    print( "Load")
    G = load_graph("full/distance_edges.tsv")
    names = load_names("full/nodes.tsv")
    
    print( "Dijkstra")
    ranks = rank_paths(G, 'P43405')
    
    for tfile in target_files:
        print( "#############################################################")
        print( "#       ", tfile)
        print( "#############################################################")
        thetargets = []
        with open( "%s.txt" % tfile ) as f:
            for line in f:
                t = line.strip('\n').split('\t')[0]
                if t in ranks:
                    thetargets.append(t)
        
        print( "Reconstruct shortests")
        se,sn = find_paths(G, ranks, thetargets, tfile, overflow, outfolder='sub_test')

