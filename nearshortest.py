from __future__ import print_function

import os
import sys
import math
import networkx as nx


############################### Creating the graph ############################
def load_graph( file_name, sep='\t'):
    with open(file_name) as f:
        arc_list = []
        for row in f:
            row = row.strip('\n').split(sep)
            w = float(row[2])
            arc_list.append( (row[0], row[1], w) )
    f.close()

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


def find_paths(G, bests, targets, tfile, overflow=0.1, title="shortests", outfolder='subnetworks'):
    selected_nodes = {}
    selected_edges = {}
    to_explore = {}
    
    for t in targets:
        if t not in bests: continue
        to_explore[t] = bests[t] * overflow
    
    while to_explore:
        next_explore = {}
        for t,o in to_explore.items():
            if t not in bests: continue
            selected_nodes[t] = o
            cur_best = bests[t]
#            print( "%s (%s)" % ( t, cur_best ))
            for n in G.predecessors(t):
                if n not in bests: continue # not reachable from the source
                prev_best = bests[n]
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
#    print( selected_nodes)


    if not os.path.exists(outfolder): os.mkdir(outfolder)
    with open( os.path.join(outfolder, "%s_%s_edges.tsv" % (title,tfile)), 'w' ) as out:
        out.write("Source\tTarget\tType\tSign\tOverflow\tWeight\n")
        for e in selected_edges:
            w = G.get_edge_data(e[0],e[1])['weight']
            out.write("%s\t%s\t%s\t%s\n" % (e[0],e[1],selected_edges[e],w))
    with open( os.path.join(outfolder, "%s_%s_nodes.tsv" % (title,tfile)), 'w' ) as out:
        out.write("UID\tLabel\tBest\tOverflow\n")
        for n in selected_nodes:
            b = bests[n]
            name = names[n]
            out.write("%s\t%s\t%s\t%s\n" % (n, name, b, overflow-selected_nodes[n]))

    return selected_edges, selected_nodes



if __name__ == '__main__':

    overflow = .2
    overflow = 0
    
    target_files = [ 'ezrin', 'cortactin' ]
    
    print( "Load")
    G = load_graph("full/distance_edges.tsv")
    names = load_names("full/nodes.tsv")
    
    print( "Dijkstra")
    bests = nx.single_source_dijkstra_path_length(G, 'P43405')
    
    for tfile in target_files:
        print( "#############################################################")
        print( "#       ", tfile)
        print( "#############################################################")
        thetargets = []
        with open( "%s.txt" % tfile ) as f:
            for line in f:
                t = line.strip('\n').split('\t')[0]
                if t in bests:
                    thetargets.append(t)
        
        print( "Reconstruct shortests")
        se,sn = find_paths(G, bests, thetargets, tfile, overflow, outfolder='sub_test')

