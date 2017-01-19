from __future__ import print_function

import os
import os.path
import math
import fisher
import networkx as nx

cachedir = os.path.join(os.getenv("HOME"), "pathwaycache", "kegg_reboot", "pathways")
keggmembers = os.path.join(cachedir, "members.txt")

def rank(targets, pathways, allmembers, outname):
    P = len(allmembers)
    realtargets = targets.intersection(allmembers)
    Pn = len(realtargets)
    Pr = float(Pn) / P
    
    # score it all
    out = open(outname, 'w')
    for uid,name,members in pathways:
        pmembers = members.intersection(realtargets)
        C  = len(members)
        Cn = len(pmembers)
        score = fisher.pvalue_population(Cn, C, Pn, P).two_tail
        log_score = math.log(score, 10)
        if (float(Cn)/C) < Pr: log_score = -log_score
        out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (uid,C,Cn,log_score, name, ','.join(pmembers)))
    out.close()

def load_pathways(membersfile):
    f = open(membersfile)
    pathways = []
    allmembers = set()
    for line in f:
        uid,name, members = line.strip('\n').split('\t')
        members = set( members.split(',') )
        allmembers.update(members)
        pathways.append( (uid,name,members) )
    f.close()
    
    return pathways, allmembers

def load_targets(targetfile):
    targets = set()
    f = open(targetfile)
    for line in f:
        targets.add(line.strip().split('\t')[0])
    return targets

def merge_ranks(outfolder, targetfiles, outname):
    all_scores = []
    for filename in targetfiles:
        f = open( os.path.join(outfolder, "%s__pathways.tsv" % filename) )
        cur_scores = {}
        all_scores.append(cur_scores)
        for line in f:
            data = line.strip('\n').split('\t')
            cur_scores[data[0]] = data[3]
        f.close()
    
    all_keys = set(all_scores[0].keys())
    for cur_scores in all_scores:
        all_keys.intersection_update(cur_scores.keys())

    # write all scores
    out = open(os.path.join(outfolder,outname), 'w')
    header = ""
    for filename in targetfiles:
        header += '\t%s' % filename
    out.write('%s\n' % header.strip())
    for key in all_keys:
        out.write(key)
        for cur_scores in all_scores:
            out.write("\t%s" % cur_scores[key])
        out.write('\n')
    out.close()


def build_network(rank_file, network_file):
    ranked = []
    f = open(rank_file)
    for line in f:
        data = line.strip('\n').split('\t')
        pid = data[0]
        members = data[5].strip()
        if not members: continue
        members = set( members.split(',') )
        score = float(data[3])
        ranked.append( (pid, score, members) )

    ranked.sort(key=lambda tup: tup[1])
    covered = set()
    net = nx.DiGraph()
    for pid,score,members in ranked:
        if score > -2:
            if not members.difference(covered):
                continue
        merge_pathway_in_network(net, pid)
        covered.update(members)
        print(pid,score,members)
    
    out = open(network_file, 'w')
    for src,tgt,data in net.edges(data=True):
        origin = ','.join(data['origin'])
        rel = ','.join(data['rel'])
        sign = data['sign']
        out.write('%s\t%s\t%s\t%s\t%s\n' % (src,tgt,rel, sign, origin))
    out.close()


def merge_pathway_in_network(net, pid):
    pathway_file = os.path.join(cachedir, '%s.txt' % pid)
    f = open(pathway_file)
    for line in f:
        src,tgt,rel,sign = line.strip('\n').split('\t')
        current = net.get_edge_data(src,tgt)
        origin = [pid,]
        if rel: rel = set((rel,))
        else:   rel = set()
        if current:
            cursign = current['sign']
            if not sign:
                sign = cursign
            elif cursign and sign != cursign:
                # TODO: improve sign merging?
                sign = '?'
            rel.update(current['rel'])
            origin += current['origin']
        net.add_edge(src,tgt, rel=rel, sign=sign, origin=origin)


pathways, allmembers = load_pathways(keggmembers)

