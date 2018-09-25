from __future__ import print_function

import os
import os.path
import math
import fisher
import networkx as nx

import converter


cachedir = os.path.join(os.getenv("HOME"), "pathwaycache", "all_pathways_KEGG&PC2", "pathways")
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
        pmember_names = [ converter.handler.to_symbol(uid) for uid in pmembers ]
        C  = len(members)
        Cn = len(pmembers)
        score = fisher.pvalue_population(Cn, C, Pn, P).two_tail
        log_score = math.log(score, 10)
        if (float(Cn)/C) < Pr: 
            log_score = -log_score
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (uid,C,Cn,log_score, name, ','.join(pmembers), ','.join(pmember_names)))
        
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
    all_data = []
    pid_name = {}
    for filename in targetfiles:
        f = open( os.path.join(outfolder, "%s__pathways.tsv" % filename) )
        cur_data = {}
        all_data.append(cur_data)
        for line in f:
            data = line.strip('\n').split('\t')
            pid = data[0]
            members = data[5].strip()
            member_names = data[6].strip()
            cur_data[pid] = (data[3],members, member_names)
            pid_name[pid] = data[4]
        f.close()
    
    all_keys = set(all_data[0].keys())
    for cur_scores in all_data:
        all_keys.intersection_update(cur_scores.keys())

    # write all scores
    out = open(os.path.join(outfolder,outname), 'w')
    header = "PID\tName"
    for filename in targetfiles:
        header += '\t%s score\t%s members\t%s names' % (filename,filename,filename)
    out.write('%s\n' % header.strip())
    for key in all_keys:
        name = pid_name[key]
        out.write("%s\t%s" % (key,name))
        for cur_data in all_data:
            score,members, member_names = cur_data[key]
            out.write("\t%s\t%s\t%s" % (score,members, member_names))
        out.write('\n')
    out.close()


def build_network(rank_file, network_file_name, targets):
    network_file = "%s.tsv" % network_file_name
    nodes_file = "%s_nodes.tsv" % network_file_name
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
        #print(pid,score,members)
    
    out = open(network_file, 'w')
    out.write('Source\tTarget\tRel\tSign\tOrigin\n')
    thenodes = set()
    for src,tgt,data in net.edges(data=True):
        thenodes.add(src)
        thenodes.add(tgt)
        origin = ','.join(data['origin'])
        rel = ','.join(data['rel'])
        sign = data['sign']
        out.write('%s\t%s\t%s\t%s\t%s\n' % (src,tgt,rel, sign, origin))
    out.close()
    
    out = open(nodes_file, 'w')
    out.write('PID\tName\tDataset\n')
    for uid in thenodes:
        name = converter.handler.to_symbol(uid)
        if uid in targets:
            out.write('%s\t%s\tx\n' % (uid, name))
        else:
            out.write('%s\t%s\t\n' % (uid, name))
    out.close()


def merge_pathway_in_network(net, pid):
    pathway_file = os.path.join(cachedir, '%s.txt' % pid)
    f = open(pathway_file)
    #print (pathway_file)
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

