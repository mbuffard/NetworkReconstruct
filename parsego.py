import os
import os.path
import sys
import gzip

# Define a set of GO terms of interest
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


# The OBO file describes the GO hierarchy
# The GOA file gives associations between GO terms and uniprot IDs
obo_filename = "/home/aurelien/go-basic.obo"
goa_filename = "/home/aurelien/Downloads/goa_human.gaf.gz"


def load_goa(filename):
    'Flatten the GOA file into a simple two column file'
    
    if not os.path.exists(filename):
        print('MISSING: GOA file ', filename)
        return None
    
    uid2go = {}
    
    # read the source annotation file
    c = 0
    f = gzip.open(filename, 'rt')
    for line in f:
        if not line.startswith("UniProtKB"): continue
        line = line.strip().split('\t')
        uid = line[1]
        term = line[4]
        qualifier = line[3]
        if qualifier.find('NOT') >= 0: c+=1 ; continue
        if uid not in uid2go:
            uid2go[uid] = set()
        uid2go[uid].add(term)
    f.close()
    
    return uid2go


def find_target(go, terms, go2target):
    if go in go2target:
        return go2target[go]

    target = set()
    if go in terms:
        for parent in terms[go][-1]:
            target.update( find_target(parent, terms, go2target) )
    go2target[go] = target
    return target


def load_targets(all_targets, target_terms, uid2go, obo_filename):
    'Find all UIDs associated to the selected GO terms'
    
    go2target = {}
    for t in target_terms:
        for go in target_terms[t]:
            thistarget = set( (t,) )
            if go not in go2target:
                go2target[go] = thistarget
            else:
                go2target[go].update( thistarget )

    terms = {}
    namespaces = set()

    in_bloc = False
    f = open(obo_filename)
    for line in f:
        line = line.strip()

        if not in_bloc:
            if line == "[Term]":
                in_bloc = True
                cur_id = None
                name = None
                namespace = False
                alternatives = []
                parents = []
            continue
        
        if not line:
            if namespace == "biological_process" or namespace == "molecular_function":
                terms[cur_id] = (namespace, name, alternatives, parents)
            in_bloc = False
            continue
        
        info = line.split(":", 1)
        key = info[0].strip()
        value = info[1].strip()

        if key == "id":
            cur_id = value
        elif key == "name":
            name = value
        elif key == "namespace":
            namespace = value
            namespaces.add(namespace)
        elif key == "alt_id":
            value = value.split(" ")[0]
            alternatives.append( value )
        elif key == "is_a":
            value = value.split(" ")[0]
            parents.append( value )
    f.close()

    for go in go2target:
        if go not in terms:
            print( "%s is missing" % go)

    # now associate targets to uniprot IDs
    uniprot2target = {}
    for uniprot in uid2go:
        for go in uid2go[uniprot]:
            targets = find_target(go, terms, go2target)
            if not targets:
                continue
            
            if uniprot not in uniprot2target:
                # make a copy to avoid messing with the source data!
                uniprot2target[uniprot] = set( targets )
            else:
                uniprot2target[uniprot].update(targets)
    
    return uniprot2target


def integrate(files, all_targets, uniprot2target):
    'Add columns with the selected GO annotations to a set of files'
    
    for filename in files:
        f = open("%s.tsv" % filename )
        out = open("%s_go.tsv" % filename, "w")
        line = f.readline()[:-1]
        out.write(line)
        for t in all_targets:
            out.write("\t%s" % t)
        out.write("\n")
        for line in f:
            line = line[:-1]
            uid = line.split("\t")[0]
            if uid not in uniprot2target:
                targets = set()
            else:
                targets = uniprot2target[uid]
            
            out.write(line)
            for t in all_targets:
                if t in targets:
                    out.write("\tx")
                else:
                    out.write("\t-")
            out.write("\n")
        f.close()
        out.close()

def integrate_folders(folders, all_targets, uniprot2target):
    'Find all files matching "*_nodes.tsv" in a group of folders and call integrate on them'
    files = []
    for thefolder in folders:
        for f in os.listdir(thefolder):
            if not f.endswith('_nodes.tsv'): continue
            files.append( os.path.join(thefolder, f[:-4]) )

    integrate(files, all_targets, uniprot2target)


if __name__ == '__main__':

    all_targets = [ t for t in target_terms ]
    
    uid2go = load_goa(goa_filename)
    uniprot2target = load_targets(all_targets, target_terms, uid2go)

if False:
    # integrate it
    folders = [
    #    '/home/aurelien/Dropbox/Work/Target_network/syk_networks',
        '/home/aurelien/Dropbox/Work/Target_network/syk_latest',
        '/home/aurelien/Dropbox/Work/Target_network/syk_latest/source_data'
    ]

    integrate_folders(folders, all_targets, uniprot2target)

