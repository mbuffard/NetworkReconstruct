"Parse downloaded KEGG pathways and save them in a simpler format"

from __future__ import print_function

import os.path
import sys
import os

import json
from Bio.KEGG.KGML import KGML_parser

import converter


cachedir = os.path.join(os.getenv("HOME"), "pathwaycache", "kegg_reboot")

kgmldir = os.path.join(cachedir, "kgml")
outdir = os.path.join(cachedir, "pathways")
if not os.path.exists(outdir):
    os.makedirs(outdir)


INFO_SIGN = {
    'activation' : "+",
    'inhibition' : "-",
}

INFO_TYPE = {
    'phosphorylation'       : 'phosphorylation',
    'dephosphorylation'     : 'phosphorylation',
    'methylation'           : 'modification',
    'ubiquitination'        : 'modification',
    'binding/association'   : 'association',
    'dissociation'          : 'dissociation',
    'expression'            : 'expression',
    'repression'            : 'repression',
    'state change'          : 'modification',
}

INFO_IGNORED = set((
'missing interaction',
'indirect effect',
'compound',
))

def uid_of(entry):
    if entry.type == "gene":
        return '%s' % entry.id

    # TODO: integrate relations involving compounds?
    if False and entry.type == "compound":
        return entry.id

f_members = open(os.path.join(outdir, 'members.txt') , 'w')
for filename in os.listdir(kgmldir):
    nodes = {}
    pathway_nodes = set()
    if not filename.endswith(".xml") or not filename.startswith("hsa"): continue
    kgid = filename[:-4]
    kgml_file = os.path.join(kgmldir, filename)
    f = open(kgml_file)
    print(filename)
    parsed = KGML_parser.read(f)
    f.close()
    
    pathway_name = parsed.title
    
    # Find and map components
    components = set()
    # Enter the parsed pathway info and extract components
    for k,entry in parsed.entries.items():
        if entry.type != "gene":
            continue

        rnames = entry.name.split(" ")
        names = [ converter.handler.to_uniprot(converter.KEGG_IDX, u) for u in rnames ]
        names = [ u for u in names if u ]
        if not names:
            # TODO: can we fix empty KEGG nodes?
            print( "empty node??", rnames)
            continue
        components.update(names)
        nodes['%s' % entry.id] = names


    # Extract relations and save in normalized format
    out = open(os.path.join(outdir, '%s.txt' % kgid), 'w')
    relations = []
    for rel in parsed.relations:
        uid1 = uid_of(rel.entry1)
        uid2 = uid_of(rel.entry2)
        if not uid1 or not uid2:
            continue

        # Extract the type and sign
        reltype = ''
        sign = ''
        for name,tag in rel.subtypes:
            if name in INFO_SIGN:
                sign = INFO_SIGN[name]
            elif name in INFO_TYPE:
                reltype = INFO_TYPE[name]
        relations.append( (uid1,uid2, reltype, sign) )
        if uid1 not in nodes or uid2 not in nodes:
            continue
        
        for n1 in nodes[uid1]:
            for n2 in nodes[uid2]:
                out.write('%s\t%s\t%s\t%s\n' % (n1, n2, reltype, sign))
                pathway_nodes.add(n1)
                pathway_nodes.add(n2)
    
    f_members.write("%s\t%s\t%s\n" % (kgid, pathway_name, ','.join(pathway_nodes)))
    out.close()


f_members.close()

