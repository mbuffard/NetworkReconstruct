from __future__ import print_function

import os
import sys

mapping_dir = os.path.join(os.getenv("HOME"), 'pathwaycache', 'mapping')
mapping_file = os.path.join(mapping_dir, "uniprot_ready.txt")

# source files
hgnc_file = os.path.join(mapping_dir, "hgnc_complete_set.txt")
uniprot_file = os.path.join(mapping_dir, "HUMAN_9606_idmapping.dat")
reviewed_filename = os.path.join(mapping_dir,"uniprot_reviewed_ids.txt")

KEGG_IDX = 0
HGNC_IDX = 1


class UniprotMapper:
    "This one loads the final clean mapping"
    
    def __init__(self):
        if not os.path.exists(mapping_file):
            handler = UniprotHandler()

        self.mapping = {}
        self.reverse = []
        
        f = open(mapping_file)
        headers = f.readline().strip().split("\t")
        for h in headers[1:]:
            self.reverse.append( {} )
        self.duplicates = self.reverse[-1]
        nb_maps = len(self.reverse)
        for line in f:
            data = line.split("\t")
            uid = data[0]
            terms = [ d.strip().split(",") for d in data[1:] ]
            for idx in range(nb_maps):
                for t in terms[idx]:
                    self.reverse[idx][t] = uid
            self.mapping[uid] = terms
        f.close()

        self.h2s = h2s = {}
        self.s2h = s2h = {}
        f = open(hgnc_file)
        f.readline()
        for line in f:
            line = line.strip().split("\t")
            if len(line) > 1:
                h2s[ line[0] ] = line[1]
                s2h[ line[1] ] = line[0]
        f.close()
        


    def to_uniprot(self, idx, external):
        mapping = self.reverse[idx]
        if external in mapping:
            return mapping[external]
    
    
    def to_external(self, idx, uniprot):
        if uniprot in self.duplicates:
            return self.to_external(idx, self.duplicates[uniprot])
        
        if uniprot in self.mapping:
            return self.mapping[uniprot][idx]
    
    def get_symbols(self, idx, uids):
        symbols = set()
        for uid in uids:
            s = self.to_external(idx, uid)
            if s:
                for v in s:
                    symbols.add(v)
        
        return symbols


    def to_symbol(self, uid):
        hgncs = self.to_external(HGNC_IDX, uid)
        if hgncs:
            hgnc = hgncs[0]
            if hgnc in self.h2s:
                return self.h2s[hgnc]

    def import_symbol(self, s):
        try:
            h = self.s2h[s]
            return self.reverse[HGNC_IDX][h]
        except:
            return None

    def import_symbols(self, symbols):
        uniprots = set()
        for s in symbols:
            u = self.import_symbol(s)
            if u: uniprots.add(u)
        return uniprots

    def clean_uid(self, uid):
        if uid in self.mapping:
            return uid
        
        if uid in self.reverse[-1]:
            return self.reverse[-1][uid]


class UniprotHandler:
    "This one is in charge of duplicate identification"
    
    def __init__(self):
        if not os.path.exists("mapping/uniprot_mapped.txt"):
            import_uniprot_mapping()
        
        f = open("mapping/uniprot_mapped.txt")
        self.mapping = {}
        self.duplicates = {}
        ext_columns = f.readline().strip().split('\t')[2:]
        nbmaps = len(ext_columns)
        
        tmpmap = {}
        reviewed = []
        unreviewed = []
        for line in f:
            info = line.strip().split('\t')
            
            uid = info[0]
            rev = info[1]
            links = []
            
            for l in info[2:]:
                if l == 'None':
                    values = ()
                else:
                    values = tuple(set(l.split(',')))
                links.append(values)
            
            links = tuple(links)
            tmpmap[uid] = links
            if rev == 'R':
                reviewed.append(uid)
            else:
                unreviewed.append(uid)
        f.close()
    
        # integrate into main and reverse mappings
        self.reverse_mappings = []
        for db in range(nbmaps):
            self.reverse_mappings.append( {} )
        self.reverse_mappings.append( {} )
        
        self.integrate_mapping(reviewed, tmpmap, nbmaps)
        self.integrate_mapping(unreviewed, tmpmap, nbmaps)
        
        out = open(mapping_file, "w")
        for uid in self.mapping:
#            if uid 
            out.write("%s" % uid)
            for terms in self.mapping[uid]:
                out.write("\t")
                if terms:
                    out.write( ",".join(terms) )
            out.write("\n")
        out.close()
    
    
    def integrate_mapping(self, uids, tmpmap, nbmaps):
        print( 'Integration of %s items' % len(uids) )
        for idx in range(nbmaps):
            cur_map = self.reverse_mappings[idx]
            for uid in uids:
                for curid in tmpmap[uid][idx]:
                    if curid in cur_map:
                        # Mark it as duplicate
                        if uid not in self.duplicates:
                            ref = cur_map[curid]
                            self.duplicates[uid] = ref
                            self.mapping[ref][nbmaps].add( uid )
                            if uid in self.mapping:
                                # Take it out of the main mapping and transfer existing links to the right place
                                deprecated = self.mapping[uid]
                                newplace   = self.mapping[ref]
                                for oldidx in range(len(deprecated)):
                                        for v in deprecated[oldidx]:
                                            newplace[oldidx].add( v )
                                            self.reverse_mappings[oldidx][v] = ref
                                #TODO: can some noise may remain in the self.duplicates map?
                                del self.mapping[uid]
                    else:
                        if uid not in self.mapping:
                            self.mapping[uid] = [ set() for a in range(nbmaps+1) ]
                        self.mapping[uid][idx].add( curid )
                        cur_map[curid] = uid
        print( 'Integrated' )


def import_uniprot_mapping(databases=["KEGG", "HGNC"]):

    mapping = {}
    dbs = {}
    idx = 0
    for dbname in databases:
        dbs[dbname] = idx
        idx += 1
    
    reviewed = set()
    f = open(reviewed_filename)
    for line in f:
        reviewed.add(line.strip())
    
    nbdupl = 0
    f = open(uniprot_file)
    for line in f:
        uid, db, extid = line.strip().split()
        try:
            idx = dbs[db] + 1
        except:
            continue
        
        uid = uid.split("-")[0]
        if uid not in mapping:
            if uid in reviewed:
                rev = "R"
            else:
                rev = "U"
            mapping[uid] = [rev, None, None]
        
        m = mapping[uid]
        if m[idx]:
            nbdupl += 1
            m[idx] += ","+extid
        else:
            m[idx] = extid
    f.close()
    
    print( "%s duplicates?" % nbdupl )
    
    # TODO: check integration of duplicate detection and cleanups
    
    # save the clean mapping
    out = open("mapping/uniprot_mapped.txt", "w")
    out.write('Uniprot\tReviewed')
    for db in dbs:
        out.write('\t%s' % db)
    out.write('\n')
    
    skel = "\t%s" * (1 + len(dbs))
    for uid in mapping:
        links = skel % tuple(mapping[uid])
        out.write("%s%s\n" % (uid, links))
    out.close()


handler = UniprotMapper()

if __name__ == "__main__":
    print( len(handler.mapping), "items" )
#        print "Reverse:"
#        for r in handler.reverse:
#            print "  ", len(r)
#        print
#        print 'hsa:6850 -->',   handler.to_uniprot(0, 'hsa:6850')
#        print 'HGNC:1067 --> ', handler.to_uniprot(1, 'HGNC:1067')

#        print 'P43405 -->', handler.to_external(0, 'P43405')
#        print 'P43405 -->', handler.to_external(1, 'P43405')
#        print 'P30450 -->', handler.duplicates['P30450']
#        print 'P30450 -->', handler.to_external(1, 'P30450')
#        print 'P30450 -->', handler.to_external(2, 'P30450')

#        print 'HGNC:4931 --> ', handler.to_uniprot(1, 'HGNC:4931')

