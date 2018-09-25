from __future__ import print_function

import os
import os.path
import networkx as nx
from rank_pathways import *
import nearshortest
import parsego
import shutil
import converter



# Set paths and other parameters
SYK_UID = 'P43405'

work_folder = '../MDA231'
outfolder = os.path.join(work_folder, "output/")

filename = 'mda231_positif_corr.csv'
#create folder for the cell line
#os.mkdir(os.path.join('../',filename.split(".")[0]))
#print (filename.split(".")[0])

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
    #"cancer": ("GO:0046222", "GO:0046223", "GO:0050813", "GO:0044598", "GO:0044597", "GO:0045122", "GO:0050814", "GO:0018983", "GO:0097681"),
    "pTyrMod": ("GO:0004725", "GO:0004713"),
    "yphosphatase": ("GO:0004725",),
    "ykinase": ("GO:0004713",),
}



#target_terms_folders={"all_targets":}
all_go_targets = [ t for t in target_terms ]
print(all_go_targets)


# Rules to assign weights
WEIGHTS = {
    "--":5,
    "-T":2,
    "-t":8,
    "d-":3,
    "dT":1,
    "dt":6,
}

# overflow defined twice 
overflow = 0.2
#overflow = 0
target_files_folder = '../targets'
#target_files = [ 'cortactin' ]


# add weights to every interaction depending on the source and target
# return arc_liste which contains source target and associated weight
def add_weights(filename, indata, targets, ptyr_mod):
    arc_list = []
	# open syk_published.tsv_network.tsv
    with open( filename ) as f:
        header = f.readline().strip('\n').split('\t')
		#extract the source and target node from the file
        for line in f:
            row = line.strip('\n').split('\t')
            source_node = row[0]
            target_node = row[1]
            if "None" == source_node or "None" == target_node: continue
            
            if source_node in indata:#
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
    print ("dans la fonction je vaut :"+str(len(ptyr_mod)))
    return arc_list

#not used?
def add_directs(arc_list, directs):
    for t in directs:
        w = WEIGHTS['dT']
        arc_list.append( ('P43405',t,w) )


targets = load_targets(os.path.join(work_folder, filename))
print ("nombre de cibles "+ str(len(targets)))
# write the pathway file with score
rank(targets, pathways, allmembers, rank_file )
build_network(rank_file, network_file_name, targets )



##############################
# Add a x to node in targets
##############################

node=open(os.path.join(outfolder, "%s__network_nodes.tsv" % filename),"r")
#f=open(work_folder+"/targets.txt",'w')
node_file={}
for line in node:
    node_file[line.strip('\n').split('\t')[0]]=line.strip('\n').split('\t')[1]

node.close()




###############################
# Add additional_edges
###############################
node=open(os.path.join(outfolder, "%s__network_nodes.tsv" % filename),"a")

net=open (os.path.join(outfolder, "%s__network.tsv" % filename),"a")
with open( os.path.join(work_folder, "additional_edges_gene.txt")) as neighb:
    #ligne=neighb.readline()
    for line in neighb:
        sourceID = line.strip('\n').split('\t')[0]
        targetID = line.strip('\n').split('\t')[2]
        rel = line.strip('\n').split('\t')[4]
        net.write(sourceID+"\t"+targetID+"\t"+rel+"\t "+"\t"+"additional_edges"+"\n")
        # add node in node file if not yet
        if sourceID not in node_file:
            node.write(sourceID+"\t"+line.strip('\n').split('\t')[1]+"\t"+" "+"\n")
            node_file[sourceID]=line.strip('\n').split('\t')[1]
        if targetID not in node_file:
            node_file[targetID]=line.strip('\n').split('\t')[3]
            node.write(targetID+"\t"+line.strip('\n').split('\t')[3]+"\t"+" "+"\n")     

node.close()
net.close()
neighb.close()

#f.close()


##########################
# Assign edge weights
##########################

# Retrieve GO processes
print('load GO')
uid2go = parsego.load_goa(goa_filename)



uniprot2process = parsego.load_targets(all_go_targets, target_terms, uid2go, obo_filename)
categorie=set()
'''for uid in uid2go:
    go_file.write(uid+"\t")
    categorie=set()
    for go in uid2go[uid]:
        if go in target_terms:
            categorie.add(target_terms[go])
            go_file.write(go+"")
    print (uid+":")
    for cat in categorie:
        print (cat+"\n")
go_file.close()'''

#classify all the targets in categories
ptyr_mod = set([ uid for uid in uniprot2process if 'pTyrMod' in uniprot2process[uid] ])
adh_motility=set([ uid for uid in uniprot2process if ('cell adhesion and motility' in uniprot2process[uid] and uid in targets)])
growth_death=set([ uid for uid in uniprot2process if ('cell growth and death' in uniprot2process[uid] and uid in targets)])
transport_metabolism=set([ uid for uid in uniprot2process if ('cell transport and metabolism"' in uniprot2process[uid] and uid in targets)])
IS_inflammation=set([ uid for uid in uniprot2process if ('immune system and inflammation' in uniprot2process[uid] and uid in targets)])
differentiation=set([ uid for uid in uniprot2process if ('cell differentiation' in uniprot2process[uid] and uid in targets)])
all_targets=[uid for uid in targets]

categorie_list=["ptyr_mod","adh_motility","growth_death","transport_metabolism","IS_inflammation","differentiation"	]
categorie_sets=[ptyr_mod,adh_motility,growth_death,transport_metabolism,IS_inflammation,differentiation,]


#Create the file for each category
i=0
for file in categorie_list:
	print (file)
	catfile=open ("../targets/"+file+".txt","w")
	

	for t in categorie_sets[i]:
		catfile.write(t+"\t")
		for cat in uniprot2process[t]:
			catfile.write(cat+",")
			catfile.write("\n")
	catfile.close()
	i+=1
ptyr_file=open ("../targets/all_targets.txt","w")
for t in all_targets:
	ptyr_file.write(t+"\n")
ptyr_file=open("../targets/ptyr.txt","w")
for t in ptyr_mod:
	ptyr_file.write(t+"\t")
	ptyr_file.write(str(converter.handler.to_symbol(t))+"\n")

#write the files and add it in target files



print ("en dehors de la fonction je vaut :"+str(len(ptyr_mod)))
print ("et adh motility vaut :"+str(len(adh_motility)))


#parsego.integrate(target_files,targets,)

# TODO: pick phospho-targets from the full target dataset
phosphotargets = targets
arc_list = add_weights('%s.tsv' % network_file_name, targets, phosphotargets, ptyr_mod)
print ("diff :"+str(len(targets))+" et "+str(len(phosphotargets)))


# Load it as a graph and compute the length of the all shortest paths from SYK
print('start shortest')
G = nx.DiGraph()
G.add_weighted_edges_from(arc_list)
bests = nx.single_source_dijkstra_path_length(G, SYK_UID)
print ("best :")
meilleur_cibles=[]
for b in bests:
	meilleur_cibles.append(b)
print (str(len(meilleur_cibles)))

# Reconstruct shortest paths to specified sets of targets, targets need to be in a file with the uniprot
for tfile in os.listdir(target_files_folder):
    print( "#############################################################")
    print( "#       ", tfile)
    print( "#############################################################")
    thetargets = []
    with open( os.path.join("../targets/", "%s" % tfile) ) as f:
        for line in f:
            t = line.strip('\n').split('\t')[0]
            if t in bests:
                thetargets.append(t)
    
    print( "Reconstruct shortests")
    #parcourir les differentes catégories et chercher + court chemin verifier nearshortest function
    se,sn = nearshortest.find_paths(G, bests, thetargets, tfile, overflow, outfolder=outfolder)


#print ("nombre de targets :",len(thetargets))
'''go_file=open(os.path.join(outfolder, "go_associated_targets.txt"), "w")
for tfile in os.listdir(target_files_folder):

    #thetargets = []
    with open( os.path.join(work_folder, "%s" % tfile) ) as f:
        for line in f:
            t = line.strip('\n').split('\t')[0]
            
            if not t in uniprot2process:
                continue
            else:
                if not uniprot2process[t]:
                    continue
                else: 
                    go_file.write(t+"\t")
                    #print (t+":")                  
                    for cat in uniprot2process[t]:
                        go_file.write(cat+",")
                        
                    go_file.write("\n")'''

'''w=0
for i in thetargets:
    if i in uniprot2process:
        w+=1
        print (i+":")
        for cat in uniprot2process[i]:
            print (cat)
print ("go target associated"+str(w))'''

#print (os.sys.argv[0])



#### copy the script used in the output folder #####
shutil.copyfile(os.sys.argv[0],os.path.join(outfolder, os.sys.argv[0]))

print (len(uniprot2process))

