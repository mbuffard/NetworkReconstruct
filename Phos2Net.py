from __future__ import print_function
import fisher2
import os
import os.path
import networkx as nx
import nearshortest
from rank_pathways import *
import parsego
import shutil
import converter
import GUI
from tkinter import *
from tkinter.messagebox import *
from tkinter import ttk
from os.path import basename
from GUI import interface
import threading

####################################################################################################################################################################
# This code is an improvement of the code from Naldi A, Larive RM, Czerwinska U, Urbach S, Montcourrier P, Roy C, et al. (2017), PLoS Comput Biol 13(3)
# new options and GUI have been added
################################################################################################################################################################""




if interface.check==0:
    interface.destroy()


 
# Set paths and other parameters
# setup paths and filename

def infowindow(message1,message2):
	progres = Tk()
	label=Label(progres,text=message1,font=("Helvetica", 24)   )
	label.grid(row=0,column=0)
	label2=Label(progres,text=message2,font=("Helvetica", 24)   )
	label2.grid(row=1,column=0)
	progres.mainloop()


progreswindow = threading.Thread(target=infowindow, args=("Job in progress : ","Please Wait"),daemon= True)
progreswindow.start()

work_folder = interface.work_folder
rpath = os.path.join( 'random_walk')
network_folder = os.path.join(work_folder,'networks')
outfolder = os.path.join(work_folder, "output/")
if not os.path.isdir(outfolder):
	os.makedirs(outfolder)
target_files_folder = './targets/'

filename_location=interface.filename
filename=basename(interface.filename)
source = converter.handler.clean_uid(interface.source.get())

###################change from Aurelien #################
rank_file = os.path.join(outfolder, "%s__pathways.tsv" % filename)
network_file_name = os.path.join(outfolder, "%s__network" % filename)
#############################################################

# The OBO file describes the GO hierarchy, the GOA files associates terms to Uniprot IDs
obo_filename = os.path.join('./pathwaycache', 'GO', "go-basic.obo")
goa_filename = os.path.join('./pathwaycache', 'GO', "goa_human.gaf.gz")


############################################################################################
# Load the custom list of modifiers to promote edges according to user selection

p_mod = set()
if interface.addKinPhos.get()==1:
	if 0 in interface.KinPhos_mod:
		f = open(os.path.join('Modifiers/Tyr_Kinase.txt'),'r')
		for line in f:
			data=line.strip('\n').split('\t')
			p_mod.add(data[0])
		f.close()
	if 1 in interface.KinPhos_mod:
		f=open(os.path.join('Modifiers/Tyr_phosphatase.txt'),'r')
		for line in f:
			data = line.strip('\n').split('\t')
			p_mod.add(data[0])
		f.close()
	if 2 in interface.KinPhos_mod:	
		f=open(os.path.join('Modifiers/ST_kinase.txt'),'r')
		for line in f:
			data = line.strip('\n').split('\t')
			p_mod.add(data[0])
		f.close()
	if 3 in interface.KinPhos_mod:	
		f=open(os.path.join('Modifiers/ST_phosphatase.txt'),'r')
		for line in f:
			data = line.strip('\n').split('\t')
			p_mod.add(data[0])
		f.close()
if interface.addSpecMod.get()==1:
	for uniprot in interface.Mod_list.get().split():
		uniprot
		uniprot=converter.handler.clean_uid(uniprot)
		p_mod.add(uniprot)
###################################################################################################




# the groups of GO terms of interest
target_terms = {
	"cell adhesion and motility": ("GO:0048870","GO:0007155", "GO:0034330",  'GO:0022610', 'GO:0060352', 'GO:0030030', ), # localisation: 'GO:0030054'
	"cell growth and death": ("GO:0008283", "GO:0007049", "GO:0008219", "GO:0019835", "GO:0000920", 'GO:0007569', 'GO:0051301', 'GO:0060242'),
	"cell transport and metabolism": ("GO:0044237", "GO:1990748", "GO:0010496"),
	"immune system and inflammation": ("GO:0002376", "GO:0001906"),
	"cell differentiation": ("GO:0030154", "GO:0036166"),
}



all_go_targets = [ t for t in target_terms ]


# Rules to assign weights
WEIGHTS = {
	"--":5,
	"-T":2,
	"-t":8,
	"d-":3,
	"dT":1,
	"dt":6,
}

################version Aurelien ############
def ptyr_score(source_node, target_node):
	if source_node in targets:
		wkey = "d"
	else:
		wkey = "-"

	if source_node in p_mod and target_node in targets:
		wkey += "T"
	elif source_node not in p_mod and target_node in targets:
		wkey += "t"
	else:
		wkey += "-"

	return WEIGHTS[wkey]

#############################################

def get_rate(weight):
	return 10 - weight

#convert unreviewed uniprot name into a reviewed one with same HGNC and load the targets
targets = load_targets( filename_location)



#select the database according to formula value and write the pathway file with score
if interface.databaseKEGG.get()==1 and interface.databasePC.get()==0:
	cachedir=os.path.join( "pathwaycache", "kegg_reboot", "pathways")
if interface.databaseKEGG.get()==0 and interface.databasePC.get()==1:
	cachedir=os.path.join( "pathwaycache", "pc2", "pathways")
if interface.databaseKEGG.get()==1 and interface.databasePC.get()==1:
	cachedir=os.path.join("pathwaycache", "all_pathways_KEGG&PC2", "pathways")

members = os.path.join(cachedir, "members.txt")
pathways, allmembers = load_pathways(members)
rank(targets, pathways, allmembers, rank_file)

#add selection mode option

if interface.selection.get()==1:
	build_network(rank_file, network_file_name, targets,cachedir,1)
if interface.selection.get()==2:
	build_network(rank_file, network_file_name, targets,cachedir,2)


###############################
# Add additional_edges
###############################
if interface.addSubstOption.get()==1:
	directSubs=(add_direct_from_list(interface.subs.get(),source,outfolder,filename))
	for subs in directSubs:
		#convert unreviewed uniprot name into a reviewed one with same HGNC
		subs=converter.handler.clean_uid(subs)
		targets.add(subs)


##########################
# Assign edge weights
##########################

# Retrieve GO processes
print('load GO')
uid2go = parsego.load_goa(goa_filename)

uniprot2process = parsego.load_targets(all_go_targets, target_terms, uid2go, obo_filename)


#classify all the targets in categories
adh_motility=set([ uid for uid in uniprot2process if ('cell adhesion and motility' in uniprot2process[uid] and uid in targets)])
growth_death=set([ uid for uid in uniprot2process if ('cell growth and death' in uniprot2process[uid] and uid in targets)])
transport_metabolism=set([ uid for uid in uniprot2process if ('cell transport and metabolism' in uniprot2process[uid] and uid in targets)])
IS_inflammation=set([ uid for uid in uniprot2process if ('immune system and inflammation' in uniprot2process[uid] and uid in targets)])
differentiation=set([ uid for uid in uniprot2process if ('cell differentiation' in uniprot2process[uid] and uid in targets)])
all_targets=[uid for uid in targets]

categorie_list=["adh_motility","growth_death","transport_metabolism","IS_inflammation","differentiation"]
categorie_sets=[adh_motility,growth_death,transport_metabolism,IS_inflammation,differentiation,]

for file in os.listdir(target_files_folder):
	os.remove(os.path.join("./targets/", "%s" % file)) 
#Create the file for each category
#parcours indices selectionnés
i=0
for i in interface.CatExtraction:
	file=categorie_list[i]
	print (file)
	catfile=open ("./targets/"+file+".txt","w")	
	for t in categorie_sets[i]:
		catfile.write(t+"\t")
		for cat in uniprot2process[t]:
			catfile.write(cat+",")
			catfile.write("\n")
	catfile.close()
	i+=1

ptyr_file=open ("./targets/all_targets.txt","w")
for t in all_targets:
	ptyr_file.write(t+"\n")
ptyr_file.close()

if interface.target_extract_option.get()==1:
	file=open ("./targets/subset_targets.txt","w")
	for t in interface.subset_targets.get().split():
		file.write(t+"\t   "+"\n")
	file.close()

# Load the list of targets from the dataset
targets = set()
f = open(os.path.join( outfolder, '%s__network_nodes.tsv' % filename))
f.readline()
for line in f:
	data = line.strip('\n').split('\t')
	if data[2]:
		targets.add(data[0])
f.close()

############################################# For all targets (all nodes....)###########################################
# Add the weights and build the graph for path search
interactions = add_weights('%s.tsv' % network_file_name, ptyr_score)
G = nearshortest.load_interactions(interactions)

ranks = nearshortest.rank_paths(G, source)
reachable = G.subgraph(ranks)


# Save the reachable graph and the shortest paths
nearshortest.save_graph(reachable, os.path.join(outfolder, "reachable.tsv"))

# Launch the random walk to refine the weights
Gs = nearshortest.random_walk(reachable, get_rate, source, rpath)

# Compute shortest paths on network after random walk
ranks_rdm = nearshortest.rank_paths(Gs, source)
reachable_rdm = G.subgraph(ranks_rdm)

############################################################################################################################
# Reconstruct shortest paths to specified sets of targets (go associated), targets need to be in a file with the uniprot

for tfile in os.listdir(target_files_folder):

	folder=tfile.split('.')[0]

	if not os.path.isdir(outfolder+folder):
		os.makedirs(outfolder+folder)
	print( "#############################################################")
	print( "#       ", tfile)
	print( "#############################################################")
	thetargets = []
	#extract the targets associated with the go term
	with open( os.path.join("./targets/", "%s" % tfile) ) as f:
		for line in f:
			t = line.strip('\n').split('\t')[0]
			if t in ranks:
				thetargets.append(t)
	#print (thetargets)
	if thetargets:
		print( "Reconstruct shortests")
		

		edges,nodes = nearshortest.find_paths(Gs, ranks_rdm, thetargets, tfile, overflow=0, title="%s_rdm_shortest"%filename, outfolder=outfolder+folder)
		if interface.OverflowOption and interface.overflowValue.get()!=0:
			edges,nodes = nearshortest.find_paths(Gs, ranks_rdm, thetargets, tfile, overflow=interface.overflowValue.get()/100, title="%s_rdm_overflow"%filename, outfolder=outfolder+folder)

############################################################################################################################

fenetre = Tk()
label=Label(fenetre,text="Job done, networks have been created in :",font=("Helvetica", 24)   )
label.grid(row=0,column=0)
label2=Label(fenetre,text=interface.work_folder,font=("Helvetica", 24)   )
label2.grid(row=1,column=0)

fenetre.mainloop()