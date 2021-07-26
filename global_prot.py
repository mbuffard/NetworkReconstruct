import converter
import os


#This script will extract the proteins expressed in a particular cell line from CCLE
#cell_line="MCF7_BREAST"

def extractData(cell_line):

	f=open("global_proteomic/summed_sn_non_normalized_formated.csv","r")
	g=open("global_proteomic/"+cell_line+"expressed_proteins.txt","w")
	h=open("global_proteomic/"+cell_line+"globalresults.txt","w")

	heads=f.readline().strip().split()


	#search for cell line index in the data
	i=0
	for head in heads:
		if str(head)==cell_line:
			cell_line_index=i
			break

		i+=1

	expressed_proteins=set()
	for line in f:
		column=line.strip().split("\t")
		h.write(column[0].split("|")[1].split("-")[0]+"\t"+str(column[cell_line_index])+"\n")
		
		if column[cell_line_index]!="NA":
			#extract uniprot ID
			prot_id=column[0].split("|")[1]
			#extract the isoforme extension
			prot_id=prot_id.split("-")[0]
			
			if converter.handler.clean_uid(prot_id) :
				prot_id=converter.handler.clean_uid(prot_id)
				expressed_proteins.add(prot_id)
				g.write(prot_id+"\t"+str(column[cell_line_index])+"\n")
			else:
				if converter.handler.import_symbol(column[1]):
					#check if it's Human for contamination elimination
					origin=column[0].split("|")[2]
					origin=origin.split("_")[1]
					
					if origin=='HUMAN':
						prot_id=converter.handler.import_symbol(column[1])
						expressed_proteins.add(prot_id)
						g.write(prot_id+"\t"+str(column[cell_line_index])+"\n")
									

	return expressed_proteins