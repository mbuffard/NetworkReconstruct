import converter
import os


#This script will extract the proteins expressed in a particular cell line from CCLE
#cell_line="MCF7_BREAST"

def extractData(cell_line):

	f=open("global_proteomic/summed_sn_non_normalized_formated.csv","r")

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
		
		if column[cell_line_index]!="NA":
			#extract uniprot ID
			prot_id=column[0].split("|")[1]
			#extract the isoforme extension
			prot_id=prot_id.split("-")[0]
			
			if converter.handler.clean_uid(prot_id) :
				prot_id=converter.handler.clean_uid(prot_id)
				expressed_proteins.add(prot_id)
			else:
				if converter.handler.import_symbol(column[1]):
					#check if it's Human for contamination elimination
					origin=column[0].split("|")[2]
					origin=origin.split("_")[1]
					
					if origin=='HUMAN':
						prot_id=converter.handler.import_symbol(column[1])
						expressed_proteins.add(prot_id)
									

	return expressed_proteins

