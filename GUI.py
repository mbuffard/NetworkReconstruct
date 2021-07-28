from tkinter import filedialog
from tkinter import *
from tkinter.messagebox import *
from tkinter import ttk


class Interface(Frame):
	
	"""Notre fenêtre principale.
	Tous les widgets sont stockés comme attributs de cette fenêtre."""
	
	def __init__(self, fenetre, **kwargs):
		self.databaseKEGG=IntVar()
		self.databasePC=IntVar()   
		self.filename=None
		self.work_folder=None
		
		self.addSubstOption=IntVar()  
		self.addKinPhos=IntVar() 
		self.addSpecMod=IntVar() 
		self.addGlobal=IntVar()
		self.DemoteCCLEOption=IntVar()
		self.OverflowOption=IntVar()
		self.overflowValue =IntVar()
		self.target_extract_option=IntVar()
		self.GOterm_extract_option=IntVar() 
		self.source=StringVar()   
		self.subs=StringVar()
		self.Mod_list=StringVar()
		self.selected_tissue=StringVar()
		self.Cell_line=[]
		self.KinPhos_mod=[]
		self.subset_targets=StringVar()
		self.selection=IntVar()
		self.check=0
		self.CatExtraction=[]
		self.target_extract_option=IntVar()
		self.GOterm_extract_option=IntVar()


		Frame.__init__(self, fenetre, width=768, height=576, **kwargs)
		self.pack(fill=BOTH)

		self.frame1=Frame(self,relief=GROOVE)
		self.frame1.grid(row=0, column=2, columnspan=3, rowspan=19)
		
		
		# function for tooltips creation
		def add_Tooltip(widget,texte='ok'):
			x=y=0
			x,y,cx,cy=self.bbox("insert")
			x+=widget.winfo_rootx()-50
			y+=widget.winfo_rooty()+20
			self.toolTip=Toplevel(self)
			self.toolTip.wm_overrideredirect(True)
			self.toolTip.wm_geometry("+%d+%d" % (x, y))
			#enter the text
			label = Label(self.toolTip, text=texte, justify='left', relief='solid', borderwidth=1,font=("times", ))
			label.pack(ipadx=1)
			
		#function for tooltip destruction
		def destroy_Tooltip(self):
			self.toolTip.destroy()


		
		#Input file selection
		self.inputFile_label = Label(self, text="Choose target file\n (with phosphorylated proteins) :",font=("Helvetica"),)
		self.inputFile_label.grid(row=1,column=0,pady=10)

		self.doc_InputFile=Label(self,text='?',font=("Helvetica", 	))
		self.doc_InputFile.grid(row=1, column=1)
		
		self.doc_InputFile.bind('<Enter>', lambda event, text='The file should be in csv format \nwith one uniprot accession number (AC).\nThis format is available in common spreadsheets.': add_Tooltip(self.doc_InputFile,text))
		self.doc_InputFile.bind('<Leave>', lambda event :destroy_Tooltip(self))

		
		
		self.inputFile_button = Button(self, text="Browse...", command=self.cliquerFile,font=("Helvetica"))
		self.inputFile_button.grid(row=1, column=2,pady=10)

		self.labelFile=Label(self,text="")
		self.labelFile.grid(row=2,column=2)

		#Directory for Output Files
		self.inputFolder_label = Label(self, text="Choose output folder :",font=("Helvetica"))
		self.inputFolder_label.grid(row=3, column=0,pady=10)

		self.inputFolder_button = Button(self, text="Browse...", command=self.cliquerFolder,font=("Helvetica"))
		self.inputFolder_button.grid(row=3, column=2 ,pady=10)

		self.labelFolder=Label(self,text='')
		self.labelFolder.grid(row=4,column=2)

		#check buttons for database selection
		self.selectionDatabase_label = Label(self, text="Choose database(s) :",font=("Helvetica"))
		self.selectionDatabase_label.grid(row=5, column=0)

		self.Check_Kegg = Checkbutton(self, text="KEGG", variable=self.databaseKEGG,onvalue=1,font=("Helvetica"))
		self.Check_Kegg.grid(row=5, column=2,sticky='w')
		
		self.Check_PC = Checkbutton(self, text="Pathway Commons", variable=self.databasePC,font=("Helvetica"))
		self.Check_PC.grid(row=6, column=2,sticky='w',)


		#radiobutton for selection Mode
		self.selectionMode_label = Label(self, text="Choose pathway's selection mode :",font=("Helvetica"))
		self.selectionMode_label.grid(row=7, column=0)

		self.label_docSelectMode=Label(self,text='?',font=("Helvetica", ), activebackground='yellow')
		self.label_docSelectMode.grid(row=7, column=0, sticky='e')
		self.label_docSelectMode.bind('<Enter>', lambda event, text='In the pathway selection mode, the "enriched pathways" option will select \nonly the most enriched pathways for each target. \nThe all pathways, will select all the pathways from the selected database(s).': add_Tooltip(self.label_docSelectMode,text))
		self.label_docSelectMode.bind('<Leave>', lambda event :destroy_Tooltip(self))


		self.Enriched_mode = Radiobutton(self, text="Enriched pathways", variable=self.selection,value=1,font=("Helvetica"))
		self.Enriched_mode.grid(row=7, column=2,sticky='w')
		
		self.All_mode = Radiobutton(self, text="All pathways", variable=self.selection,value=2,font=("Helvetica")   )
		self.All_mode.grid(row=8, column=2,sticky='w')

		#option to add extra direct substrates
		self.add_dirrectSubst = Checkbutton(self, variable=self.addSubstOption,text="Add source direct substrates :", font=("Helvetica"),command=self.Able_subst)
		self.add_dirrectSubst.grid(row=9, column=0)

		self.label_docAddSubs=Label(self,text='?',font=("Helvetica",))
		self.label_docAddSubs.grid(row=9, column=0, sticky='e')
		self.label_docAddSubs.bind('<Enter>', lambda event, text='The option "add source direct substrates" will add direct edges coming from source to substrate,\n and promote the passage via substrates as it do for targets. \nSubstrates should be added with their uniprot AC and separated by spaces.': add_Tooltip(self.label_docAddSubs,text))
		self.label_docAddSubs.bind('<Leave>', lambda event :destroy_Tooltip(self))

		self.directSubs=Entry(self,state=DISABLED,textvariable=self.subs,font=("Helvetica"))
		self.directSubs.grid(row=9,column=2)

		#Source Node
		self.SourceLabel=Label(self,text="Enter source node :",font=("Helvetica"))
		self.SourceLabel.grid(row=10,column=0)

		self.SourceInput=Entry(self,textvariable=self.source,font=("Helvetica"))
		self.SourceInput.grid(row=10,column=2)  

		#Add weights
		self.WeightLabel=Label(self,text="Weight the edges, enter list of :",font=("Helvetica"))
		self.WeightLabel.grid(row=11,column=0,columnspan=1,pady=5)

		self.label_docWeight=Label(self,text='?',font=("Helvetica", ))
		self.label_docWeight.grid(row=11, column=0, sticky='e')
		self.label_docWeight.bind('<Enter>', lambda event, text='To weight the edges one can promote the passage via particular proteins: \nfrom the established list (tyrosine kinase/tyrosine phosphatase/...)\nand/or a personalized list of uniprot AC (separated by spaces)': add_Tooltip(self.label_docWeight,text))
		self.label_docWeight.bind('<Leave>', lambda event :destroy_Tooltip(self))


		self.KinPhos=Checkbutton(self,text="Kinases/phosphatases",command=self.Able_KinPhos,variable=self.addKinPhos,font=("Helvetica"))
		self.KinPhos.grid(row=12,column=0)  

		self.yDefilB = Scrollbar(self, orient='vertical')
		self.yDefilB.grid(row=12, column=3, sticky='ns' )
		self.KinPhos_list=Listbox(self,font=("Helvetica"),selectmode=MULTIPLE,yscrollcommand=self.yDefilB.set,exportselection=0)


		self.yDefilB['command'] = self.KinPhos_list.yview
		for item in ["Tyrosine kinase","Tyrosine phosphatase","Serine/threonine kinase", "Serine/threonine phosphatase"]:
			self.KinPhos_list.insert(END, item)
		self.KinPhos_list.config(width=0,height=2,state=DISABLED)
		self.KinPhos_list.grid(row=12,column=2)
		

		self.SpecModOption=Checkbutton(self,text="Add specific proteins to promote",command=self.Able_Go_term,variable=self.addSpecMod,font=("Helvetica"))
		self.SpecModOption.grid(row=13,column=0)  

		self.SpecModOption_list=Entry(self,state=DISABLED,font=("Helvetica"),textvariable=self.Mod_list)
		self.SpecModOption_list.grid(row=13,column=2)

		listCCLE = ['autonomic ganglia', 'biliary tract', 'bone', 'breast', 'central nervous system', 'endometrium', 'haematopoietic and lymphoid tissue', 'kidney', 'large intestine', 'liver', 'lung', 'oesophagus', 'ovary', 'pancreas', 'pleura', 'prostate', 'skin', 'soft tissue', 'stomach', 'thyroid', 'upper aerodigestive tract', 'urinary tract']
		cellLineCCLEdico = {'prostate': ['22RV1', 'DU145', 'LNCAPCLONEFGC', 'PC3', 'VCAP'], 'haematopoietic and lymphoid tissue': ['697', 'CMK', 'EOL1', 'F36P', 'HDMYZ', 'HEL', 'HEL9217', 'JEKO1', 'JM1', 'JURKAT', 'K562', 'KARPAS299', 'KARPAS422', 'KASUMI1', 'KASUMI2', 'KMS11', 'KMS12BM', 'KMS27', 'L428', 'MOLM13', 'MOLM16', 'MONOMAC1', 'MONOMAC6', 'NALM6', 'NCIH929', 'NUDHL1', 'OCIAML5', 'OCILY3', 'OPM2', 'RCHACV', 'REC1', 'REH', 'RL', 'RPMI8226', 'SEM', 'SUDHL4', 'SUDHL6', 'TF1', 'THP1', 'U937'], 'kidney': ['769P', '786O', 'A498', 'A704', 'CAKI1', 'CAKI2', 'KMRC1', 'KMRC20', 'OSRC2', 'TUHR4TKB', 'VMRCRCW'], 'thyroid': ['8305C', '8505C'], 'skin': ['A101D', 'A2058', 'A375', 'C32', 'COLO679', 'COLO741', 'COLO829', 'HS294T', 'HS695T', 'HS944T', 'IGR1', 'IGR37', 'IGR39', 'IPC298', 'K029AX', 'LOXIMVI', 'MELJUSO', 'MEWO', 'RPMI7951', 'RVH421', 'SH4', 'SKMEL2', 'SKMEL28', 'SKMEL3', 'SKMEL30', 'SKMEL5', 'UACC257', 'UACC62', 'WM115', 'WM1799', 'WM2664', 'WM793', 'WM88'], 'central nervous system': ['A172', 'DAOY', 'GAMG', 'GB1', 'KNS42', 'KNS81', 'LN18', 'LN229', 'SF295', 'SNB19', 'SNU1105', 'SW1783', 'U118MG', 'U87MG'], 'soft tissue': ['A204', 'G401', 'G402', 'HT1080', 'KYM1', 'RD', 'RH30', 'RH41'], 'ovary': ['A2780', 'CAOV3', 'COV362', 'FUOV1', 'HEYA8', 'IGROV1', 'JHOS2', 'KURAMOCHI', 'NIHOVCAR3', 'OV56', 'OV90', 'OVCAR4', 'OVCAR8', 'OVSAHO', 'RMUGS', 'SNU119', 'TYKNU'], 'lung': ['A549', 'ABC1', 'CALU1', 'CALU6', 'CHAGOK1', 'CORL105', 'CORL23', 'CORL47', 'CORL88', 'DMS114', 'DMS273', 'DV90', 'EBC1', 'HCC15', 'HCC1833', 'HCC44', 'HCC827', 'HCC95', 'IALM', 'LCLC103H', 'LK2', 'LU65', 'LUDLU1', 'LXF289', 'NCIH1048', 'NCIH1155', 'NCIH1299', 'NCIH1355', 'NCIH1435', 'NCIH1437', 'NCIH146', 'NCIH1568', 'NCIH1573', 'NCIH1581', 'NCIH1650', 'NCIH1666', 'NCIH1693', 'NCIH1703', 'NCIH1792', 'NCIH1793', 'NCIH1944', 'NCIH196', 'NCIH1975', 'NCIH2009', 'NCIH2030', 'NCIH2066', 'NCIH2110', 'NCIH2122', 'NCIH2126', 'NCIH2170', 'NCIH2172', 'NCIH2228', 'NCIH226', 'NCIH2286', 'NCIH2291', 'NCIH23', 'NCIH292', 'NCIH3255', 'NCIH358', 'NCIH441', 'NCIH446', 'NCIH460', 'NCIH520', 'NCIH522', 'NCIH647', 'NCIH650', 'NCIH661', 'NCIH838', 'PC14', 'RERFLCMS', 'RERFLCSQ1', 'SBC5', 'SHP77', 'SKLU1', 'SQ1', 'SW1271', 'SW1573'], 'bone': ['A673', 'SAOS2', 'SJSA1', 'SKES1', 'SKNMC', 'TC71', 'U2OS'], 'stomach': ['AGS', 'HGC27', 'HUG1N', 'IM95', 'KATOIII', 'LMSU', 'MKN1', 'MKN45', 'MKN7', 'NCIN87', 'NUGC3', 'OCUM1', 'SNU1', 'SNU719'], 'pancreas': ['ASPC1', 'BXPC3', 'CFPAC1', 'DANG', 'HUPT3', 'HUPT4', 'KP2', 'KP4', 'L33', 'MIAPACA2', 'PANC0203', 'PANC0403', 'PANC1', 'PATU8988T', 'PL45', 'QGP1', 'SU8686', 'SUIT2', 'SW1990', 'TCCPAN2'], 'breast': ['AU565', 'BT20', 'BT549', 'CAL120', 'CAL51', 'CAL851', 'CAMA1', 'EFM19', 'EFM192A', 'HCC1143', 'HCC1187', 'HCC1395', 'HCC1500', 'HCC1806', 'HCC1937', 'HCC1954', 'HCC2218', 'HCC38', 'HCC70', 'HDQP1', 'JIMT1', 'KPL1', 'MCF7', 'MDAMB157', 'MDAMB231', 'MDAMB436', 'MDAMB453', 'MDAMB468', 'T47D', 'ZR751'], 'upper aerodigestive tract': ['BICR22', 'BICR6', 'CAL27', 'CAL33', 'DETROIT562', 'FADU', 'HSC3', 'HSC4', 'PECAPJ34CLONEC12', 'SCC25'], 'large intestine': ['CCK81', 'CL34', 'COLO205', 'COLO320', 'COLO678', 'HCC56', 'HCT116', 'HCT15', 'HT115', 'HT29', 'HT55', 'LS180', 'LS411N', 'LS513', 'MDST8', 'NCIH716', 'NCIH747', 'RKO', 'SKCO1', 'SNU61', 'SNUC1', 'SNUC2A', 'SNUC5', 'SW1417', 'SW403', 'SW48', 'SW480', 'SW620', 'SW837', 'SW948'], 'endometrium': ['HEC108', 'HEC1A', 'HEC251', 'HEC265', 'HEC50B', 'HEC59', 'HEC6', 'ISHIKAWAHERAKLIO02ER', 'JHUEM2', 'MFE280', 'MFE296', 'MFE319', 'SNGM', 'SNU685'], 'liver': ['HEP3B217', 'HEPG2', 'HLF', 'HUH1', 'HUH6', 'HUH7', 'JHH1', 'JHH4', 'JHH5', 'JHH6', 'JHH7', 'SKHEP1', 'SNU423', 'SNU449'], 'urinary tract': ['HT1197', 'HT1376', 'J82', 'JMSU1', 'KU1919', 'RT112', 'RT4', 'T24', 'TCCSUP', 'UBLC1', 'UMUC3'], 'oesophagus': ['KYSE150', 'KYSE180', 'KYSE30', 'KYSE410', 'KYSE450', 'KYSE510', 'KYSE70', 'OE33', 'TE1', 'TE10', 'TE11', 'TE14', 'TE4', 'TE6'], 'pleura': ['MSTO211H', 'NCIH2052'], 'autonomic ganglia': ['SKNAS'], 'biliary tract': ['SNU1079']}

	
		#option to add global proteomic data for a specified cell line
		self.addGlobalProteomic=Checkbutton(self,text="Add global proteomic data from CCLE :",command=self.Able_cellLine,variable=self.addGlobal,font=("Helvetica"))
		self.addGlobalProteomic.grid(row=14,column=0)

		# Documentation
		self.label_docCCLEoption=Label(self,text='?',font=("Helvetica", ))
		self.label_docCCLEoption.grid(row=14, column=0, sticky='e')
		self.label_docCCLEoption.bind('<Enter>', lambda event, text='The CCLE option will take in account the protein expression in the selected cell line to reconstruct the network': add_Tooltip(self.label_docCCLEoption,text))
		self.label_docCCLEoption.bind('<Leave>', lambda event :destroy_Tooltip(self))

		# Definition of list of cell lines in listbox once tissue has been selected
		def SpecCellLineUpdate(event):
			if self.selectTissue_combobox.get() in listCCLE:
				self.SpecCellLine_entry.configure(state=NORMAL)
				self.SpecCellLine_list.configure(state=NORMAL)
				self.SpecCellLine_entry_label.configure(state=NORMAL)
				self.SpecCellLine_list_label.configure(state=NORMAL)
				self.SpecCellLine_list.delete(0, 'end')
				for cellLine in cellLineCCLEdico[self.selectTissue_combobox.get()]:
					self.SpecCellLine_list.insert(END,cellLine)

		# Modificationn of list of cell line in listbox according to submitted value
		def SpecCellLine_keyrelease(event):
			value = event.widget.get()
			value = value.strip().lower()
			self.tissue = self.selectTissue_combobox.get()
			if value == '':
				data = cellLineCCLEdico[self.tissue]
			else:
				data = []
				for item in cellLineCCLEdico[self.tissue]:
					if item.lower().startswith(value):
						data.append(item)   
			# update data in listbox
			SpecCellLine_listbox_update(data)

		# Modification of list of cell lines in listbox according to submitted value in entrybox
		def SpecCellLine_listbox_update(data):
			# delete previous data
			self.SpecCellLine_list.delete(0, 'end')
			if data:
				# sorting data 
				data = sorted(data, key=str.lower)
				# put new data
				for item in data:
					self.SpecCellLine_list.insert('end', item)
			else:
				self.SpecCellLine_list.insert('end','Not found')

		# Activation of check button Demote proteins once selected cell line has been validated
		def SpecCellLineDemote(event):
			self.tissue = self.selectTissue_combobox.get()
			if self.SpecCellLine_list.get(0) in cellLineCCLEdico[self.tissue]:
				self.demote_option.configure(state=NORMAL)


		# Entry widget to submit cell line
		self.SpecCellLine_entry=Entry(self,state=DISABLED,font=("Helvetica"))
		self.SpecCellLine_entry.grid(row=16,column=2)
		self.SpecCellLine_entry.bind('<KeyRelease>', SpecCellLine_keyrelease)
		self.SpecCellLine_entry_label=Label(self,text='Type to help to find a cell line:',font=("Helvetica"),state=DISABLED)
		self.SpecCellLine_entry_label.grid(row=16,column=0)

		self.SpecCellDefilB = Scrollbar(self, orient='vertical')
		self.SpecCellDefilB.grid(row=17, column=3, sticky='ns' )

		# Listbox widget to select cell line
		self.SpecCellLine_list = Listbox(self,font=("Helvetica"),selectmode=SINGLE,yscrollcommand=self.SpecCellDefilB.set,exportselection=0)
		self.SpecCellLine_list.grid(row=17,column=2)
		self.SpecCellLine_list.bind('<<ListboxSelect>>', SpecCellLineDemote)
		self.SpecCellLine_list_label=Label(self,text='Choose a cell line :',font=("Helvetica"),state=DISABLED)
		self.SpecCellLine_list_label.grid(row=17,column=0)
		self.SpecCellLine_list.config(width=20,height=2,state=DISABLED)

		# Combobox to select tissue
		self.selectTissue_combobox = ttk.Combobox(self,state=DISABLED,textvariable = self.selected_tissue)
		self.selectTissue_combobox['values'] = listCCLE
		self.selectTissue_combobox.bind('<<ComboboxSelected>>', SpecCellLineUpdate)
		self.selectTissue_combobox.grid(row=15,column=2)

		self.selectTissue_label=Label(self,text='Choose a tissue :',font=("Helvetica"),state=DISABLED)
		self.selectTissue_label.grid(row=15,column=0)

		# Check button to demote or not proteins
		self.demote_option= Checkbutton(self, variable=self.DemoteCCLEOption,text="Demote unidentified proteins from CCLE", font=("Helvetica"),state=DISABLED)
		self.demote_option.grid(row=18, column=0)


		self.add_overflow = Checkbutton(self, variable=self.OverflowOption,text="Add overflow (in % of the shortest path length) :", font=("Helvetica"),command=self.Able_scaleOverflow)
		self.add_overflow.grid(row=19, column=0)

		self.label_docOverflow=Label(self,text='?',font=("Helvetica"))
		self.label_docOverflow.grid(row=19, column=1)
		self.label_docOverflow.bind('<Enter>', lambda event, text='The "overflow option" will display extra networks\n with alternative paths with a percentage of shortest path extra length.': add_Tooltip(self.label_docOverflow,text))
		self.label_docOverflow.bind('<Leave>', lambda event :destroy_Tooltip(self))

		self.overflow=Scale(self,state=DISABLED,orient=HORIZONTAL,width=30,length=300,variable=self.overflowValue)
		self.overflow.grid(row=19,column=2)

		#Network extraction
		self.ExtractionLabel=Label(self,text="Network extraction, enter list of :",font=("Helvetica"))
		self.ExtractionLabel.grid(row=20,column=0,columnspan=1,pady=5)

		self.label_docExtraction=Label(self,text='?',font=("Helvetica", ))
		self.label_docExtraction.grid(row=20, column=0, sticky='e')
		self.label_docExtraction.bind('<Enter>', lambda event, text='The "extraction option" will display the sub-networks focusing on a subset of targets. \nThese targets can be selected based on GO term associated categories \nand or a specified subset of targets.': add_Tooltip(self.label_docExtraction,text))
		self.label_docExtraction.bind('<Leave>', lambda event :destroy_Tooltip(self))



		self.targets=Checkbutton(self,text="Subset of target(s) :",command=self.Able_target_extraction,variable=self.target_extract_option,font=("Helvetica"))
		self.targets.grid(row=21,column=0)  

		self.target_list_extraction=Entry(self,state=DISABLED,font=("Helvetica"),textvariable=self.subset_targets)
		self.target_list_extraction.grid(row=21,column=2)

		self.Go_extract=Checkbutton(self,text="GO terms associated categories:",command=self.Able_GOterm_extraction,variable=self.GOterm_extract_option,font=("Helvetica"))
		self.Go_extract.grid(row=22,column=0)  

		self.yDefilCat = Scrollbar(self, orient='vertical')
		self.yDefilCat.grid(row=22, column=3, sticky='ns' )

		self.GO_list_extraction=Listbox(self,selectmode=MULTIPLE,font=("Helvetica"),yscrollcommand=self.yDefilCat.set,exportselection=0)
		self.GO_list_extraction.grid(row=22,column=2)

		self.yDefilCat['command'] = self.GO_list_extraction.yview
		for item in ["cell adhesion and motility","cell growth and death","cell transport and metabolism", "immune system and inflammation","cell differentiation"]:
			self.GO_list_extraction.insert(END, item)
		self.GO_list_extraction.config(width=0,height=2,state=DISABLED)

		
		#submit button
		self.submit=Button(self, text="Submit", font=("Helvetica"),command=self.cliquerSubmit,)
		self.submit.grid(row=23,column=0)

	#Function to select file
	def cliquerFile(self):
		self.filename =  filedialog.askopenfilename(initialdir = "/HOME",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
		self.labelFile.configure(text=self.filename)
	#function to select folder	
	def cliquerFolder(self):
		self.work_folder =  filedialog.askdirectory(initialdir = "/HOME",title = "Select folder")   
		self.labelFolder.configure(text=self.work_folder)	

	# able or disable widget if option is selected or not
	def Able_subst(self):
		if self.addSubstOption.get()==1:
			self.directSubs.configure(state=NORMAL) 
		else:
			self.directSubs.configure(state=DISABLED)

	def Able_KinPhos(self):
		if self.addKinPhos.get()==1:
			self.KinPhos_list.configure(state=NORMAL) 
		else:
			self.KinPhos_list.configure(state=DISABLED)

	def Able_Go_term(self):
		if self.addSpecMod.get()==1:
			self.SpecModOption_list.configure(state=NORMAL) 
		else:
			self.SpecModOption_list.configure(state=DISABLED)
	
	def Able_cellLine(self):
		if self.addGlobal.get()==1:
			self.selectTissue_combobox.configure(state=NORMAL)
			self.selectTissue_label.configure(state=NORMAL)
		else:
			self.selectTissue_combobox.configure(state=DISABLED)
			self.selectTissue_label.configure(state=DISABLED)
			self.SpecCellLine_entry.configure(state=DISABLED)
			self.SpecCellLine_entry_label.configure(state=DISABLED)
			self.SpecCellLine_list.configure(state=DISABLED)
			self.SpecCellLine_list_label.configure(state=DISABLED)
			self.demote_option.configure(state=DISABLED)
			self.selectTissue_combobox.selection_clear()
			self.SpecCellLine_list.selection_clear(0)
			self.demote_option.deselect()
			self.selected_tissue.get()==""
			self.Cell_line=[]
			self.DemoteCCLEOption.get()==0

	def Able_scaleOverflow(self):
		if self.OverflowOption.get()==1:
			self.overflow.configure(state=NORMAL) 
		else:
			self.overflow.configure(state=DISABLED)

	def Able_target_extraction(self):
		if self.target_extract_option.get()==1:
			self.target_list_extraction.configure(state=NORMAL) 
		else:
			self.target_list_extraction.configure(state=DISABLED)
	
	def Able_GOterm_extraction(self):
		if self.GOterm_extract_option.get()==1:
			self.GO_list_extraction.configure(state=NORMAL) 
		else:
			self.GO_list_extraction.configure(state=DISABLED)   
	###############################################################################""



	#implement check function that return an error if formulare isn't complete
	def cliquerSubmit(self):
		self.KinPhos_mod=self.KinPhos_list.curselection()
		self.CatExtraction=self.GO_list_extraction.curselection()
		self.Cell_line=self.SpecCellLine_list.curselection()
		if self.checkFormular() ==1:
			fenetre.destroy()
		else:
			showwarning('Missing values ', 'Some information must be filled', )


	# Verify the formular 
	def checkFormular(self):
		self.check=1
		if self.filename==None:
			self.inputFile_label.config(fg='red')
			self.check=0
		else:
			self.inputFile_label.config(fg='black') 

		if self.work_folder==None:
			self.inputFolder_label.config(fg='red')
			self.check=0
		else:
			self.inputFolder_label.config(fg='black') 

		if self.databaseKEGG.get()==0 and self.databasePC.get()==0:
			self.selectionDatabase_label.config(fg='red')
			self.check=0
		else:
			self.selectionDatabase_label.config(fg='black')

		if self.selection.get()!=1 and self.selection.get()!=2:
			self.selectionMode_label.config(fg='red')
			self.check=0
		else:
			self.selectionMode_label.config(fg='black')

		if self.source.get()=="":
			self.SourceLabel.config(fg="red")
			self.check=0
		else:
			self.SourceLabel.config(fg="black")

		if self.addSubstOption.get()==1 and self.subs.get()=="":
			self.add_dirrectSubst.config(fg="red")
			self.check=0
		else:
			self.add_dirrectSubst.config(fg="black")

		if self.addKinPhos.get()==1 and not self.KinPhos_mod:
			self.KinPhos.config(fg="red")
			self.check=0
		else:
			self.KinPhos.config(fg="black")

		if self.addSpecMod.get()==1 and self.Mod_list.get()=="":
			self.SpecModOption.config(fg='red')
			self.check=0
		else:
			self.SpecModOption.config(fg='black')

		if self.addGlobal.get()==1 and self.selected_tissue.get()=="":
			self.addGlobalProteomic.config(fg='red')
			self.selectTissue_label.config(fg='red')
			self.check=0
		else:
			self.addGlobalProteomic.config(fg='black')

		if self.addGlobal.get()==1 and not self.Cell_line:
			self.addGlobalProteomic.config(fg='red')
			self.SpecCellLine_list_label.config(fg='red')
			self.check=0
		else:
			self.addGlobalProteomic.config(fg='black')

		if self.addGlobal.get()==1 and self.Cell_line and self.selected_tissue.get()!="":
			self.cell_line=self.SpecCellLine_list.get(self.Cell_line[0])+"_"+self.selected_tissue.get().upper()

		if self.OverflowOption.get()==1 and self.overflowValue.get()==0:
			self.add_overflow.config(fg='red')
			self.check=0
		else:
			self.add_overflow.config(fg='black')

		if self.target_extract_option.get()==1 and self.subset_targets.get()=="":
			self.targets.config(fg='red')
			self.check=0
		else:
			self.targets.config(fg='black')

		if self.GOterm_extract_option.get()==1 and not self.CatExtraction:
			self.Go_extract.config(fg='red')
			self.check=0
		else:
			self.Go_extract.config(fg='black')

		return self.check 

fenetre = Tk()
fenetre.title("Phos2Net")
#fenetre.config(background='#41B77F')
#img = fenetre.Image("photo", file="LogoPhos2Net.png")
img = Image("photo", file="LogoPhos2Net.gif")
fenetre.tk.call('wm','iconphoto',fenetre._w, img)
interface = Interface(fenetre)



interface.mainloop()
'''
if __name__ == "__main__":
	fenetre = Tk()
	interface = Interface(fenetre)
	interface.mainloop()
	interface.destroy()
'''
