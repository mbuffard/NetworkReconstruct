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
		self.OverflowOption=IntVar()
		self.overflowValue =IntVar()
		self.target_extract_option=IntVar()
		self.GOterm_extract_option=IntVar() 
		self.source=StringVar()   
		self.subs=StringVar()
		self.Mod_list=StringVar()
		self.selected_tissue=StringVar()
		self.Cell_line=StringVar()
		self.KinPhos_mod=[]
		self.subset_targets=StringVar()
		self.selection=IntVar()
		self.check=0
		self.CatExtraction=[]
		self.target_extract_option=IntVar()
		self.GOterm_extract_option=IntVar()


		Frame.__init__(self, fenetre, width=768, height=576, **kwargs)
		self.pack(fill=BOTH)
		#todo nicer windows
		#self.config(background='#41B77F')

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
		self.WeightLabel.grid(row=11,column=0,columnspan=2,pady=5)

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
		cellLineCCLEdico = {'prostate': ['22RV1_PROSTATE', 'DU145_PROSTATE', 'LNCAPCLONEFGC_PROSTATE', 'PC3_PROSTATE', 'VCAP_PROSTATE'], 'haematopoietic and lymphoid tissue': ['697_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'CMK_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'EOL1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'F36P_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'HDMYZ_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'HEL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'JEKO1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'JM1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'KARPAS299_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'KARPAS422_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'KASUMI1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'KASUMI2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'KMS11_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'KMS12BM_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'KMS27_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'L428_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'MOLM13_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'MOLM16_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'MONOMAC1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'MONOMAC6_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'NALM6_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'NCIH929_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'NUDHL1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'OCIAML5_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'OCILY3_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'OPM2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'RCHACV_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'REC1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'REH_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'RL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'RPMI8226_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'SEM_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'SUDHL4_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'SUDHL6_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'TF1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'U937_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'], 'kidney': ['769P_KIDNEY', '786O_KIDNEY', 'A498_KIDNEY', 'A704_KIDNEY', 'CAKI1_KIDNEY', 'CAKI2_KIDNEY', 'KMRC1_KIDNEY', 'KMRC20_KIDNEY', 'OSRC2_KIDNEY', 'TUHR4TKB_KIDNEY', 'VMRCRCW_KIDNEY'], 'thyroid': ['8305C_THYROID', '8505C_THYROID'], 'skin': ['A101D_SKIN', 'A2058_SKIN', 'A375_SKIN', 'C32_SKIN', 'COLO679_SKIN', 'COLO741_SKIN', 'COLO829_SKIN', 'HS294T_SKIN', 'HS695T_SKIN', 'HS944T_SKIN', 'IGR1_SKIN', 'IGR37_SKIN', 'IGR39_SKIN', 'IPC298_SKIN', 'K029AX_SKIN', 'LOXIMVI_SKIN', 'MELJUSO_SKIN', 'MEWO_SKIN', 'RPMI7951_SKIN', 'RVH421_SKIN', 'SH4_SKIN', 'SKMEL2_SKIN', 'SKMEL28_SKIN', 'SKMEL3_SKIN', 'SKMEL30_SKIN', 'SKMEL5_SKIN', 'UACC257_SKIN', 'UACC62_SKIN', 'WM115_SKIN', 'WM1799_SKIN', 'WM2664_SKIN', 'WM793_SKIN', 'WM88_SKIN'], 'central nervous system': ['A172_CENTRAL_NERVOUS_SYSTEM', 'DAOY_CENTRAL_NERVOUS_SYSTEM', 'GAMG_CENTRAL_NERVOUS_SYSTEM', 'GB1_CENTRAL_NERVOUS_SYSTEM', 'KNS42_CENTRAL_NERVOUS_SYSTEM', 'KNS81_CENTRAL_NERVOUS_SYSTEM', 'LN18_CENTRAL_NERVOUS_SYSTEM', 'LN229_CENTRAL_NERVOUS_SYSTEM', 'SF295_CENTRAL_NERVOUS_SYSTEM', 'SNB19_CENTRAL_NERVOUS_SYSTEM', 'SNU1105_CENTRAL_NERVOUS_SYSTEM', 'SW1783_CENTRAL_NERVOUS_SYSTEM', 'U118MG_CENTRAL_NERVOUS_SYSTEM', 'U87MG_CENTRAL_NERVOUS_SYSTEM'], 'soft tissue': ['A204_SOFT_TISSUE', 'G401_SOFT_TISSUE', 'G402_SOFT_TISSUE', 'HT1080_SOFT_TISSUE', 'KYM1_SOFT_TISSUE', 'RD_SOFT_TISSUE', 'RH30_SOFT_TISSUE', 'RH41_SOFT_TISSUE'], 'ovary': ['A2780_OVARY', 'CAOV3_OVARY', 'COV362_OVARY', 'FUOV1_OVARY', 'HEYA8_OVARY', 'IGROV1_OVARY', 'JHOS2_OVARY', 'KURAMOCHI_OVARY', 'NIHOVCAR3_OVARY', 'OV56_OVARY', 'OV90_OVARY', 'OVCAR4_OVARY', 'OVCAR8_OVARY', 'OVSAHO_OVARY', 'RMUGS_OVARY', 'SNU119_OVARY', 'TYKNU_OVARY'], 'lung': ['A549_LUNG', 'ABC1_LUNG', 'CALU1_LUNG', 'CALU6_LUNG', 'CHAGOK1_LUNG', 'CORL105_LUNG', 'CORL23_LUNG', 'CORL47_LUNG', 'CORL88_LUNG', 'DMS114_LUNG', 'DMS273_LUNG', 'DV90_LUNG', 'EBC1_LUNG', 'HCC15_LUNG', 'HCC1833_LUNG', 'HCC44_LUNG', 'HCC827_LUNG', 'HCC95_LUNG', 'IALM_LUNG', 'LCLC103H_LUNG', 'LK2_LUNG', 'LU65_LUNG', 'LUDLU1_LUNG', 'LXF289_LUNG', 'NCIH1048_LUNG', 'NCIH1155_LUNG', 'NCIH1299_LUNG', 'NCIH1355_LUNG', 'NCIH1435_LUNG', 'NCIH1437_LUNG', 'NCIH146_LUNG', 'NCIH1568_LUNG', 'NCIH1573_LUNG', 'NCIH1581_LUNG', 'NCIH1650_LUNG', 'NCIH1666_LUNG', 'NCIH1693_LUNG', 'NCIH1703_LUNG', 'NCIH1792_LUNG', 'NCIH1793_LUNG', 'NCIH1944_LUNG', 'NCIH196_LUNG', 'NCIH1975_LUNG', 'NCIH2009_LUNG', 'NCIH2030_LUNG', 'NCIH2066_LUNG', 'NCIH2110_LUNG', 'NCIH2122_LUNG', 'NCIH2126_LUNG', 'NCIH2170_LUNG', 'NCIH2172_LUNG', 'NCIH2228_LUNG', 'NCIH226_LUNG', 'NCIH2286_LUNG', 'NCIH2291_LUNG', 'NCIH23_LUNG', 'NCIH292_LUNG', 'NCIH3255_LUNG', 'NCIH358_LUNG', 'NCIH441_LUNG', 'NCIH446_LUNG', 'NCIH460_LUNG', 'NCIH520_LUNG', 'NCIH522_LUNG', 'NCIH647_LUNG', 'NCIH650_LUNG', 'NCIH661_LUNG', 'NCIH838_LUNG', 'PC14_LUNG', 'RERFLCMS_LUNG', 'RERFLCSQ1_LUNG', 'SBC5_LUNG', 'SHP77_LUNG', 'SKLU1_LUNG', 'SQ1_LUNG', 'SW1271_LUNG', 'SW1573_LUNG'], 'bone': ['A673_BONE', 'SAOS2_BONE', 'SJSA1_BONE', 'SKES1_BONE', 'SKNMC_BONE', 'TC71_BONE', 'U2OS_BONE'], 'stomach': ['AGS_STOMACH', 'HGC27_STOMACH', 'HUG1N_STOMACH', 'IM95_STOMACH', 'KATOIII_STOMACH', 'LMSU_STOMACH', 'MKN1_STOMACH', 'MKN45_STOMACH', 'MKN7_STOMACH', 'NCIN87_STOMACH', 'NUGC3_STOMACH', 'OCUM1_STOMACH', 'SNU1_STOMACH', 'SNU719_STOMACH'], 'pancreas': ['ASPC1_PANCREAS', 'BXPC3_PANCREAS', 'CFPAC1_PANCREAS', 'DANG_PANCREAS', 'HUPT3_PANCREAS', 'HUPT4_PANCREAS', 'KP2_PANCREAS', 'KP4_PANCREAS', 'L33_PANCREAS', 'MIAPACA2_PANCREAS', 'PANC0203_PANCREAS', 'PANC0403_PANCREAS', 'PANC1_PANCREAS', 'PATU8988T_PANCREAS', 'PL45_PANCREAS', 'QGP1_PANCREAS', 'SU8686_PANCREAS', 'SUIT2_PANCREAS', 'SW1990_PANCREAS', 'TCCPAN2_PANCREAS'], 'breast': ['AU565_BREAST', 'BT20_BREAST', 'BT549_BREAST', 'CAL120_BREAST', 'CAL51_BREAST', 'CAL851_BREAST', 'CAMA1_BREAST', 'EFM19_BREAST', 'EFM192A_BREAST', 'HCC1143_BREAST', 'HCC1187_BREAST', 'HCC1395_BREAST', 'HCC1500_BREAST', 'HCC1806_BREAST', 'HCC1937_BREAST', 'HCC1954_BREAST', 'HCC2218_BREAST', 'HCC38_BREAST', 'HCC70_BREAST', 'HDQP1_BREAST', 'JIMT1_BREAST', 'KPL1_BREAST', 'MCF7_BREAST', 'MDAMB157_BREAST', 'MDAMB231_BREAST', 'MDAMB436_BREAST', 'MDAMB453_BREAST', 'MDAMB468_BREAST', 'T47D_BREAST', 'ZR751_BREAST'], 'upper aerodigestive tract': ['BICR22_UPPER_AERODIGESTIVE_TRACT', 'BICR6_UPPER_AERODIGESTIVE_TRACT', 'CAL27_UPPER_AERODIGESTIVE_TRACT', 'CAL33_UPPER_AERODIGESTIVE_TRACT', 'DETROIT562_UPPER_AERODIGESTIVE_TRACT', 'FADU_UPPER_AERODIGESTIVE_TRACT', 'HSC3_UPPER_AERODIGESTIVE_TRACT', 'HSC4_UPPER_AERODIGESTIVE_TRACT', 'PECAPJ34CLONEC12_UPPER_AERODIGESTIVE_TRACT', 'SCC25_UPPER_AERODIGESTIVE_TRACT'], 'large intestine': ['CCK81_LARGE_INTESTINE', 'CL34_LARGE_INTESTINE', 'COLO205_LARGE_INTESTINE', 'COLO320_LARGE_INTESTINE', 'COLO678_LARGE_INTESTINE', 'HCC56_LARGE_INTESTINE', 'HCT116_LARGE_INTESTINE', 'HCT15_LARGE_INTESTINE', 'HT115_LARGE_INTESTINE', 'HT29_LARGE_INTESTINE', 'HT55_LARGE_INTESTINE', 'LS180_LARGE_INTESTINE', 'LS411N_LARGE_INTESTINE', 'LS513_LARGE_INTESTINE', 'MDST8_LARGE_INTESTINE', 'NCIH716_LARGE_INTESTINE', 'NCIH747_LARGE_INTESTINE', 'RKO_LARGE_INTESTINE', 'SKCO1_LARGE_INTESTINE', 'SNU61_LARGE_INTESTINE', 'SNUC1_LARGE_INTESTINE', 'SNUC2A_LARGE_INTESTINE', 'SNUC5_LARGE_INTESTINE', 'SW1417_LARGE_INTESTINE', 'SW403_LARGE_INTESTINE', 'SW48_LARGE_INTESTINE', 'SW480_LARGE_INTESTINE', 'SW620_LARGE_INTESTINE', 'SW837_LARGE_INTESTINE', 'SW948_LARGE_INTESTINE'], 'endometrium': ['HEC108_ENDOMETRIUM', 'HEC1A_ENDOMETRIUM', 'HEC251_ENDOMETRIUM', 'HEC265_ENDOMETRIUM', 'HEC50B_ENDOMETRIUM', 'HEC59_ENDOMETRIUM', 'HEC6_ENDOMETRIUM', 'ISHIKAWAHERAKLIO02ER_ENDOMETRIUM', 'JHUEM2_ENDOMETRIUM', 'MFE280_ENDOMETRIUM', 'MFE296_ENDOMETRIUM', 'MFE319_ENDOMETRIUM', 'SNGM_ENDOMETRIUM', 'SNU685_ENDOMETRIUM'], 'liver': ['HEP3B217_LIVER', 'HEPG2_LIVER', 'HLF_LIVER', 'HUH1_LIVER', 'HUH6_LIVER', 'HUH7_LIVER', 'JHH1_LIVER', 'JHH4_LIVER', 'JHH5_LIVER', 'JHH6_LIVER', 'JHH7_LIVER', 'SKHEP1_LIVER', 'SNU423_LIVER', 'SNU449_LIVER'], 'urinary tract': ['HT1197_URINARY_TRACT', 'HT1376_URINARY_TRACT', 'J82_URINARY_TRACT', 'JMSU1_URINARY_TRACT', 'KU1919_URINARY_TRACT', 'RT112_URINARY_TRACT', 'RT4_URINARY_TRACT', 'T24_URINARY_TRACT', 'TCCSUP_URINARY_TRACT', 'UBLC1_URINARY_TRACT', 'UMUC3_URINARY_TRACT'], 'oesophagus': ['KYSE150_OESOPHAGUS', 'KYSE180_OESOPHAGUS', 'KYSE30_OESOPHAGUS', 'KYSE410_OESOPHAGUS', 'KYSE450_OESOPHAGUS', 'KYSE510_OESOPHAGUS', 'KYSE70_OESOPHAGUS', 'OE33_OESOPHAGUS', 'TE1_OESOPHAGUS', 'TE10_OESOPHAGUS', 'TE11_OESOPHAGUS', 'TE14_OESOPHAGUS', 'TE4_OESOPHAGUS', 'TE6_OESOPHAGUS'], 'pleura': ['MSTO211H_PLEURA', 'NCIH2052_PLEURA'], 'autonomic ganglia': ['SKNAS_AUTONOMIC_GANGLIA'], 'biliary tract': ['SNU1079_BILIARY_TRACT']}

		#option to add global proteomic data for a specified cell line
		self.addGlobalProteomic=Checkbutton(self,text="Add global proteomic data from CCLE",command=self.Able_cellLine,variable=self.addGlobal,font=("Helvetica"))
		self.addGlobalProteomic.grid(row=14,column=0)

		def SpecCellLineUpdate(event):
			if self.selectTissue_combobox.get() in listCCLE:
				self.SpecCellLine_combobox.configure(state=NORMAL)
				self.SpecCellLine_combobox['values'] = cellLineCCLEdico[self.selectTissue_combobox.get()]


		self.SpecCellLine_combobox= ttk.Combobox(self,state=DISABLED,textvariable = self.Cell_line)
		self.SpecCellLine_combobox.grid(row=14,column=3)

		self.selectTissue_combobox = ttk.Combobox(self,state=DISABLED,textvariable = self.selected_tissue)
		self.selectTissue_combobox['values'] = listCCLE
		self.selectTissue_combobox.bind('<<ComboboxSelected>>', SpecCellLineUpdate)
		self.selectTissue_combobox.grid(row=14,column=2)


		self.add_overflow = Checkbutton(self, variable=self.OverflowOption,text="Add overflow (in % of the shortest path length) :", font=("Helvetica"),command=self.Able_scaleOverflow)
		self.add_overflow.grid(row=15, column=0)

		self.label_docOverflow=Label(self,text='?',font=("Helvetica", ))
		self.label_docOverflow.grid(row=15, column=0, sticky='e')
		self.label_docOverflow.bind('<Enter>', lambda event, text='The "overflow option" will display extra networks\n with alternative paths with a percentage of shortest path extra length.': add_Tooltip(self.label_docOverflow,text))
		self.label_docOverflow.bind('<Leave>', lambda event :destroy_Tooltip(self))

		self.overflow=Scale(self,state=DISABLED,orient=HORIZONTAL,width=30,length=300,variable=self.overflowValue)
		self.overflow.grid(row=15,column=2)

		#Network extraction
		self.ExtractionLabel=Label(self,text="Network extraction, enter list of :",font=("Helvetica"))
		self.ExtractionLabel.grid(row=16,column=0,columnspan=2,pady=5)

		self.label_docExtraction=Label(self,text='?',font=("Helvetica", ))
		self.label_docExtraction.grid(row=16, column=0, sticky='e')
		self.label_docExtraction.bind('<Enter>', lambda event, text='The "extraction option" will display the sub-networks focusing on a subset of targets. \nThese targets can be selected based on GO term associated categories \nand or a specified subset of targets.': add_Tooltip(self.label_docExtraction,text))
		self.label_docExtraction.bind('<Leave>', lambda event :destroy_Tooltip(self))



		self.targets=Checkbutton(self,text="Subset of target(s) :",command=self.Able_target_extraction,variable=self.target_extract_option,font=("Helvetica"))
		self.targets.grid(row=17,column=0)  

		self.target_list_extraction=Entry(self,state=DISABLED,font=("Helvetica"),textvariable=self.subset_targets)
		self.target_list_extraction.grid(row=17,column=2)

		self.Go_extract=Checkbutton(self,text="GO terms associated categories:",command=self.Able_GOterm_extraction,variable=self.GOterm_extract_option,font=("Helvetica"))
		self.Go_extract.grid(row=18,column=0)  

		self.yDefilCat = Scrollbar(self, orient='vertical')
		self.yDefilCat.grid(row=18, column=3, sticky='ns' )

		self.GO_list_extraction=Listbox(self,selectmode=MULTIPLE,font=("Helvetica"),yscrollcommand=self.yDefilCat.set,exportselection=0)
		self.GO_list_extraction.grid(row=18,column=2)

		self.yDefilCat['command'] = self.GO_list_extraction.yview
		for item in ["cell adhesion and motility","cell growth and death","cell transport and metabolism", "immune system and inflammation","cell differentiation"]:
			self.GO_list_extraction.insert(END, item)
		self.GO_list_extraction.config(width=0,height=2,state=DISABLED)

		
		#submit button
		self.submit=Button(self, text="Submit", font=("Helvetica"),command=self.cliquerSubmit,)
		#
		self.submit.grid(row=19,column=0)

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
		else:
			self.selectTissue_combobox.configure(state=DISABLED)
			self.SpecCellLine_combobox.configure(state=DISABLED)

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
		print(self.selected_tissue.get())
		print(self.Cell_line.get())
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
			self.check=0
		else:
			self.addGlobalProteomic.config(fg='black')

		if self.addGlobal.get()==1 and self.Cell_line.get()=="":
			self.addGlobalProteomic.config(fg='red')
			self.check=0
		else:
			self.addGlobalProteomic.config(fg='black')

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