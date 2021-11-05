from tkinter import filedialog
from tkinter import *
from tkinter.messagebox import *
from tkinter import ttk
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
		self.OverflowOption=IntVar()
		self.overflowValue =IntVar()
		self.target_extract_option=IntVar()
		self.GOterm_extract_option=IntVar() 
		self.source=StringVar()   
		self.subs=StringVar()
		self.Mod_list=StringVar()
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

		self.Go_extract=Checkbutton(self,text="GO terms associated categories    :",command=self.Able_GOterm_extraction,variable=self.GOterm_extract_option,font=("Helvetica"))
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