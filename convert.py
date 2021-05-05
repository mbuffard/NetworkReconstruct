import converter

f=open("liste5.txt","r")
g=open("liste5_clean.txt",'w')

for line in f:
	#print converter.handler.clean_uid(line.strip().split()[0])
	uid_clean=converter.handler.clean_uid(line.strip().split()[0])
	if uid_clean!=None:
		g.write(uid_clean+"\n")
	
 
f.close()




g.close()
