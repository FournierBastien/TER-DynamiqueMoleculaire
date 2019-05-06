
x,v,name=[5, 0,5,9],[7,4,1,5],[8 ,2,0,4]

tab1=[]
tab2=[]
tab3=[]
def tri_ins(t1,t2,t3):
	for k in range(1,len(t1)):
		temp1=t1[k]
		temp2=t2[k]
		temp3=t3[k]
		j=k
		while j>0 and temp1<t1[j-1]:
			t1[j]=t1[j-1]
			t2[j]=t2[j-1]
			t3[j]=t3[j-1]
			j-=1
		t1[j]=temp1
		t2[j]=temp2
		t3[j]=temp3

def tri(t1,t2,t3):
	print (t1,t2,t3)
	n=len(t1)
	if n<2:
		t=(t1,t2,t3)
		#print ("n<2")
		print (type(t))
		return t
	else:
		m=n//2
		print(t1[:m],t2[:m],t3[:m],"  ",t1[m:],t2[m:],t3[m:])
		#print("type t1=",type(t1))
		return fusion(tri(t1[:m],t2[:m],t3[:m]),tri(t1[m:],t2[m:],t3[m:]))

def fusion(l1,l2):
	t1=[]
	t2=[]
	t3=[]
		
	print ("tuple l1,l2=",l1,"	",l2)
	if l1 is not None:
		"""print("fusion")
		print("f",type(l1))
		print (l1)
		print (l2)"""
		global tab1
		global tab2
		global tab3
		#print ("l1[0]=",l1[0])
		#print ("l2[0]=",l2[0])
		if l1[0]==[]:
			tab1+=l2[0]
			tab2+=l2[1]
			tab3+=l2[2]
			tri_ins(tab1,tab2,tab3)
			return l2[0],l2[1],l2[2]
		elif l2[0]==[]:
			print("l2[0] Vide Avnat  tab1=",tab1,"tab2=",tab2,"tab3=",tab3)
			tab1+=l1[0]
			tab2+=l1[1]
			tab3+=l1[2]
			tri_ins(tab1,tab2,tab3)
			print("l2[0] Vide tab1=",tab1,"tab2=",tab2,"tab3=",tab3)
			return l1[0],l1[1],l1[2]

		elif l1[0][0]<l2[0][0]:
			#print(" IF l1[0][0]<l2[0][0]",l1[0][0],"	",l2[0][0])

			#print ("tab1=",tab1,"l1[0]=",l1[0])
			tab1=l1[0]+tab1
			#print ("tab1 apres fusuion=",tab1)
			tab2=l1[1]+tab2
			tab3=l1[2]+tab3
			tri_ins(tab1,tab2,tab3)
			print ("l1[0][0]<l2[0][0] tab1=",tab1)
			fusion((l1[0][1:],l1[1][1:],l1[2][1:]),(l2[0],l2[1],l2[2]))
		else:
			#print ("Else tab1=",tab1,"l1[0]=",l1[0])
			tab1=l2[0]+tab1
			tab2=l2[1]+tab2
			tab3=l2[2]+tab3
			tri_ins(tab1,tab2,tab3)
			print ("l1[0][0]>l2[0][0] tab1=",tab1,"tab2=",tab2,"tab3=",tab3)
			fusion((l1[0],l1[1],l1[2]),(l2[0][1:],l2[1][1:],l2[2][1:]))

def fusionAux(t1,t2):
	if t1==[]:
		return t2
	elif t2==[]:
		return t1
	elif t1[0]<t2[0]:
		return [t1[0]]+fusion(t1[1:],t2)
	else:
		return [t2[0]]+fusion(t1,t2[1:])
tri(x,v,name)

print(tab1,tab2,tab3)

#[5, 2, 0, 4],[3, 1, 6, 2],[9, 7, 1, 8]

"""
x,v,name=[5, 2, 0, 4],[3, 1, 6, 2],[9, 7, 1, 8]

tab1=[]
tab2=[]
tab3=[]

def tri(t1,t2,t3):
	print (t1,t2,t3)
	n=len(t1)
	if n<2:
		t=(t1,t2,t3)
		#print ("n<2")
		print (type(t))
		return t
	else:
		m=n//2
		print(t1[:m],t2[:m],t3[:m],"  ",t1[m:],t2[m:],t3[m:])
		#print("type t1=",type(t1))
		return fusion(tri(t1[:m],t2[:m],t3[:m]),tri(t1[m:],t2[m:],t3[m:]))

def fusion(l1,l2):
	print ("tuple l1,l2=",l1,"	",l2)
	if l1 is not None:
		print("fusion")
		print("f",type(l1))
		print (l1)
		print (l2)
		global tab1
		global tab2
		global tab3
		print ("l1[0]=",l1[0])
		if l1[0]==[]:
			tab1+=l2[0]
			tab2+=l2[1]
			tab3+=l2[2]
			return l2[0],l2[1],l2[2]
		elif l2[0]==[]:
			print("l2[0] Vide Avnat  tab1=",tab1,"tab2=",tab2,"tab3=",tab3)
			tab1+=l1[0]
			tab2+=l1[1]
			tab3+=l1[2]
			print("l2[0] Vide tab1=",tab1,"tab2=",tab2,"tab3=",tab3)
			return l1[0],l1[1],l1[2]

		elif l1[0][0]<l2[0][0]:
			#print(" IF l1[0][0]<l2[0][0]",l1[0][0],"	",l2[0][0])

			#print ("tab1=",tab1,"l1[0]=",l1[0])
			tab1=l1[0]+tab1
			#print ("tab1 apres fusuion=",tab1)
			tab2=l1[1]+tab2
			tab3=l1[2]+tab3
			print ("l1[0][0]<l2[0][0] tab1=",tab1)
			return fusion((l1[0][1:],l1[1][1:],l1[2][1:]),(l2[0],l2[1],l2[2]))

		else:
			#print ("Else tab1=",tab1,"l1[0]=",l1[0])
			tab1=l2[0]+tab1
			tab2=l2[1]+tab2
			tab3=l2[2]+tab3
			print ("l1[0][0]>l2[0][0] tab1=",tab1,"tab2=",tab2,"tab3=",tab3)
			return fusion((l1[0],l1[1],l1[2]),(l2[0][1:],l2[1][1:],l2[2][1:]))


tri(x,v,name)

print(tab1,tab2,tab3)"""
