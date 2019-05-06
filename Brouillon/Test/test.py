# coding: utf-8


import sys
import os
import matplotlib.pyplot as plt
from math import sqrt, exp
import numpy
import time
import random
import pf

# Specification du dossier de travail
#os.chdir("/media/mcd/MCD/M1Informatique/S2/TER/DynamiqueMoleculaire/Code_&_Resultat/Test")

class Univexp():

	fichier="univexp01"
	tecr=0.
	tsor=0.
	n=100

	def __init__(self):

		self.init()
		self.f = open(self.fichier, "w")
		self.wbande()

		self.tecr = self.tecr + self.dti
		self.tsor = self.tecr + self.dtsor
		#print("tecr=",self.tecr)
		while abs(self.tecr -self.tstop) > self.dtsor / 2.:
			#print("abs(tecr-tstop)=",abs(self.tecr -self.tstop),"dtsor/2.=",self.dtsor / 2.)
			while abs(self.tecr -self.tsor) > self.dti / 2.:
				self.avance()
				self.ordonne()
				self.tecr = self.tecr + self.dti
				#plt.scatter(self.x, self.v)
				#plt.pause(0.01)
			self.wbande()
			self.tsor = self.tsor + self.dtsor

		self.f.close()
		print ("fin de programme.")
        
    #************************************************************************************************

	def avance(self):
		r2 = -sqrt(2.)
		tier = 1./3.
		#AA=BB=E=0.
		
		#calcul de l'acceleration
		a =[0] * self.m# numpy.zeros((self.m), dtype='d')
		self.eav = self.epolar + 0.5 * self.n
		for i in range(self.ifirst, self.ilast):
			#eap = eav - self.gravplas * self.ma[i]
			#a[i] = 0.5 * (eav + eap) / self.mi[i]
			self.eap = self.eav - 1
			a[i] = 0.5 * (self.eav + self.eap)
			self.eav = self.eap

		#les particules avancent de dti
		#for i in self.x:
			#print("x0=",i)
		#print("E0=",E,"AA0=",AA,"BB0=",BB,"dti0=",self.dti,"r20=",r2)
		for i in range(self.ifirst, self.ilast):
			#plt.scatter(self.x, self.v)
			#plt.pause(0.01)
			E = a[i]
			AA = (self.x[i] + E + r2 * self.v[i]) * exp(r2 * self.dti)
			BB = 2. * (self.x[i] + E - self.v[i]/r2) * exp(-self.dti/r2)
			self.x[i] = ((AA + BB) * tier) - E
			self.v[i] = (AA * r2 - BB/r2) * tier
			#plt.scatter(self.x, self.v)
			#plt.pause(0.01)
		    
			if self.x[i] > self.x[self.m-1]:
				#self.epolar self.epolar + self.gravplas * self.ma[i]
				self.epolar = self.epolar + 1
				self.x[i] = self.x[0] + self.x[i] - self.x[self.m-1]

			if self.x[i] < self.x[0] :
				#self.epolar self.epolar - self.gravplas * self.ma[i]
				self.epolar = self.epolar - 1
				self.x[i] = self.x[self.m-1] + self.x[i] - self.x[0]

                


    #************************************************************************************************

	def ordonne(self):
		"""pf.univexp.ordonne
		#print('ordonne ok')"""
		
		i=j=k=np=0

		xp=vp=mip=map=0.

		for i in range(self.ifirst+1, self.ilast):
			j = i 
			xp=self.x[i]
			vp=self.v[i]
			#mip=self.mi[i]
			#map=self.ma[i]
			np=self.name[i]

			while j!=0 and self.x[i]<self.x[j-1]:
				j-=1
		
			for k in range(i,j,-1):
				self.x[k]=self.x[k-1]
				self.v[k]=self.v[k-1]
				#self.mi[k]=self.mi[k-1]
				#self.ma[k]=self.ma[k-1]
				self.name[k]=self.name[k-1]

			self.x[j]=xp
			self.v[j]=vp
			#self.mi[j]=mip
			#self.ma[j]=map
			self.name[j]=np
    
   #************************************************************************************************

	def wbande(self):
		print("enregistrement a tecr = {:7.3f}".format(self.tecr))

		self.f.write("#   n={:<5d} m={:<5d}  ifirst={:<5d} ilast={:<5d} tecr={:<7.3f} \n".format(self.m, self.n, self.ifirst, self.ilast, self.tecr))
        
		for i in range(self.ifirst, self.ilast):
			self.f.write("   {:>2.9f}    \t {:>4.9f}     \t {:>5}\n".format(self.x[i],self.v[i],self.name[i]))

		self.f.write("\n\n")


    #************************************************************************************************

	def init(self):
		self.m = self.n + 2
		self.ifirst = 1
		self.ilast=self.n + 1
		self.dti = 0.001
		self.dtsor = 1.
		self.tstop = 15.
		self.gravplas=1.
		self.pvit = 100.

		#allocation pour les tableaux
		"""self.x = numpy.zeros((self.m), dtype='d')
		self.v = numpy.zeros((self.m), dtype='d')
		self.mi = numpy.zeros((self.m), dtype='d')
		self.ma = numpy.zeros((self.m), dtype='d')
		self.name = numpy.zeros((self.m), dtype='l')"""
		self.x = [0] * self.m
		self.v = [0] * self.m
		self.mi = [0] * self.m
		self.ma = [0] * self.m
		self.name = [0] * self.m


		print("simulation :",self.fichier)
		print("{:19.15f} : pas de temps ".format(self.dti))
		print("n={:8d} pvit={:7.1f}".format(self.n, self.pvit))
		print("tstop={:7.1f} dtsor ={:7.1f}".format(self.tstop, self.dtsor))

		# initialisation des "particules-mur"
		self.x[0] = -0.5*self.n
		self.v[0] = 0
		self.mi[0] = 0
		self.ma[0] = 0
		self.name[0] = 0

		self.x[self.m-1] = 0.5*self.n
		self.v[self.m-1] = 0
		self.mi[self.m-1] = 0
		self.ma[self.m-1] = 0
		self.name[self.m-1] = 0

		# initialisation des particules
		
		with open("univexp01d.d",'r') as fich:
			fich.readline()
			i=1
			while i < self.ilast:
				ligne = fich.readline()
				self.v[i] = float(ligne.split()[1])
				#self.x[i] = float(ligne.split()[0]) 
				i+=1

		for i in range(self.ifirst, self.ilast):
			self.name[i]=i#Pour la camparaison des deux fichiers sinon i-1
			self.x[i]=self.x[0]+0.5+i-self.ifirst
			#self.v[i]=2.0*self.pvit*(0.5-random.random())
			self.mi[i]=1.0
			self.ma[i]=1.0
		self.epolar=0


			
		#calage du barycentre à zéro avec une vitesse moyenne  nulle
		"""for v in self.v:
			print("v=",v)
		s=0
		for i in range(self.ifirst,self.ilast):
			s+=self.v[i]
			print("s=",s)"""
		vmoy=sum(self.v[self.ifirst:self.ilast])
		print(("vmoyen = {:2.8f}").format(vmoy))


def main() :
	# Debut du decompte du temps
	start_time = time.time()
	Univexp()
	#plt.show()
	# Affichage du temps d execution
	print("Temps d execution : %s secondes" % (time.time() - start_time))

main()
"""import cProfile
import re
cProfile.run('main()')"""

#cProfile.run('re.compile("foo|bar")', 'restats')
exit()
