# coding: utf-8


import sys
import os
import matplotlib.pyplot as plt
from math import sqrt, exp
import numpy
import time
import random
from threading import *
lock = RLock()
# Specification du dossier de travail
#os.chdir("/media/mcd/MCD/M1Informatique/S2/TER/DynamiqueMoleculaire/Code_&_Resultat/Test")

class Univexp(Thread):

	fichier="univexp01"
	tecr=0.
	tsor=0.
	n=100

	def __init__(self):
		Thread.__init__(self)
		self.init()
		self.f = open(self.fichier, "w")
		self.wbande()

		self.tecr = self.tecr + self.dti
		self.tsor = self.tecr + self.dtsor

	def run(self):
		while True:
			if abs(self.tecr -self.tstop) > self.dtsor / 2.:
				if abs(self.tecr -self.tsor) > self.dti / 2.:
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
		a = numpy.zeros((self.m), dtype='d')
		self.eav = self.epolar + 0.5 * self.n
		for i in range(self.ifirst, self.ilast):
			self.eap = self.eav - 1
			a[i] = 0.5 * (self.eav + self.eap)
			self.eav = self.eap

		#les particules avancent de dti
		for i in range(self.ifirst, self.ilast):
			E = a[i]
			AA = (self.x[i] + E + r2 * self.v[i]) * exp(r2 * self.dti)
			BB = 2. * (self.x[i] + E - self.v[i]/r2) * exp(-self.dti/r2)
			self.x[i] = ((AA + BB) * tier) - E
			self.v[i] = (AA * r2 - BB/r2) * tier

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
		T=[]
		for i in range(self.ifirst, self.ilast):
			T.append([self.x[i],self.v[i],self.name[i]])

		T = self.tri_fusion(T)
		
		for i in range(self.ifirst, self.ilast):
			self.x[i],self.v[i],self.name[i]=T[i-1][0],T[i-1][1],T[i-1][2]


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
		self.x = [0] * self.m
		self.v = [0] * self.m
		self.mi = [0] * self.m
		self.ma = [0] * self.m
		self.name = [0] * self.m
		"""self.x = numpy.zeros((self.m), dtype='d')
		self.v = numpy.zeros((self.m), dtype='d')
		self.mi = numpy.zeros((self.m), dtype='d')
		self.ma = numpy.zeros((self.m), dtype='d')
		self.name = numpy.zeros((self.m), dtype='l')"""


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

	thread_count = 4
	for i in range(thread_count):
		univex = Univexp()
		univex.start()
	univex.join()
	#plt.show()
	# Affichage du temps d execution
	print("Temps d execution : %s secondes" % (time.time() - start_time))
main()
"""import cProfile
import re
cProfile.run('main()')"""

exit()
