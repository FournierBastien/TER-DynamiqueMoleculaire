# coding: utf-8


import sys
import os
from multiprocessing import Process
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from math import sqrt, exp
import numpy as np
import time
import random
import cProfile
import re

sys.setrecursionlimit(99999999)


class Univexp():

	fichier="univexp01"
	tecr=0.
	tsor=0.
	n=100

	def __init__(self):

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
		vmoy=sum(self.v[self.ifirst:self.ilast])
		print(("vmoyen = {:2.8f}").format(vmoy))

		self.tecr = self.tecr + self.dti
		self.tsor = self.tecr + self.dtsor


		self.prog_principal()


     #************************************************************************************************


	def prog_principal(self):

		self.f = open(self.fichier, "w")
		self.wbande()
		self.i=1
		while abs(self.tecr -self.tstop) > self.dtsor / 2.:
			while abs(self.tecr -self.tsor) > self.dti / 2.:
				self.avance()

				#Pour afficher le mouvement des particules
				#self.display(self.x, self.v)

				self.ordonne()
				self.tecr = self.tecr + self.dti
				
				self.i+=1
			self.wbande()
			self.tsor = self.tsor + self.dtsor
		self.f.close()

		print ("fin de programme.")

    #************************************************************************************************
	#Pour plot
	def display(self,x,v):
		plt.scatter(x, v)
		plt.savefig('images_plot/plot'+str(self.i)+'.png')
		plt.draw()
		plt.pause(0.01)


    #************************************************************************************************


	def avance(self):
		r2 = -sqrt(2.)
		tier = 1./3.
		#AA=BB=E=0.
		
		#calcul de l'acceleration
		a = np.zeros((self.m), dtype='d')
		self.eav = self.epolar + 0.5 * self.n
		for i in range(self.ifirst, self.ilast):
			#eap = eav - self.gravplas * self.ma[i]
			#a[i] = 0.5 * (eav + eap) / self.mi[i]
			self.eap = self.eav - 1
			a[i] = 0.5 * (self.eav + self.eap)
			self.eav = self.eap

		#les particules avancent de dti
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

                



	def tri_ins(self,t1,t2,t3):
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
		return [t1,t2,t3]

	

	def fusion(self,T1,T2):
		if T1==[] :
			return T2
		if T2==[]:
			return T1
		if T1[0][0]<T2[0][0]:
			return [T1[0]]+self.fusion(T1[1 :],T2)
		else:
			return [T2[0]]+self.fusion(T1,T2[1 :])
				
	def tri_fusion(self,T):
		if len(T)<2:
			return T
		n=len(T)

		return self.fusion(self.tri_fusion(T[:n//2]),self.tri_fusion(T[n//2:]))
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
#Animation
"""u=Univexp()
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,xlim=(-50, 50), ylim=(-50, 50))
particles, = ax.plot([], [])
def init():
	#initialize animation
	global u
	particles.set_data([], [], 'bo', ms=6)
	return particles

def animate(i):
	#perform animation step
	global u,ax, fig
	ms = int(fig.dpi * 2 * box.size * fig.get_figwidth()/ np.diff(ax.get_xbound())[0])
	# update pieces of the animation
	particles.set_data(u.x,u.y)
	particles.set_markersize(ms)
	return particles

ani = animation.FuncAnimation(fig, animate, frames=600,interval=10, blit=True, init_func=init)

plt.show()"""



def main() :

	# Debut du decompte du temps
	start_time = time.time()

	u=Univexp()
	#p = Process(target=u.display, args=(u.x, u.v))
	#plt.savefig('images_plot/plot'+str(self.i)+'.png')
	#p.start()
	
	# Affichage du temps d execution
	print("Temps d execution : %s secondes \n" % (time.time() - start_time))
	#plt.show()
	#p.join()
#main()

cProfile.run('main()')

#exit()
