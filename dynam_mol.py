# -*- coding: utf-8 -*
from math import sqrt, exp
import random
import os
import sys


class Univexp :
    fichier ="univexp01"
    n=10000
    tecr=0.
    tsor=0.

    def __init__(self) :

        self.init()
        self.file_1 = open(self.fichier, "w")
       
        self.wbande()

        self.tecr = self.tecr + self.dti
        self.tsor = self.tsor + self.dtsor

        while(abs(self.tecr-self.tstop) > self.dtsor/2.):
            while(abs(self.tecr - self.tsor) > self.dti/2.):
                self.avance()
                self.ordonne()
                self.tecr=self.tecr+self.dti
            self.wbande()
            self.tsor=self.tsor+self.dtsor
        self.file_1.close()
       
        print ("\nfin du programme.\n")
        
    def init(self):
        self.m = self.n+2
        self.ifirst = 2
        self.ilast = self.n+1
        self.dti = 0.001
        self.dtsor = 1.
        self.tstop = 15.
        self.gravplas = 1.
        self.pvit = 100.

       

        self.x = [0] * self.m
        self.v = [0] * self.m
        self.mi = [0] * self.m
        self.ma = [0] * self.m
        self.name = [0] * self.m

       


        print("simulation : " + self.fichier)
        print(self.dti , " : pas de temps")
        print(" n = " , self.n, " pvit = " , self.pvit)
        print(" tstop = " , self.tstop , " dtsor = " , self.dtsor)

        # initialisation des "particules-mur"

        self.x[0] =  - 0.5*self.n
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
        random.seed(100)
        for i in range(self.ifirst,self.ilast) :
            self.name[i] = i-1
            self.x[i] = self.x[0] + 0.5+i-(self.ifirst)
            self.v[i] = 2. * self.pvit * (0.5 - random.random())
            self.mi[i] = 1.
            self.ma[i] = 1.
            
	self.epolar=0.
        # calage du barycentre à zero avec une vitesse moyenne nulle

        vmoy = sum(self.v[self.ifirst:self.ilast]) 
        vmoy= vmoy/ self.n

        self.v[self.ifirst:self.ilast]=[val - vmoy for val in self.v[self.ifirst:self.ilast]]
        vmoy = sum(self.v[self.ifirst:self.ilast])
	print("vmoyen={:7.3f}".format(vmoy))


    def avance(self):
         a = [0] * self.m
         r2= - sqrt(2.)
         tier = 1/3
         self.eav=self.epolar+0.5*self.n
         for i in range(self.ifirst, self.ilast):
            self.eap=self.eav-1
            a[i]=0.5*(self.eav+self.eap)
            self.eav=self.eap

         for i in range(self.ifirst, self.ilast):
            E= a[i]
            AA=(self.x[i]+E+r2*self.v[i])*exp(r2*self.dti)
            BB=2*(self.x[i]+E-self.v[i]/r2)*exp(-self.dti/r2)
            self.x[i]=((AA+BB)*tier)-E
            self.v[i]=(AA*r2-BB/r2)*tier
            if(self.x[i] > self.x[self.m-1]):
                self.epolar=self.epolar+1
                self.x[i]=self.x[0]+self.x[i]-self.x[self.m-1]
            if (self.x[i] < self.x[0]):
                #self.epolar = self.epolar - self.gravplas*self.ma[i]
                self.epolar = self.epolar - 1
                self.x[i] = self.x[self.m-1] + self.x[i] - self.x[0]




    def wbande(self) :
        print(" enregistrement a tecr= ", self.tecr)
        self.file_1.write("# {:<5d} {:<5d}  {:<5d} {:<5d} {:<7.3f} \n".format(self.m, self.n, self.ifirst, self.ilast, self.tecr))

        for i in range(self.ifirst,self.ilast) :
            self.file_1.write(" {:5.13f}    \t {:5.13f}     \t {}\n".format(self.x[i],self.v[i],self.name[i]))
        
        self.file_1.write("\n")

  
    def ordonne(self) :
        j=i=k=np=0
	xp=0.
	vp=0.
        for i in range(self.ifirst+1,self.ilast) :
            j=i
            xp = self.x[i]
            vp = self.v[i]
            np = self.name[i]

            while(j != 0 and self.x[i] < self.x[j-1]) :
                j = j-1
            for k in range(i,j,-1) :
                self.x[k] = self.x[k-1]
                self.v[k] = self.v[k-1]
                self.name[k] = self.name[k-1]

            self.x[j] = xp
            self.v[j] = vp
            self.name[j] = np





def main() :
    univ = Univexp()
  


main()
