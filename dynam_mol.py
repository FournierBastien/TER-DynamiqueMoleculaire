# -*- coding: utf-8 -*
from math import sqrt
from math import exp
import random


class Univexp :
    fichier ="uniexp01"
    tecr=0.0
    tsor=0.0
    n=10000

    def __init__(self) :

        self.init()
        self.file_1 = open(self.fichier, "w")
        self.file_2 = open(self.fichier, "w")

        self.wbande()

        self.tecr = self.tecr + self.dti
        self.tsor = self.tsor + self.dtsor

        while(abs(self.tecr-self.tstop) > self.dtsor/2.0):
            while(abs(self.tecr - self.tsor) > self.dti/2.0):
                self.avance()
                self.ordonne()
                self.tecr=self.tecr+self.dti
            self.wbande()
            self.tsor=self.tsor+self.dtsor
        self.file_1.close()
        self.file_2.close()
        print ("\nfin du programme.\n")
        
    def init(self):
        self.m = self.n+2
        self.ifirst = 1
        self.ilast = self.n
        self.dti = 0.001
        self.dtsor = 1.0
        self.tstop = 15.0
        self.gravplas = 1.0
        self.pvit = 100.0

        self.epolar = 0

        self.x = [0] * self.m
        self.v = [0] * self.m
        self.mi = [0] * self.m
        self.ma = [0] * self.m
        self.name = [0] * self.m

        self.tecr  = 0.0
        self.tsor = 0.0


        print("simulation : " + self.fichier)
        print(self.dti , " : pas de temps")
        print(" n = " , self.n, " pvit = " , self.pvit)
        print(" tstop = " , self.tstop , " dtsor = " , self.dtsor)

        # initialisation des "particules-mur"

        self.x[0] = self.x[0] - 0.5*float(self.n)
        self.x[self.m-1] = 0.5*float(self.n)

        # initialisation des particules
        
        for i in range(self.ifirst,self.ilast+1) :
            #self.name[i] = i-1
            self.x[i] = self.x[1] + 0.5+i-(self.ifirst+1)
            self.v[i] = 2.0 * self.pvit * (0.5 - random.random())
            self.mi[i] = 1.0
            self.ma[i] = 1.0

        # calage du barycentre à zero avec une vitesse moyenne nulle

        vmoy = sum(self.v[self.ifirst:self.ilast]) 
        vmoy= vmoy/ float(self.n)

        for i in range(self.ifirst,self.ilast+1) :
            self.v[i] = self.v[i] - vmoy
        vmoy = sum(self.v[self.ifirst:self.ilast])
        print(" vmoyen = ", vmoy)


    def avance(self):
         a = [0] * self.m
         AA=0
	 BB=0
	 E=0
         r2= - sqrt(2.0)
         tier = 1.0/3.0
         self.eav=self.epolar+0.5*self.n
         for i in range(self.ifirst, self.ilast+1):
            self.eap=self.eav-1
            a[i]=0.5*(self.eav+self.eap)
            self.eav=self.eap

         for i in range(self.ifirst, self.ilast+1):
            E= a[i]
            AA=(self.x[i]+E+r2*self.v[i])*exp(r2*self.dti)
            BB=2.0*(self.x[i]+E-self.v[i]/r2)*exp(-self.dti/r2)
            self.x[i]=((AA+BB)*tier)-E
            self.v[i]=(AA*r2-BB/r2)*tier
            if(self.x[i] > self.x[self.m-1]):
                self.epolar=self.epolar+1.0
                self.x[i]=self.x[1]+self.x[i]-self.x[self.m-1]
            if (self.x[i] < self.x[1]):
                #self.epolar = self.epolar - self.gravplas*self.ma[i]
                self.epolar = self.epolar - 1.0
                self.x[i] = self.x[self.m-1] + self.x[i] - self.x[1]




    def wbande(self) :
        print(" enregistrement à tecr= ", self.tecr)
        self.file_1.write("#  %s  %s  %s  %s  %s" %(self.m ,self.n,self.ifirst,self.ilast,self.tecr))

        for i in range(self.ifirst,self.ilast) :
            self.file_1.write(" %s  %s  %s" %(self.x[i],self.v[i],self.name[i]))
        
        self.file_1.write("/")

  
    def ordonne(self) :
        j=0
	xp=0.0
	vp=0.0
        for i in range(self.ifirst+1,self.ilast) :
            j=i
            xp = self.x[i]
            vp = self.v[i]
            np = self.name[i]

            while(j != 1 and self.x[i] < self.x[j-1]) :
                j = j-1
            for k in range(i,j+2,-1) :
                self.x[k] = self.x[k-1]
                self.v[k] = self.v[k-1]
                self.name[k] = self.name[k-1]

            self.x[j] = xp
            self.v[j] = vp
            self.name[j] = np





def main() :
    univ = Univexp()
  


main()