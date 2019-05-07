# -*- coding: utf-8 -*
from math import sqrt, exp
from decimal import *
import random
import os
import sys
import time

start_time = time.time()

getcontext().prec = 9

class Univexp :
    fichier ="univexp01.txt"
    n=3
    tecr=Decimal(0.)
    tsor=Decimal(0.)

    def __init__(self) :

        self.init()
        self.file_1 = open(self.fichier, "w")
       
        self.wbande()

        self.tecr = self.tecr + self.dti
        self.tsor = self.tsor + self.dtsor

        while(abs(self.tecr-self.tstop) > self.dtsor/Decimal(2.)):
            while(abs(self.tecr - self.tsor) > self.dti/Decimal(2.)):
                self.avance()
                self.ordonne()
                self.tecr=self.tecr+self.dti
            self.wbande()
            self.tsor=self.tsor+self.dtsor
        self.file_1.close()
       
        print ("\nfin du programme.\n")
        
    def init(self):
        self.m = self.n+2
        self.ifirst = 1
        self.ilast = self.n+1
        self.dti = Decimal(0.001)
        self.dtsor = Decimal(.1)
        self.tstop = Decimal(1.)
        self.gravplas = Decimal(1.)
        self.pvit = Decimal(100.)

       

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

        self.x[0] =  Decimal(- 0.5*self.n)
        self.v[0] = 0
        self.mi[0] = 0
        self.ma[0] = 0
        self.name[0] = 0
        self.x[self.m-1] = Decimal(0.5)*self.n
        self.v[self.m-1] = 0
        self.mi[self.m-1] = 0
        self.ma[self.m-1] = 0
        self.name[self.m-1] = 0

        # initialisation des particules
#        random.seed(100)
#        for i in range(self.ifirst,self.ilast) :
#            self.name[i] = i-1
#            self.x[i] = self.x[0] + 0.5+i-(self.ifirst)
#            self.v[i] = 2. * self.pvit * (0.5 - random.random())
#            self.mi[i] = 1.
#            self.ma[i] = 1.

        self.name[1] = 1
        self.name[2] = 2
        self.name[3] = 3

        self.x[1] = Decimal(-1.00000000)
        self.x[2] = Decimal(0.00000000)
        self.x[3] = Decimal(1.00000000)

        self.v[1] = Decimal(59.1418419)
        self.v[2] = Decimal(32.8358345)
        self.v[3] = Decimal(-91.9776764)



        self.epolar=Decimal(0.)
        # calage du barycentre à zero avec une vitesse moyenne nulle

        vmoy = Decimal(sum(self.v[self.ifirst:self.ilast])) 
        vmoy= Decimal(vmoy/ self.n)

        self.v[self.ifirst:self.ilast]=[val - vmoy for val in self.v[self.ifirst:self.ilast]]
        vmoy = sum(self.v[self.ifirst:self.ilast])
        print("vmoyen={:7.3f}".format(vmoy))


    def avance(self):
         a = [0] * self.m
         r2= Decimal(2).sqrt()
         tier = Decimal(1/3)
         self.eav=self.epolar+Decimal(0.5*self.n)
         for i in range(self.ifirst, self.ilast):
            self.eap=self.eav-1
            a[i]=Decimal(0.5)*(self.eav+self.eap)
            self.eav=self.eap

         for i in range(self.ifirst, self.ilast):
            E= a[i]
            AA=(Decimal(self.x[i])+Decimal(E)+Decimal(r2)*Decimal(self.v[i]))*Decimal(r2*self.dti).exp()
            BB=2*(Decimal(self.x[i])+Decimal(E)-Decimal(self.v[i])/Decimal(r2))*Decimal(-self.dti/r2).exp()
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
        xp=Decimal(0.)
        vp=Decimal(0.)
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

print("éxécuté en %s seconds " % (time.time() - start_time))
