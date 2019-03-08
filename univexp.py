#!/usr/bin/env python
# coding: utf-8


import sys
import os
from math import sqrt, exp
import numpy
from random import random
# Specification du dossier de travail
os.chdir("/media/mcd/MCD/M1Informatique/S2/TER/DynamiqueMoleculaire")

class Univexp():

    fichier="univexp01"
    tecr=0.
    tsor=0.
    n=10000

    def __init__(self):

        self.init()
        self.f = open(self.fichier, "w")

        self.wbande()

        self.tecr = self.tecr + self.dti
        self.tsor = self.tecr + self.dtsor

        while abs(self.tecr -self.tstop) > self.dtsor / 2.:
            
            #print("abs(self.tecr-self.tstop)>self.dtsor/2.0: {}-{}={}>{}/2.0".format(self.tecr,self.tstop,self.tecr-self.tstop,self.dtsor))
            while abs(self.tecr -self.tsor) > self.dti / 2.:
                self.avance()
                self.ordonne()
                self.tecr = self.tecr + self.dti
            self.wbande()
            self.tsor = self.tsor + self.dtsor
        self.f.close()
        print ("fin de programme.")
        

    #************************************************************************************************

    def avance(self):
        r2 = -sqrt(2.)
        tier = 1/3
        
        AA=BB=E=0.
        
        #calcul de l'acceleration
        a = numpy.zeros((self.m), dtype='d')
        self.eav = self.epolar + 0.5 * self.n
        for i in range(self.ifirst, self.ilast +1):
            #eap = eav - self.gravplas * self.ma[i]
            #a[i] = 0.5 * (eav + eap) / self.mi[i]
            self.eap = self.eav - 1
            a[i] = 0.5 * (self.eav + self.eap)
            self.eav = self.eap

        #les particules avancent de dti
        for i in range(self.ifirst, self.ilast +1):
            E = a[i]
            AA = (self.x[i] + E + r2 * self.v[i]) * exp(r2 * self.dti)
            BB = 2 * (self.x[i] + E - self.v[i]/r2) * exp(-self.dti/r2)
            self.x[i] = ((AA + BB) * tier) - E
            self.v[i] = (AA * r2 - BB/r2) * tier
            if self.x[i] > self.x[self.m-1]:
                self.epolar = self.epolar + 1
                self.x[i] = self.x[1] + self.x[i] - self.x[self.m-1]

            if self.x[i] < self.x[1] :
                #self.epolar self.epolar - self.gravplas * self.ma[i]
                self.epolar = self.epolar - 1
                self.x[i] = self.x[self.m-1] + self.x[i] - self.x[1]

                


    #************************************************************************************************

    def ordonne(self):
        
        i=j=k=np=0

        xp=vp=mip=map=0.

        for i in range(self.ifirst+1, self.ilast+1):
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
        #print(max(self.x) == self.x[self.ilast])
        print(self.x == sorted(self.x))
        print(self.x)
        print(sorted(self.x))
    
    
    #************************************************************************************************

    def wbande(self):
        print("enregistrement a tecr = {:7.3f}".format(self.tecr))

        self.f.write("# m = {:<5d} \t n = {:<5d} \t ifirst = {:<5d} \t ilast = {:<5d} \t tecr = {:<7.3f} \n".format(
                self.m, self.n, self.ifirst, self.ilast, self.tecr))
        self.f.write("\n")
            #"# m="+str(self.m)+ " \t n= "+str( self.n)+ " \t ifirst= "+str(self.ifirst)+ " \t ilast= "+str(self.ilast)+ " \t tecr= "+str(self.tecr)+"\n")
        
        for i in range(self.ifirst, self.ilast+1):
            self.f.write("x[{}]= {:<25}  v[{}]= {:<25}  name[{}]= {}\n".format(i,self.x[i],i,self.v[i],i,self.name[i]))

        self.f.write("\n\n")
                # self.f.write(("x["+str(i)+"]= "+str(self.x[i])+ " \t v["+str(i)+"]= "+str(self.v[i])+ " \t name["+str(i)+"]= "+str(self.name[i])+"\n")
            
        """with open('self01', 'a') as fichier510:

            fichier510.write("# m={:5d} \t n={:5d} \t ifirst={:5d} \t ilast={:5d} \t tecr={:7.3f} \n".format(
                self.m, self.n, self.ifirst, self.ilast, self.tecr))
            #"# m="+str(self.m)+ " \t n= "+str( self.n)+ " \t ifirst= "+str(self.ifirst)+ " \t ilast= "+str(self.ilast)+ " \t tecr= "+str(self.tecr)+"\n")

            for i in range(self.ifirst, self.ilast+1):
                fichier510.write("x["+str(i)+"]= "+str(self.x[i])+ " \t v["+str(i)+"]= "+str(self.v[i])+ " \t name["+str(i)+"]= "+str(self.name[i])+"\n")
            """

    #************************************************************************************************

    def init(self):
        self.m = self.n + 2
        self.ifirst = 1
        self.ilast=self.n+1
        self.dti = 0.001
        self.dtsor = 1.
        self.tstop = 15.
        self.gravplas=1.
        self.pvit = 100.

        #allocation pour les tableaux
        self.x = numpy.zeros((self.m), dtype='d')
        self.v = numpy.zeros((self.m), dtype='d')
        self.mi = numpy.zeros((self.m), dtype='d')
        self.ma = numpy.zeros((self.m), dtype='d')
        self.name = numpy.zeros((self.m), dtype='l')


        print("simulation :",self.fichier)

        print("pas de temps {:19.15f}".format(self.dti))
        print("n={:8d} pvit={:7.1f}".format(self.n, self.pvit))
        print("tstop={:7.1f} dtsor ={:7.1f}".format(self.tstop, self.dtsor))

        # initialisation des "particules-mur"
        self.x[1] = -0.5*self.n
        self.v[1] = 0.
        self.mi[1] = 0.
        self.ma[1] = 0.
        self.name[1] = 0
        self.x[self.m-1] = 0.5*self.n
        self.v[self.m-1] = 0.
        self.mi[self.m-1] = 0.
        self.ma[self.m-1] = 0.
        self.name[self.m-1] = 0.

        # initialisation des particules
        for i in range(self.ifirst -1, self.ilast+1):# pourquoi ?
            self.name[i] = i#i-1
            self.x[i] = self.x[1] + 0.5+1 - self.ifirst
            self.v[i] = 2. * self.pvit * (0.5 - random())# rand() de fortran =7.62939453E-06 sans graine
            self.mi[i] = 1.
            self.ma[i] = 1.

        self.epolar = 0.

        # calage du barycentre à zero avec une vitesse moyenne nulle
        vmoy = sum(self.v)
        vmoy = vmoy/self.n
 
        self.v=[i - vmoy for i in self.v]

        print("vmoyen={:7.3f}".format(vmoy))
        print("x init=",self.x)
        print("v init",self.v)
        print("n init",self.name)

def main() :
    Univexp()

main()

exit()
