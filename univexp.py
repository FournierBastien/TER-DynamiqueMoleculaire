#!/usr/bin/env python
# coding: utf-8


import sys
import os
from math import sqrt, exp
import numpy
# Specification du dossier de travail
os.chdir("/media/mcd/MCD/M1Informatique/S2/TER/DynamiqueMoleculaire")

class Univexp():
    def __init__(self):

        self.fichier = open("univexp01", "w")
        self.fichier.write(" dvcsdcd")
        self.fichier.close()
        #self.fichier = open("univexp01", "a")
        
        self.n = 10000
        self.m = self.n + 2
        self.ifirst = 2
        self.ilast=self.n +  1

        self.tecr = 0.
        self.tsor = 0.
        self.tstop = 15.
        self.dti = 0.001
        self.dtsor = 1.

        self.pvit = 100.
        self.epolar = self.eav = self.eap = 0.
        self.gravplas = 1.
        #allocation pour les tableaux
        self.x = numpy.zeros((self.m), dtype='d')
        self.v = numpy.zeros((self.m), dtype='d')
        self.mi = numpy.zeros((self.m), dtype='d')
        self.ma = numpy.zeros((self.m), dtype='d')
        self.name = numpy.zeros((self.m), dtype='l')
        print(type(self.name))

        #self.fichier500 = argv[1]
        #self.fichier510 = argv[2]




    #************************************************************************************************

    def avance(self):
        r2 = -sqrt(2.)
        tier = 1/3
        
        #AA=BB=E=0.
        

        #calcul de l'acceleration
        a = numpy.zeros((self.m), dtype='d')
        eav = self.epolar + 0.5 * self.n
        for i in range(self.ifirst, self.ilast +1):
            #eap = eav - self.gravplas * self.ma[i]
            #a[i] = 0.5 * (eav + eap) / self.mi[i]
            eap = eav - 1
            a[i] = 0.5 * (eav + eap)
            eav = eap

        #les particules avancent de dti
        for i in range(self.ifirst, self.ilast +1):
            E = a[i]
            AA = (self.x[i] + E + r2 * self.v[i]) * exp(r2 * self.dti)
            BB = 2 * (self.x[i] + E - self.v[i]/r2) * exp(-self.dti/r2)
            self.x[i] = ((AA + BB) * tier) - E
            self.v[i] = (AA * r2 - BB/r2) * tier
            if self.x[i] > self.x[self.m-1]:
                #self.epolar self.epolar + self.gravplas * self.ma[i]
                self.epolar = self.epolar + 1
                self.x[i] = self.x[1] + self.x[i] - self.x[self.m-1]
                #print("")
            #print("fin premier if")
            if self.x[i] < self.x[1] :
                #self.epolar self.epolar - self.gravplas * self.ma[i]
                self.epolar = self.epolar - 1
                self.x[i] = self.x[self.m-1] + self.x[i] - self.x[1]

                


    #************************************************************************************************

    def ordonne(self):
        """ 
        i=0
        j=0
        k=0
        np=0

        xp=0.0
        vp=0.0
        mip=0.0
        map=0.0 
        """

        for i in range(self.ifirst+1, self.ilast+1):
            j = i 
            xp=self.x[i]
            vp=self.v[i]
            #mip=self.mi[i]
            #map=self.ma[i]
            np=self.name[i]

            while j!=1 and self.x[i]<self.x[j-1]:
                j=j-1
        
            for k in range(i,j+2,-1):
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
        i=0
        print("enregistrement a tecr =", self.tecr)

        with open('univexp01', 'a') as fichier510:

            fichier510.write("# m= "+str(self.m)+ " \t n= "+str( self.n)+ " \t ifirst= "+str(self.ifirst)+ " \t ilast= "+str(self.ilast)+ " \t tecr= "+str(self.tecr)+"\n")

            for i in range(self.ifirst, self.ilast+1):
                fichier510.write("x["+str(i)+"]= "+str(self.x[i])+ " \t v["+str(i)+"]= "+str(self.v[i])+ " \t name["+str(i)+"]= "+str(self.name[i])+"\n")

            print("\n")

    #************************************************************************************************

    def initialisation(self):
        i=0
        vmoy=0.0


        

        # affichage sur l'ecran

        print("simulation :")

        print("pas de temps", self.dti)
        print("n=", self.n," pvit=", self.pvit)
        print("tstop=", self.tstop, " dtsor = ",self.dtsor)

        # initialisation des "particules-mur"
        self.x[1] = -0.5*self.n
        self.v[1] = 0.0
        self.mi[1] = 0.0
        self.ma[1] = 0.0
        self.name[1] = 0
        self.x[self.m-1] = 0.5*self.n
        self.v[self.m-1] = 0.0
        self.mi[self.m-1] = 0.0
        self.ma[self.m-1] = 0.0
        self.name[self.m-1] = 0.0

        # initialisation des particules
        for i in range(self.ifirst, self.ilast+1):
            self.name[i] = i-1
            self.x[i] = self.x[1] + 0.5+1 - self.ifirst
            self.v[i] = 2.0 * self.pvit * (0.5 - 7.62939453E-06)# rand() de fortran =7.62939453E-06
            self.mi[i] = 1.0
            self.ma[i] = 1.0

        self.epolar = 0.0

        # calage du barycentre à zero avec une vitesse moyenne nulle
        vmoy = sum(self.v)
        vmoy = vmoy/self.n
 
        self.v=[i - vmoy for i in self.v]

        print("vmoyen=", vmoy)