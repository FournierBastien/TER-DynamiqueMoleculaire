from threading import Thread
import random
import math


class Univexp(Thread) :

    def __init__(self, n) :

        Thread.__init__(self)
        self.n = n
        self.fichier = "univexp01"

        self.file_1 = None
        self.file_2 = None

        # initialisation
        self.m = self.n+2
        self.ifirst = 2
        self.ilast = n+1
        self.dti = .001
        self.dtsor = 1.0
        self.tstop = 15.0
        self.gravplas = 1.0
        self.pvit = 100.0

        self.epolar = 0

        self.x = [0.0] * self.m
        self.v = [0.0] * self.m
        self.mi = [0.0] * self.m
        self.ma = [0.0] * self.m
        self.names = [0.0] * self.m

        self.tecr  = 0.0
        self.tsor = 0.0

        # fichier = open(self.fichier + ".txt", "a")
        # fichier.write("Simulation :")
        # fichier.close()

        print("simulation : " + self.fichier)
        print(str(self.dti) + " : pas de temps")
        print(" n = " + str(self.n) + " pvit = " + str(self.pvit))
        print(" tstop = " + str(self.tstop) + " dtsor = " + str(self.dtsor))

        # initialisation des "particules-mur"

        self.x[1] = self.x[1] - 0.5*self.n
        self.x[0] = 0.5*float(self.n)

        # initialisation des particules
        
        for i in range(self.ifirst,self.ilast):
            self.names[i] = i-1
            self.x[i] = self.x[1] + 0.5+i-self.ifirst
            self.v[i] = 2.0 * self.pvit * (0.5 - random.random())
            self.mi[i] = 1.0
            self.ma[i] = 1.0

        # calage du barycentre à zero avec une vitesse moyenne nulle

        vmoy = sum(self.v[self.ifirst:self.ilast]) / self.n

        for i in range(self.ifirst,self.ilast) :
            self.v[i] = self.v[i] - vmoy
        vmoy = sum(self.v[self.ifirst:self.ilast])
        print(" vmoyen = " + str(vmoy))

        self.run()

    def run(self) :

        self.file_1 = open(self.fichier + "a.d", "w")
        # self.file_2 = open(self.fichier + "d.d", "w")

        self.wbande()

        self.tecr = self.tecr + self.dti
        self.tsor = self.tsor + self.dtsor

        while(abs(self.tecr-self.tstop)>self.dtsor/2.0):
            print("abs(self.tecr-self.tstop)>self.dtsor/2.0: {}-{}={}>{}/2.0".format(self.tecr,self.tstop,self.tecr-self.tstop,self.dtsor))
            while(abs(self.tecr-self.tsor)>self.dti/2.0):
                self.avance()
                self.ordonne()
                self.tecr=self.tecr+self.dti
            self.wbande()
            self.tsor=self.tsor+self.dtsor
            
        self.file_1.close()	
        print("\nfin du programme.\n")


    def wbande(self) :
        # print(" enregistrement à tecr= " + str(self.tecr))
        # self.file_1.write("# " + str(self.m) + " " + str(self.n) + " " + str(self.ifirst) + " " + str(self.ilast) + " " + str(self.tecr))

        # for i in range(self.ifirst,self.ilast) :
        #     self.file_1.write(" " + str(self.x[i]) + " " + str(self.v[i]) + " " + str(self.names[i]))
        
        # self.file_1.write("\n")

        print("enregistrement a tecr={:7.3f}".format(self.tecr))
		
        self.file_1.write(("#{:5d} {:5d} {:5d} {:5d} {:7.3f}\n").format(self.m, self.n, self.ifirst, self.ilast, self.tecr))
        for i in range(self.ifirst, self.ilast):
            self.file_1.write(("{} {} {}\n".format(self.x[i],self.v[i],self.names[i])))
        self.file_1.write("\n")


    # blocage dans ordonne, indice out of range ?
    def ordonne(self) :

        # xp = [0] * (self.ilast - self.ifirst+1)
        

        #vp = [0] * (self.ilast - self.ifirst+1)

        j = 0

        for i in range(self.ifirst+1,self.ilast) :

            xp = self.x[i]
            vp = self.v[i]
            np = self.names[i]

            while(j != 1.0 and self.x[i] < self.x[j-1]) :
                j = j-1
            for k in range(j+1,i) :
                self.x[k] = self.x[k-1]
                self.v[k] = self.v[k-1]
                self.names[k] = self.names[k-1]

            self.x[j] = xp
            self.v[j] = vp
            self.names[j] = np

    def avance(self) :
        
        AA = 0.0
        BB = 0.0

        r2 = -math.sqrt(2.0)
        tier = 1.0/3.0

        a = [0] * self.m

        # calcul de l'accélération
        eav = self.epolar + 0.5 * self.n
        for i in range(self.ifirst,self.ilast) :
            eap = eav - 1
            a[i] = 0.5 * (eav+eap)
            eav = eap
        # les particules avances de dti

        for i in range(self.ifirst,self.ilast) :
            E = a[i]
            AA = (self.x[i] + E + r2 * self.v[i]) * math.exp(r2 * self.dti)
            BB = 2.0 * (self.x[i] + E - self.v[i] / r2) * math.exp(r2*self.dti)
            self.x[i] = ((AA + BB) * tier) - E
            self.v[i] = (AA * r2 - BB / r2) * tier

            if(self.x[i] > self.x[self.m-1]) :
                self.epolar = self.epolar + 1
                self.x[i] = self.x[1] + self.x[i] - self.x[self.m-1]
            if(self.x[i] < self.x[i]) :
                self.epolar = self.epolar - 1
                self.x[i] = self.x[self.m] + self.x[i] - self.x[1]


def main() :
    univ = Univexp(10000)


main()
exit()