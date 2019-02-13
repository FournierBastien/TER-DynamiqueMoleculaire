from threading import Thread
import random


class Univexp(Thread) :

    def __init__(self, n) :

        Thread.__init__(self)
        self.n = n
        self.fichier = "univexp01"

        self.file_1 = None
        self.file_2 = None

        # initialisation
        self.m = n+2
        self.ifirst = 2
        self.ilast = n+1
        self.dti = .001
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

        # fichier = open(self.fichier + ".txt", "a")
        # fichier.write("Simulation :")
        # fichier.close()

        print("simulation : " + self.fichier)
        print(self.dti + " : pas de temps")
        print(" n = " + self.n + " pvit = " + self.pvit)
        print(" tstop = " + self.tstop + " dtsor = " + self.dtsor)

        # initialisation des "particules-mur"

        self.x[1] = self.x[1] - 0.5*self.n
        self.x[self.m] = 0.5*self.n

        # initialisation des particules
        
        for i in range(self.ifirst,self.ilast) :
            self.name[i] = i-1
            self.x[i] = self.x[1] + 0.5+i-self.ifirst
            self.v[i] = 2.0 * self.pvit * (0.5 - random.random())
            self.mi[i] = 1.0
            self.ma[i] = 1.0

        # calage du barycentre à zero avec une vitesse moyenne nulle

        vmoy = sum(self.v[self.ifirst:self.ilast]) / self.n

        for i in range(self.ifirst,self.ilast) :
            self.v[i] = self.v[i] - vmoy
        vmoy = sum(self.v[self.ifirst:self.ilast])
        print(" vmoyen = " + vmoy)

    def run(self) :

        self.file_1 = open(self.fichier + "a.d", "w")
        self.file_2 = open(self.fichier + "d.d", "w")

        self.wbande()

        self.tecr = self.tecr + self.dti
        self.tsor = self.tsor + self.dtsor

        # while(abs(self.tecr-self.tstop))


    def wbande(self) :
        print(" enregistrement à tecr= " + self.tecr)
        self.file_1.write("# " + self.m + " " + self.n + " " + self.ifirst + " " + self.ilast + " " + self.tecr)

        for i in range(self.ifirst,self.ilast) :
            self.file_1.write(" " + self.x[i] + " " + self.v[i] + " " + self.name[i])
        
        self.file_1.write("/")



def main() :
    univ = Univexp(10000)


main()
exit()