#!/usr/bin/env python
#-*- coding: utf_8 -*-

import univexp
from univexp import *

#cd /media/mcd/MCD/M1Informatique/S2/TER/DynamiqueMoleculaire/

def main() :
    univexp = Univexp()

    univexp.initialisation()

    print("Fin methode intialisation")
    univexp.wbande()
    print("Fin methode wbande")
    univexp.tecr = univexp.tecr + univexp.dti
    univexp.tsor = univexp.tecr + univexp.dtsor
    i=0
    while abs(univexp.tecr -univexp.tstop) > univexp.dtsor / 2:
        while abs(univexp.tecr -univexp.tsor) > univexp.dti / 2:
            print("While",i)
            univexp.avance()
            print("While",i)
            univexp.ordonne()
            print("While",i)
            univexp.tecr = univexp.tecr + univexp.dti
            print("While")
        univexp.wbande()
        univexp.tsor = univexp.tsor + univexp.dtsor
        i=i+1
        print("While",i)
    print ("fin de programme.")

main()

exit()