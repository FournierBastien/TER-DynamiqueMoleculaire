OpenMP -> l'implémentation en OpenMP est plus lente et si valeur du tableau plus grande, l'utilisation de la mémoire/procs devient trop grande et fait planter l'ordinateur

MPI -> encore qu'une ébauche, les deux boucles while ont été remplacées par deux conditions IF permettant de paralléliser le tri en MPI, mais le programme ne s'arrête pas et écrit un fichier trop volumineux.