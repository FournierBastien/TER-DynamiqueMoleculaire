#include <iostream> 
#include <iomanip>
#include <fstream> 
#include <string> 
#include <cstdio> 
#include <cstdlib>
#include <ctime>
#include <cmath> 
#include <chrono>
#include <random>
#include<algorithm>
#include<limits.h>
#include <mpi.h>
#include <tuple>  

using namespace std; 

#define SMALL 32

inline void avance(int n, int m, int ifirst, int ilast, float dti, float * eav, float * eap, float * epolar, float * x, float * v) {

  float r2 = -sqrt(2.0);
  float tier = 1.0 / 3.0;
  float AA = 0.0 , BB = 0.0 , E = 0.0;

  float* a = new float[m];
  
  // calcul de l'accélération
  * eav =(* epolar) + 0.5 * n;

  for (int i = ifirst; i <= ilast; i++) {
    * eap = (* eav) - 1.0;
    a[i] = 0.5 * ( (* eav) + (* eap));
    * eav = (* eap);
  }

  for (int i = ifirst; i <= ilast; i++) {
    E = a[i];
    AA = ((x[i] + E + r2 * v[i]) * exp(r2 * dti));
    BB = 2.0 * (x[i] + E - v[i] / r2) * exp((-dti / r2));

    x[i] = ((AA + BB) * tier) - E;
    v[i] = (AA * r2 - BB / r2) * tier;

    if (x[i] > x[m-1]) {
      * epolar = (* epolar) + 1.0;
      x[i] = x[1] + x[i] - x[m-1];
    }

    if (x[i] < x[0]) {
      * epolar = (* epolar) - 1.0;
      x[i] = x[m-1] + x[i] - x[0];
    }

  }

  delete[] a;
  a = NULL;
}

inline void wbande(ofstream * flux_fichier, int ifirst, int ilast, int tecr, int m, int n, float * x, float * v, int * name) {
  * flux_fichier << "enregistrement a tecr=" << tecr;
  * flux_fichier << "#" << m << " " << n << " " << ifirst << " " << tecr;

  for (int i = ifirst; i <= ilast; i++) {
    * flux_fichier << fixed << setprecision(10) << x[i] << " " << v[i] << " " << name[i] << "\n";
  }

   * flux_fichier << "\n";
}

inline void init(ofstream * flux_fichier, int m, int n, int ifirst, int ilast, float pvit, float * x, float * v, float * mi, float * ma, int * name) {
  
  default_random_engine generator;
  uniform_real_distribution<float> distribution(0.0,1.0);

  x[0] = -0.5 * n;
  // x[0] = INT_MIN;
  x[m-1] = 0.5 * n;
  // x[m-1] = INT_MAX;
  v[0] = 0;
  name[0] = 0;

  v[m-1] = 0;
  name[m-1] = 0;

  for (int i = ifirst; i <= ilast; i++) {
    name[i] = i - 1;
    x[i] = x[0] + 0.5 + i - ifirst;
    v[i] = 2.0 * pvit * (0.5e0 - distribution(generator));
    //( * v)[i] = 2.0 * pvit * (0.5e0 - rand());
    mi[i] = 1.0;
    ma[i] = 1.0;
  }

  // calage du barycentre a zero avec une vitesse moyenne nulle

  float vmoy = 0.0;

  for (int i = ifirst; i <= ilast; i++) {
    vmoy += v[i];
  }
  vmoy = vmoy / n;

  for (int i = ifirst; i <= ilast; i++) {
    v[i] -= vmoy;
  }

  for (int i = ifirst; i <= ilast; i++) {
    vmoy += v[i];
  }

  // ( * flux_fichier) << " vmoyen = " << vmoy;
}


void fusion(float * x, float * v ,int * name, int size, float * tp_x,float * tp_v,int * tp_n) {
    int i1 = 0;
    int i2 = size/2;
    int temp = 0;

    while (i1 < size/2 && i2 < size) {
        if (x[i1] < x[i2]) {
            tp_x[temp] = x[i1];
            tp_v[temp] = v[i1];
            tp_n[temp] = name[i1];
            i1++;
        } else {
            tp_x[temp] = x[i2];
            tp_v[temp] = v[i2];
            tp_n[temp] = name[i2];
            i2++;
        }
        temp++;
    }
    while (i1 < size/2) {
        tp_x[temp] = x[i1];
        tp_v[temp] = v[i1];
        tp_n[temp] = name[i1];
        i1++;
        temp++;
    }
    while (i2 < size) {
        tp_x[temp] = x[i2];
        tp_v[temp] = v[i2];
        tp_n[temp] = name[i2];
        i2++;
        temp++;
    }
    // Copy sorted temp array into main array, a
    memcpy(x, tp_x, size*sizeof(float));
    memcpy(v, tp_v, size*sizeof(float));
    memcpy(name, tp_n, size*sizeof(int));
}

// permet de s'affranchir des problèmes de distributions avec des tableaux de taille impaire 
void tri_insertion(float * x,float * v,int * name, int size) {
  int i;
  for (i=0; i < size; i++) {
    int j, t3 = name[i];
    float t1 = x[i], t2 = v[i];
    for (j = i - 1; j >= 0; j--) {
      if (x[j] <= t1) break; 
      x[j + 1] = x[j];
      v[j + 1] = v[j];
      name[j + 1] = name[j];
    }
    x[j + 1] = t1;
    v[j + 1] = t2;
    name[j + 1] = t3;
  }
}

void tri_fusion(float * x, float * v ,int * name, int size, float * tp_x,float * tp_v,int * tp_n) {

  if (size < SMALL) {
      tri_insertion(x,v,name, size);
      return;
  }    
  tri_fusion(x,v,name, size/2, tp_x,tp_v,tp_n);
  tri_fusion(x + size/2,v + size/2, name + size/2, size - size/2, tp_x,tp_v,tp_n);

  fusion(x,v,name, size, tp_x,tp_v,tp_n);
}

int get_level_mpi(int my_rank) {
  int level = 0;
  while (pow(2, level) <= my_rank) level++;
  return level;    
}

// MPI tri fusion
void tri_fusion_mpi(float * x, float * v ,int * name, int size, float * tp_x,float * tp_v,int * tp_n, int level, int my_rank, int max_rank, int tag, MPI_Comm comm) {
  
  int helper_rank = my_rank + pow(2, level);
  if (helper_rank > max_rank) {
    tri_fusion(x,v,name, size, tp_x,tp_v,tp_n);
  } else {
    MPI_Request request; MPI_Status status;

    MPI_Isend(x + size/2, size - size/2, MPI_FLOAT, helper_rank, tag, comm, &request);
    MPI_Isend(v + size/2, size - size/2, MPI_FLOAT, helper_rank, tag+1, comm, &request);
    MPI_Isend(name + size/2, size - size/2, MPI_INT, helper_rank, tag+2, comm, &request);

    tri_fusion_mpi(x,v,name, size/2, tp_x,tp_v,tp_n, level+1, my_rank, max_rank, tag, comm);

    MPI_Recv(x + size/2, size - size/2, MPI_FLOAT, helper_rank, tag, comm, &status);
    MPI_Recv(v + size/2, size - size/2, MPI_FLOAT, helper_rank, tag+1, comm, &status);
    MPI_Recv(name + size/2, size - size/2, MPI_INT, helper_rank, tag+2, comm, &status);

    fusion(x,v,name, size, tp_x,tp_v,tp_n);
  } 
  return;
}

// executer par les procs non root
void tri_non_root_mpi(int my_rank, int max_rank, int tag, MPI_Comm comm) {
    int level = get_level_mpi(my_rank);

    MPI_Status status; int size;
    MPI_Probe(MPI_ANY_SOURCE, tag, comm, &status);
    MPI_Get_count(&status, MPI_INT, &size);
    cout << size << endl;

    // recuperer le rang de la source du message
    int parent_rank = status.MPI_SOURCE;

    float * x = new float[size];
    float * v = new float[size];
    int * name = new int[size];

    float *tp_x = new float[size];
    float *tp_v = new float[size];
    int *tp_n = new int[size];

    //
    MPI_Recv(x, size, MPI_FLOAT, parent_rank, tag, comm, &status);
    MPI_Recv(v, size, MPI_FLOAT, parent_rank, tag+1, comm, &status);
    MPI_Recv(name, size, MPI_INT, parent_rank, tag+2, comm, &status);
    tri_fusion_mpi(x,v,name, size, tp_x,tp_v,tp_n, level, my_rank, max_rank, tag, comm);  
    
    // On envoie les tableaux triés au proc parent
    MPI_Send(x, size, MPI_FLOAT, parent_rank, tag, comm);
    MPI_Send(v, size, MPI_FLOAT, parent_rank, tag+1, comm);
    MPI_Send(name, size, MPI_INT, parent_rank, tag+2, comm);
    return;
}


// executer par le proc root
void tri_root_mpi (float * x, float * v ,int * name, int size, float * tp_x,float * tp_v,int * tp_n, int max_rank, int tag, MPI_Comm comm) {
  int my_rank;
  MPI_Comm_rank(comm, &my_rank);
  if (my_rank != 0) {
      printf("Error: run_root_mpi called from process %d; must be called from process 0 only\n", my_rank);
      MPI_Abort(MPI_COMM_WORLD, 1);
  }
  tri_fusion_mpi(x,v,name, size, tp_x,tp_v,tp_n, 0, my_rank, max_rank, tag, comm); 
  return;
}

inline void algorithme(float * x, float * v,int * name, int size){
    int initialized, finalized;

    MPI_Initialized(&initialized);
    if (!initialized){
        MPI_Init(NULL, NULL);
    }
   
    int comm_size; MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int my_rank; MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int max_rank = comm_size - 1;
    int tag = 123, root_rank = 0;
    if (my_rank == 0) { 
	
        float *tp_x = new float[size];
        float *tp_v = new float[size];
        int *tp_n = new int[size];
        tri_root_mpi(x,v,name, size, tp_x,tp_v,tp_n, max_rank, tag, MPI_COMM_WORLD);
	  }
    else {
      tri_non_root_mpi(my_rank, max_rank, tag, MPI_COMM_WORLD);     
    }

    MPI_Finalized(&finalized);
    if (!finalized){
        MPI_Finalize();
    }
}

inline void run(string fichier, int n) {
  
  int initialized, finalized;

  MPI_Initialized(&initialized);
  if (!initialized){
      MPI_Init(NULL, NULL);
  }

  int comm_size; MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  int my_rank; MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int max_rank = comm_size - 1;
  int root = 0;
  int tag = 123, root_rank = 0;

  int continuer = 1;


  int m = n + 2;
  int ifirst = 1;
  int ilast = n;
  float dti = 0.001;
  float dtsor = 0.1;
  float tstop = 1.0;
  //float gravplas = 1.0;
  float pvit = 100.0;
  float epolar = 0.0;
  float eav = 0.0;
  float eap = 0.0;

  float tecr = 0.0;
  float tsor = 0.0;

  ofstream flux_fichier;

  // float * x;
  // float * v;
  // float * mi;
  // float * ma;
  // int * name;

  float * x = new float[m];
  float * v = new float[m];
  float * mi = new float[m];
  float * ma = new float[m];
  int * name = new int[m];

  if (my_rank == 0) { 
    flux_fichier.open(fichier);

    // flux_fichier << "simulation : " << fichier << " ";
    // flux_fichier << dti << " : pas de temps";
    // flux_fichier << " n = " << n << " pvit = " << pvit;
    // flux_fichier << " tstop = " << tstop << " dtsor = " << dtsor;

    // int continuer = 0;
    // float * x = new float[m];
    // float * v = new float[m];
    // float * mi = new float[m];
    // float * ma = new float[m];
    // int * name = new int[m];

    init( & flux_fichier, m, n, ifirst, ilast, pvit, x, v, mi, ma, name);
  }

  tecr=tecr+dti;
  tsor=tecr+dtsor;

  while(continuer == 1){
    cout << "boucle" << " " << continuer << endl;
    if (my_rank == 0) { 
      cout << "root" << endl;
      if(abs(tecr - tstop) > dtsor / 2.0){
        if(abs(tecr - tsor) > dti / 2.0){
          avance(n, m, ifirst, ilast, dti, & eav, & eap, & epolar, x, v);

          float *tp_x = new float[m];
          float *tp_v = new float[m];
          int *tp_n = new int[m];
          // distribuer le tableau
          run_root_mpi(x, v, name, m, tp_x, tp_v, tp_n, max_rank,tag, MPI_COMM_WORLD);
          // trier

          // recupérer 
          tecr += dti;
        }
        wbande( &flux_fichier, ifirst, ilast, tecr, m, n, x, v, name);
        tsor += dtsor;
      }else{
        continuer = 0;
      }


    }else{
      cout << "pas_root" << endl;
      run_helper_mpi(my_rank,max_rank,tag,MPI_COMM_WORLD);
    }

    cout << "après if" << endl;
    MPI_Bcast(&continuer,1,MPI_INT,root,MPI_COMM_WORLD);
    // if (my_rank == root) {
    //   int i;
    //   for (i = 0; i < max_rank; i++) {
    //     if (i != my_rank) {
    //       MPI_Send(&continuer, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    //     }
    //   }
    // } else {
    //   MPI_Recv(&continuer, 1, MPI_INT, root, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // }
  }

  // // init
  // init( & flux_fichier, m, n, ifirst, ilast, pvit, x, v, mi, ma, name);

  // while (abs(tecr - tstop) > dtsor / 2.0) {
  //   while (abs(tecr - tsor) > dti / 2.0) {
  //     algorithme( &flux_fichier, n, m, ifirst, ilast, & eav, & eap, & epolar, x, v, name, &tecr, &tstop, &dtsor, &dti, &tsor);
  //     tecr += dti;
  //   }
  //   wbande( & flux_fichier, ifirst, ilast, tecr, m, n, x, v, name);
  //   tsor = tecr+dtsor;

  // }
  
  if(my_rank == 0){
    flux_fichier << "\nFin du programme.\n";
    flux_fichier.close();

    delete[] x;
    delete[] v;
    delete[] mi;
    delete[] ma;
    delete[] name;
  }


  MPI_Finalized(&finalized);
  if (!finalized){
      MPI_Finalize();
  }

}

int main(void) {

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  string fichier = "resultat_test_tri_fusion_mpi.txt";
  int n = 10000;

  run(fichier, n);

  end = std::chrono::system_clock::now();
  int temps_ecoule = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
  cout << "Temps écoulé : " << temps_ecoule << " secondes"<< endl;

  return 0;
}