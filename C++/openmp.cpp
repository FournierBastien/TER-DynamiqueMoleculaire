#include <iostream> 
#include <iomanip>
#include <fstream> 
#include <cstring>
#include <cstdio> 
#include <cstdlib>
#include <cmath> 
#include <chrono>
#include <random>

#include <omp.h>

using namespace std; 

#define MIN 512

// fusion de deux tableaux
inline void fusion(float * x, float * v ,int * name, int size, float * tp_x,float * tp_v, int * tp_n) {
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

// utilisation du tri insertion en fonction de la taille du tableau 
inline void tri_insertion(float * x,float * v,int * name, int size) {
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

// implémantation du tri fusion
inline void tri_fusion(float * x, float * v ,int * name, int size, float * tp_x,float * tp_v, int * tp_n) {

    if (size < MIN) {
       tri_insertion(x,v,name, size);
       return;
    }    
    tri_fusion(x,v,name, size/2, tp_x,tp_v,tp_n);
    tri_fusion(x + size/2,v + size/2,name + size/2 , size - size/2, tp_x,tp_v,tp_n);

    fusion(x,v,name, size, tp_x,tp_v,tp_n);
}

// Tri fusion OpenMP
inline void tri_fusion_omp(float * x, float * v ,int * name, int size, float * tp_x,float * tp_v, int * tp_n, int nb_procs) {
  if ( nb_procs == 1) {
    tri_fusion(x,v,name, size, tp_x,tp_v,tp_n);
  } else if (nb_procs > 1) {

    // 

    #pragma omp parallel sections
    {
      #pragma omp section
      { 
        tri_fusion_omp(x,v,name, size/2, tp_x,tp_v,tp_n, nb_procs/2);
      }

      #pragma omp section
      { 
        tri_fusion_omp(x + size/2,v + size/2,name + size/2, size - size/2, tp_x + size/2,tp_v + size/2,tp_n + size/2, nb_procs - nb_procs/2);
      }

    }
    fusion(x,v,name, size, tp_x,tp_v,tp_n); 
  } else {
      return;
  }
}

inline void avance(int n, int m, int ifirst, int ilast, float dti, float * eav, float * eap, float * epolar, float * x, float * v) {

  float r2 = -sqrt(2.0), tier = 1.0 / 3.0;
  float AA = 0.0 , BB = 0.0 , E = 0.0;

  float* a = new float[m];
  
  // calcul de l'accélération
  * eav = * epolar + 0.5 * n;

  for (int i = ifirst; i <= ilast; i++) {
    * eap = * eav - 1.0;
    a[i] = 0.5 * ( * eav + * eap);
    * eav = * eap;
  }

  for (int i = ifirst; i <= ilast; i++) {
    E = a[i];
    AA = ((x[i] + E + r2 * v[i]) * exp(r2 * dti));
    BB = 2.0 * (x[i] + E - v[i] / r2) * exp((-dti / r2));

    x[i] = ((AA + BB) * tier) - E;
    v[i] = (AA * r2 - BB / r2) * tier;

    if (x[i] > x[m-1]) {
      * epolar = * epolar + 1.0;
      x[i] = x[0] + x[i] - x[m-1];
    }

    if (x[i] < x[0]) {
      * epolar = * epolar - 1.0;
      x[i] = x[m-1] + x[i] - x[0];
    }

  }

  delete[] a;
  a = NULL;
}

inline void wbande(ofstream * flux_fichier, int ifirst, int ilast, int tecr, int m, int n, float * x, float * v, int * name) {
  // * flux_fichier << "enregistrement a tecr=" << tecr;
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
  x[m-1] = 0.5 * n;
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


inline void run(string fichier, int n) {
  
  int m = n + 2;
  int ifirst = 1;
  int ilast = n;
  float dti = 0.001;
  float dtsor = 1;
  float tstop = 15.0;
  //float gravplas = 1.0;
  float pvit = 100.0;
  float epolar = 0.0;
  float eav = 0.0;
  float eap = 0.0;

  float * x = new float[m];
  float * v = new float[m];
  float * mi = new float[m];
  float * ma = new float[m];
  int * name = new int[m];
  float tecr = 0.0;
  float tsor = 0.0;

  ofstream flux_fichier;
  flux_fichier.open(fichier);

  // flux_fichier << "simulation : " << fichier << " ";
  // flux_fichier << dti << " : pas de temps";
  // flux_fichier << " n = " << n << " pvit = " << pvit;
  // flux_fichier << " tstop = " << tstop << " dtsor = " << dtsor;

  init( & flux_fichier, m, n, ifirst, ilast, pvit, x, v, mi, ma, name);
  wbande( & flux_fichier, ifirst, ilast, tecr, m, n, x, v, name);


  omp_set_num_threads(2);

  tecr += dti;
  tsor = tecr+dtsor;

  while (abs(tecr - tstop) > dtsor / 2.0) {
    while (abs(tecr - tsor) > dti / 2.0) {
      avance(n, m, ifirst, ilast, dti, & eav, & eap, & epolar, x, v);

      float* tp_x = new float[m];
      float* tp_v = new float[m];
      int* tp_n = new int[m];

      // run_omp(x,v,name, m-1, tp_x,tp_v,tp_n,8);
      tri_fusion_omp(x,v,name, m-1, tp_x,tp_v,tp_n, 8);

      tecr += dti;
    }
    wbande( & flux_fichier, ifirst, ilast, tecr, m, n, x, v, name);
    tsor += dtsor;
  }

  flux_fichier << "\nFin du programme.\n";
  flux_fichier.close();

  delete[] x;
  delete[] v;
  delete[] mi;
  delete[] ma;
  delete[] name;
}


int main(void) {

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  string fichier = "resultat_test_tri_fusion_openmp.txt";
  int n = 10000;

//   on initialise le nombre de threads
//   omp_set_num_threads(atoi(argv[1]));

  run(fichier, n);
  
  end = std::chrono::system_clock::now();
  int temps_ecoule = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
  cout << "Temps écoulé : " << temps_ecoule << " secondes"<< endl;

  return 0;
}