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

using namespace std; 

inline void fusion(float * x, float * v, int * name, int start, int mid, int end) { 

  int start2 = mid + 1; 
  if (x[mid] <= x[start2]) return; 

  while (start <= mid && start2 <= end) { 

    if (x[start] <= x[start2]) start++; 
    else { 
      float x2 = x[start2]; 
      float v2 = v[start2];
      float n2 = name[start2];
      int index = start2; 

      while (index != start) { 
          x[index] = x[index - 1];
          v[index] = v[index - 1] ;
          name[index] = name[index - 1];
          index--; 
      } 
      x[start] = x2; 
      v[index] = v2 ;
      name[index] = n2;

      start++; 
      mid++; 
      start2++; 
    } 
  } 
} 

// version récursive du tri fusion, attention au nombre d'appel sur la pile ?
// ~ 8 secondes 
inline void triFusion(float * x, float * v, int * name, int start, int end) 
{ 

  // si le tableau est supérieur à un élément
  if (start < end) { 

    // calcule de l'indice du milieu
    int m = start + (end - start) / 2; 

    // partie gauche
    triFusion(x,v,name,start, m); 

    // partie droite
    triFusion(x,v,name, m + 1, end);

    // on fusionne les deux parties
    fusion(x,v,name, start, m, end); 
  } 
} 

// trie à bulle ~ 14 secondes
inline void ordonne(int ifirst, int ilast, float * x, float * v, int * name) {
  int np = 0;
  float xp = 0.0, vp = 0.0;
  // float mip, map;

  // tri à bulle
  // pour chaque valeurs du tableau, noté i
  // Tant qu'il existe une valeur plus grande que i, on déplace i de -1
  // Ensuite, on déplace i jusqu'à l'indice j
  
  for (int i = ifirst+1; i <= ilast; i++) {
    int j = i;
    xp = x[i];
    vp = v[i];

    np = name[i];

    while (j != 1.0 && x[i] < x[j - 1]) {
      j = j - 1.0;
    }

    // on décale de +1 toutes les valeurs comprises entre j et i
    for (int k = i; k > j + 1; k--) {
      x[k] = x[k - 1];
      v[k] = v[k - 1];
      name[k] = name[k - 1];
    }

    // on insert la valeur trié ensuite
    x[j] = xp;
    v[j] = vp;
    name[j] = np;
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

  ( * flux_fichier) << " vmoyen = " << vmoy;
}

inline void run(string fichier, int n) {
  int m = n + 2;
  int ifirst = 1;
  int ilast = n + 1;
  float dti = 0.001;
  float dtsor = 1.0;
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

  flux_fichier << "simulation : " << fichier << " ";
  flux_fichier << dti << " : pas de temps";
  flux_fichier << " n = " << n << " pvit = " << pvit;
  flux_fichier << " tstop = " << tstop << " dtsor = " << dtsor;

  // init
  init( & flux_fichier, m, n, ifirst, ilast, pvit, x, v, mi, ma, name);
  wbande( & flux_fichier, ifirst, ilast, tecr, m, n, x, v, name);

  tecr += dti;
  tsor = tecr+dtsor;

  while (abs(tecr - tstop) > dtsor / 2.0) {
    while (abs(tecr - tsor) > dti / 2.0) {
      avance(n, m, ifirst, ilast, dti, & eav, & eap, & epolar, x, v);
      //triFusion(x,v,name,0,m-1);
      ordonne(ifirst,ilast,x,v,name);
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

inline void generate(){

  int n = 10000;

  ofstream flux_fichier;
  flux_fichier.open("random_values.txt");

  default_random_engine generator;
  uniform_real_distribution<double> distribution(0.0,1.0);

  for( int i = 0; i < n; i++){
    flux_fichier << distribution(generator) << "\n";
  }

  flux_fichier.close();
}

int main(void) {

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  string fichier = "resultat_test_tri_bulle.txt";
  int n = 10000;

  run(fichier, n);

  end = std::chrono::system_clock::now();
  int temps_ecoule = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
  cout << "Temps écoulé : " << temps_ecoule << " secondes"<< endl;
  
  return 0;
}