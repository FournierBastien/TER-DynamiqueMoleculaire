#include <iostream> 
#include <fstream> 
#include <string> 
#include <cstdio> 
#include <cstdlib>
#include <ctime>
#include <cmath> 
#include <vector>
#include<algorithm>

using namespace std;

void ordonne(int ifirst, int ilast, vector < float > * x, vector < float > * v, vector < int > * name) {
  int np;
  float xp, vp;
  //float mip, map;

  for (int i = ifirst+1; i <= ilast; i++) {
    int j = i;
    xp = ( * x)[i];
    vp = ( * v)[i];

    np = ( * name)[i];

    while (j != 1.0 && ( * x)[i] < ( * x)[j - 1]) {
      j = j - 1.0;
    }
    //cout << "T" << endl;
    for (int k = i; k > j + 1; k--) {
      ( * x)[k] = ( * x)[k - 1];
      ( * v)[k] = ( * v)[k - 1];
      ( * name)[k] = ( * name)[k - 1];
    }

    ( * x)[j] = xp;
    ( * v)[j] = vp;
    ( * name)[j] = np;
    
  }

}

void avance(int n, int m, int ifirst, int ilast, float dti, float * eav, float * eap, float * epolar, vector < float > * x, vector < float > * v) {

  float r2 = -sqrt(2.0), tier = 1.0 / 3.0;
  float AA, BB, E;

  vector < float > a(m);

  // calcul de l'accélération
  * eav = * epolar + 0.5 * n;

  for (int i = ifirst; i <= ilast; i++) {
    * eap = * eav - 1.0;
    a[i] = 0.5 * ( * eav + * eap);
    * eav = * eap;
  }

  for (int i = ifirst; i <= ilast; i++) {
    E = a[i];
    AA = (( * x)[i] + E + r2 * ( * v)[i]) * exp(r2 * dti);
    BB = 2.0 * (( * x)[i] + E - ( * v)[i] / r2) * exp((-dti / r2));

    ( * x)[i] = ((AA + BB) * tier) - E;
    ( * v)[i] = (AA * r2 - BB / r2) * tier;

    if (( * x)[i] > ( * x)[m-1]) {
      * epolar = * epolar + 1.0;
      ( * x)[i] = ( * x)[0] + ( * x)[i] - ( * x)[m-1];
    }

    if (( * x)[i] < ( * x)[0]) {
      * epolar = * epolar - 1.0;
      ( * x)[i] = ( * x)[m-1] + ( * x)[i] - ( * x)[0];
    }

  }

}

void wbande(ofstream * flux_fichier, int ifirst, int ilast, int tecr, int m, int n, vector < float > * x, vector < float > * v, vector < int > * name) {
  * flux_fichier << "enregistrement a tecr=" << tecr;
  * flux_fichier << "#" << m << " " << n << " " << ifirst << " " << tecr;

  for (int i = ifirst; i <= ilast; i++) {
    * flux_fichier << ( * x)[i] << " " << ( * v)[i] << " " << ( * name)[i] << "\n";
  }

   * flux_fichier << "\n" << endl;
}

void init(ofstream * flux_fichier, int m, int n, int ifirst, int ilast, float pvit, vector < float > * x, vector < float > * v, vector < float > * mi, vector < float > * ma, vector < int > * name) {
  
  /*( * x)[0] = 0.5 * n;
  ( * x)[1] += -0.5 * n;*/

  ( * x)[0] = -0.5 * n;
  ( * x)[m-1] = 0.5 * n;

  for (int i = ifirst; i <= ilast; i++) {
    ( * name)[i] = i - 1;
    ( * x)[i] = ( * x)[0] + 0.5 + i - ifirst;
    ( * v)[i] = 2.0 * pvit * (0.5e0 - rand());
    ( * mi)[i] = 1.0;
    ( * ma)[i] = 1.0;
  }

  // calage du barycentre a zero avec une vitesse moyenne nulle

  float vmoy;

  for (int i = ifirst; i <= ilast; i++) {
    vmoy += ( * v)[i];
  }
  vmoy = vmoy / n;

  for (int i = ifirst; i <= ilast; i++) {
    ( * v)[i] -= vmoy;
  }

  for (int i = ifirst; i <= ilast; i++) {
    vmoy += ( * v)[i];
  }
  //for_each((*v)[ifirst],(*v)[ilast], [&](int n){vmoy += n;});

  ( * flux_fichier) << " vmoyen = " << vmoy;

}

void run(string fichier, int n) {
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
  vector < float > * x = new vector < float > (m);
  vector < float > * v = new vector < float > (m);
  vector < float > * mi = new vector < float > (m);
  vector < float > * ma = new vector < float > (m);
  vector < int > * name = new vector < int > (m);
  float tecr = 0.0;
  float tsor = 0.0;

  ofstream flux_fichier;
  flux_fichier.open(fichier);

  flux_fichier << "simulation : " << fichier;
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
      ordonne(ifirst, ilast, x, v, name);
      tecr += dti;
      //cout << tecr << " " << abs(tecr - tstop) << " " << dti / 2.0 << endl;
    }
    //cout << "Temps" << endl;
    wbande( & flux_fichier, ifirst, ilast, tecr, m, n, x, v, name);
    tsor += dtsor;
    

  }

  flux_fichier << "\nFin du programme.\n";
  flux_fichier.close();
}

int main() {

  srand(1089);
  clock_t start = clock();
  string fichier = "resultat_test_cpp.txt";
  int n = 10000;

  run(fichier, n);


  double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  //t = clock() - t;
  cout << "Temps écoulé : " << duration << " secondes"<< endl;

  return 0;
}

// Temps écoulé : 241.908 secondes soit 4 minutes