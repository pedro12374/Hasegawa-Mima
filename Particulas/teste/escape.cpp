#include <stdio.h>  // biblioteca padrão
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <gsl/gsl_rng.h> // biblioteca p numeros aleatorios
#include <iomanip>

int phi1[4];// array that can hold 100 numbers for 1st column 
int phi2[4];// array that can hold 100 numbers for 2nd column 
int phi3[4];// array that can hold 100 numbers for 3rd column  
int t_phi[4];


double dydt2(double t,double x,double y,double *w, double *An, double *kx, double *ky){

  return An[0]*kx[0]*cos(kx[0]*x)*cos(ky[0]*y) + An[1]*kx[1]*cos(kx[1]*x)*cos(ky[1]*(y-(w[1]/ky[1] - w[0]/ky[0])*t))+ An[2]*kx[2]*cos(kx[2]*x)*cos(ky[2]*(y-(w[2]/ky[2] - w[0]/ky[0])*t));
}

double dxdt2(double t,double x,double y,double *w, double *An, double *kx, double *ky){

  return An[0]*ky[0]*sin(kx[0]*x)*sin(ky[0]*y) + An[1]*ky[1]*sin(kx[1]*x)*sin(ky[1]*(y-(w[1]/ky[1] - w[0]/ky[0])*t)) + An[2]*ky[2]*sin(kx[2]*x)*sin(ky[2]*(y-(w[2]/ky[2] - w[0]/ky[0])*t));
}


using namespace std;

int main(int argc, char const *argv[]) {


  double An[3], kx[3], ky[3], w[3]; // frequencais das demais ondas
  



   w[0] = 1.31e-3;
   w[1] = 1.31e-3;
   w[2] = 1.31e-3;
  
  kx[0] = 17.5;
  kx[1] = -7.5;
  kx[2] = -10.0;
  
  ky[0] = 5.0;
  ky[1] = -2.5;
  ky[2] = -2.5;
  
  
  //int c_s = 8;

  double t = 0;
  double step = 1e-3;      // passo temporal já normalizado


  double k1,k2,k3,k4;
  double l1,l2,l3,l4;


  int i = 0;

  // parametros pra cada particula
  double x = atof(argv[1]); // x' inicial (já normalizado)
  double y = atof(argv[2]); // y' inicial (já normalizadtmax = atof(argv[4]);  // numero de interaçõe
  double tmax = atof(argv[3]);
  //string out_folder = string(argv[4]); // saida do arquivo
  double x0 = x;
  double y0 = y;
  

  ifstream infile;   

  infile.open("phi_-0.211.dat");// file containing numbers in 3 columns 
     if(infile.fail()) // checks to see if file opended 
    { 
      cout << "error" << endl; 
      return 1; // no point continuing if the file didn't open...
    } 

double tfile;
  //  Loop de integracao
  while (t <= tmax) {
    if (x > 1.05) {
      cout << x0 << "\t" << y0 << "\t" << t << "\t" << fabs(fmod(y,2*M_PI)) << "\n";
      return 0;
    }
    /// Depois da quali, arrumar a normalização pelo fator B
    
    infile>> tfile;
    infile>> An[0];
    infile>> An[1];
    infile>> An[2];
    double k1 = dxdt2(t,x,y,w,An,kx,ky);
    double l1 = dydt2(t,x,y,w,An,kx,ky);

    double k2 = dxdt2(t+step/2,x + k1*step/2, y + l1*step/2,w,An,kx,ky);
    double l2 = dydt2(t+step/2,x + k1*step/2, y + l1*step/2,w,An,kx,ky);

    double k3 = dxdt2(t+step/2,x + k2*step/2, y + l2*step/2,w,An,kx,ky);
    double l3 = dydt2(t+step/2,x + k2*step/2, y + l2*step/2,w,An,kx,ky);

    double k4 = dxdt2(t+step,x + k3*step, y + l3*step,w,An,kx,ky);
    double l4 = dydt2(t+step,x + k3*step, y + l3*step,w,An,kx,ky);

    x +=  (k1 +  2*k2 + 2*k3 + k4)*step/6;
    y +=  (l1 +  2*l2 + 2*l3 + l4)*step/6;
    t += step;
    //cout<<t<<"\t"<<x<<"\t"<<fabs(fmod(y,2*M_PI))<<endl;
  }
  
  //cout << x0 << "\t" << y0 << "\t" << t << "\t" << y << "\n";

  return 0;
}
