#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/version.hpp>
#include <omp.h>

using namespace std;
int N = 5;
typedef array< complex<double>, 5> state_type;
using namespace boost::numeric::odeint;


struct mima
{
    double *m_gamma,*m_omega,*m_lambda,*m_kx, *m_ky;
	
    mima( double omega[3],double lambda[3],double gamma[3] ,double kx[3],double ky[3])
    : m_omega( omega ),m_lambda( lambda ),m_gamma( gamma ), m_kx(kx), m_ky(ky) { }

    void operator()( const state_type &H , state_type &dHdt , double t ) const
    {
        const complex< double > I( 0.0 , 1.0 );
       	dHdt[0] = -I*m_omega[0]*H[0]+m_lambda[0]*conj(H[1])*conj(H[2])+m_gamma[0]*H[0];
		dHdt[1] = -I*m_omega[1]*H[1]+m_lambda[1]*conj(H[0])*conj(H[2])+m_gamma[1]*H[1];
	   	dHdt[2] = -I*m_omega[2]*H[2]+m_lambda[2]*conj(H[0])*conj(H[1])+m_gamma[2]*H[2];
		dHdt[3] = (m_ky[0]*H[0].real()*sin(m_kx[0]*H[3]+m_ky[0]*H[4])+m_ky[0]*H[0].imag()*cos(m_kx[0]*H[3]+m_ky[0]*H[4])+ m_ky[1]*H[1].real()*sin(m_kx[1]*H[3]+m_ky[1]*H[4])+m_ky[1]*H[1].imag()*cos(m_kx[1]*H[3]+m_ky[1]*H[4]) +  m_ky[2]*H[2].real()*sin(m_kx[2]*H[3]+m_ky[2]*H[4])+m_ky[2]*H[2].imag()*cos(m_kx[2]*H[3]+m_ky[2]*H[4])      ) ;
		dHdt[4]= -5.0*H[3]+ (-1.0)*(m_kx[0]*H[0].real()*sin(m_kx[0]*H[3]+m_ky[0]*H[4])+m_kx[0]*H[0].imag()*cos(m_kx[0]*H[3]+m_ky[0]*H[4])+ m_kx[1]*H[1].real()*sin(m_kx[1]*H[3]+m_ky[1]*H[4])+m_kx[1]*H[1].imag()*cos(m_kx[1]*H[3]+m_ky[1]*H[4]) +  m_kx[2]*H[2].real()*sin(m_kx[2]*H[3]+m_ky[2]*H[4])+m_kx[2]*H[2].imag()*cos(m_kx[2]*H[3]+m_ky[2]*H[4])      ) ;
			
    }
};

template<class system>
struct observer
{
    state_type m_dxdt;
    system m_odefun;
	ofstream &m_file;
  
	observer(system odefun,ofstream &file) : m_odefun(odefun),m_file(file) { }
		
    void operator()( state_type &x, double t)
    {
				if(fmod(t,1)==0){
	
        m_file << t<<","; 
		for (int i =0; i<N; i++) {
		m_file <<x[i].real()<<","<<x[i].imag()<<",";
		m_odefun(x, m_dxdt, t);
		m_file<<m_dxdt[i].real()<<","<<m_dxdt[i].imag()<<",";
		}
		m_file <<endl;
		} 
    }
};

string doubleToString(double valor, int casas_decimais) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(casas_decimais) << valor;
    return ss.str();
}
int main(int argc,char* argv[]) {

    
    double x_i,x_f,y_i,y_f;
	x_i = 0;
	x_f = M_PI;
	y_i = 0;
	y_f = 2*M_PI;
	double gamma[3]={-0.211,0.01,-0.211};
	double lambda[3]={1,2,3};
	double omega[3] = {0,0,0};	
	
	double kx1 = 17.5;//atof(argv[1]);
	double kx2 = -7.5;//atof(argv[2]);
	double kx3 = -10;//atof(argv[3]);

	double ky1 = 5;//atof(argv[4]);
	double ky2 = -2.5;//atof(argv[5]);
	double ky3 = -2.5;//atof(argv[6]);

	double kx[3] = {kx1,kx2,kx3};
	double ky[3] = {ky1,ky2,ky3};
	double x_step = (x_f-x_i)/10;
	double y_step = (y_f-y_i)/10;
    

	omega[0]=1.31e-3;//(ky[0]-kx[0])*0.17/(1+kx[0]*kx[0]+ky[0]*ky[0]);
	omega[1]=1.31e-3;//(ky[1]-kx[1])*0.17/(1+kx[1]*kx[1]+ky[1]*ky[1]);
	omega[2]=1.31e-3;//(ky[2]-kx[2])*0.17/(1+kx[2]*kx[2]+ky[2]*ky[2]);

	lambda[0]=0.04;//(omega[0]/(2*0.17))*((kx[2]*kx[2]+ky[2]*ky[2]-kx[1]*kx[1]+ky[1]*ky[1])*(kx[1]*ky[2]-kx[2]*ky[1])/(ky[0]-kx[0]));
	lambda[1]=-0.5;//(omega[1]/(2*0.17))*((kx[0]*kx[0]+ky[0]*ky[0]-kx[2]*kx[2]+ky[2]*ky[2])*(kx[2]*ky[0]-kx[0]*ky[2])/(ky[1]-kx[1]));
	lambda[2]=0.4;//(omega[2]/(2*0.17))*((kx[1]*kx[1]+ky[1]*ky[1]-kx[0]*kx[0]+ky[0]*ky[0])*(kx[0]*ky[1]-kx[1]*ky[0])/(ky[2]-kx[2]));
	
	double x_p = M_PI/2;
	double y_p = M_PI+0.2;
	adams_bashforth_moulton< 8 , state_type > stepper;    
	string filename="kx_";
	for (double ks : kx) {
		filename+=doubleToString(ks, 1)+"_";
	}
	filename+="ky_";
	for (double ks : ky) {
		filename+=doubleToString(ks, 1)+"_";
	}
	filename+=".csv";
	ofstream myfile(filename);
	myfile<<"time"<<","<<"phi1_r"<<","<<"phi1_i"<<","<<"dphi1_r"<<","<<"dphi1_i"<<","<<"phi2_r"<<","<<"phi2_i"<<","<<"dphi2_r"<<","<<"dphi2_i"<<","<<"phi3_r"<<","<<"phi3_i"<<","<<"dphi3_r"<<","<<"dphi3_i"<<","<<"x_r"<<","<<"x_i"<<","<<"dx_r"<<","<<"dx_i"<<","<<"y_r"<<","<<"y_i"<<","<<"dy_r"<<","<<"dy_i"<<","<<endl;

	for (x_p=x_i; x_p<=x_f; x_p+=x_step) {
		for (y_p=y_i; y_p<=y_f; y_p+=y_step) {
	   		state_type x;   
			x[0] = complex<double>(0.1,0.0);
        	x[1] = complex<double>(0.1,0.0);
			x[2] = complex<double>(0.1,0.0);
			x[3] = complex<double>(x_p,0.0);
			x[4] = complex<double>(y_p,0.0);
			mima system(omega,lambda,gamma,kx,ky);
			observer< mima> obs(system,myfile);
			integrate_adaptive( stepper , system , x , 0.0 ,1000.0 , 1e-3, boost::ref(obs));
		}
	}
   	myfile.close();

 
    
    
}
