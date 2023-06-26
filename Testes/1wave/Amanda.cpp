#include <iostream>
#include <cmath>
#include <complex>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/version.hpp>
#include <omp.h>

using namespace std;
typedef array<double,2> state_type;
using namespace boost::numeric::odeint;

const double omega1 = 0.2;
const double omega2 = 1.2;
const double kx1 = 1.0;
const double kx2 = 1.0;
const double ky1 = 1.0;
const double ky2 = 1.0;
const double u = 1.0;
const double A1 = 1.0;
const double A2 = 0.5;
const double U = 0.0;
const double u1 = omega1/ky1;
struct mima
{
    double m_x,m_y;

    mima( double x= 0.0,double y=0.0 )
    : m_x( x),m_y( y) { }

    void operator()( const state_type &H , state_type &dHdt , double t ) const
    {
        dHdt[0] = ky1*A1*sin(kx1*H[0])*sin(ky1*H[1])+A2*ky2*sin(kx2*H[0])*sin(ky2*(H[1]-u*t)); 
        dHdt[1] = kx1*A1*cos(kx1*H[0])*cos(ky1*H[1])+A2*kx2*cos(kx2*H[0])*cos(ky2*(H[1]-u*t)); 
    }
};
void write_func( const state_type &x , const double t )
{
    if(fmod(t,10)==0){
	cout<<t<<"\t"<<x[0]<<"\t"<<x[1]<<endl;
	}
}

int main() {

    vector<state_type> x_vec;
    vector<double> times;
    state_type x;
    vector<double> psi2;
	double x_i,x_f,y_i,y_f;
	x_i = 0;
	x_f = M_PI;
	y_i = 0;
	y_f = 2*M_PI;
	double x_step = (x_f-x_i)/10;
	double y_step = (y_f-y_i)/10;
	adams_bashforth_moulton< 4 , state_type > stepper;
	   
	for (double x_p=x_i; x_p<=x_f; x_p+=x_step) {
		for (double y_p=y_i; y_p<=y_f; y_p+=y_step) {
			x[0] =x_p; 
   			x[1] =y_p; 
			integrate_adaptive( stepper , mima() , x , 0.0 ,100000.0 , 1e-3, write_func);
		
		}
	}    		



    
 
    
    
}
