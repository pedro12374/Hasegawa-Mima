#include <iostream>
#include <cmath>
#include <complex>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/version.hpp>

using namespace std;
typedef array< complex<double>, 3> state_type;
using namespace boost::numeric::odeint;

const double omega1 = 1.5;
const double omega2 = 1.31e-3;
const double omega3 = 1.31e-3;

const double lambda1 = 0.041;
const double lambda2 = -0.5;
const double lambda3 = 0.4;

const double kx1 = 2.0;
const double ky1 = 2.0;

struct mima
{
    

    mima(  )
      { }

    void operator()( const state_type &H , state_type &dHdt , double t ) const
    {
        const complex< double > I( 0.0 , 1.0 );
       	dHdt[0] = -I*omega1*H[0];

	   	dHdt[1] = (0.5)*ky1*(H[0].real()*sin(kx1*H[1]+ky1*H[2])+ H[0].imag()*cos(kx1*H[1]+ky1*H[2])    )  ; 
        dHdt[2] = (-0.5)*kx1*(H[0].real()*sin(kx1*H[1]+ky1*H[2])+ H[0].imag()*cos(kx1*H[1]+ky1*H[2])    ); 

	           //dxdt[1] = -I*omega2*x[1]+lambda2*conj(x[2])*conj(x[0])+m_gamma2*x[1];
        //dxdt[2] = -I*omega3*x[2]+lambda3*conj(x[0])*conj(x[1])+m_gamma3*x[2];
    }
};


void write_lorenz( const state_type &x , const double t )
{
    cout << t << '\t' << x[0].real() << '\t' << x[0].imag() << '\t' << x[1].real() << '\t' << x[1].imag() << '\t' << x[2].real()  << '\t' << x[2].imag()  << endl;
}
void write_func( const state_type &x , const double t )
{
   // if(fmod( t, 2*M_PI)>=-0.001 && fmod(t,2*M_PI      )<=0.001){
	cout<<t<<"\t"<<x[1].real()<<"\t"<<x[2].real()<<endl;
	//}
}
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

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
	double x_step = (x_f-x_i)/5;
	double y_step = (y_f-y_i)/5;
double x_p = M_PI/2;
double y_p = M_PI+0.5;
        vector<double> temp_psi;
        adams_bashforth_moulton< 4 , state_type > stepper;
                //for (double x_p=x_i; x_p<=x_f; x_p+=x_step) {
		//for (double y_p=y_i; y_p<=y_f; y_p+=y_step) {
		x[0] = complex<double>(0.1,0.0);
        x[1] = complex<double>(x_p,0.0);
        x[2] = complex<double>(y_p,0.0);

			 
			integrate_adaptive( stepper , mima() , x , 0.0 ,100.0 , 1e-3, write_func);
		
		//}
	//}
   	    
 
    
    
}
