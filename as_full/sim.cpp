/*
 * =====================================================================================
 *
 *       Filename:  simulation.cpp
 *
 *    Description:  Contains all calculations for the simulation
 *
 *        Version:  1.0
 *        Created:  10/21/2011 09:25:10 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Pablo Marchant
 *        Company:  Pontificia Universidad Católica de Chile
 *
 * =====================================================================================
 */
#include <iostream>
#include <stdlib.h>
#include "math.h"
#include "sim.h"
#include "io.h"
#include "initial.h"
#include <gsl/gsl_sf.h>
namespace sim{
//Values given as cli arguments.
int rNum;
int thNum;
double dt;
int tNum;
int plotSteps;
double thtd;
//Values of the function beta on the 2-d grid.
double **B;
#ifndef TOROIDAL
//Values of the function alpha on the 2-d grid
double **A;
#endif
//Arrays that contain the precalculated quantities neccesary to solve evolution of A.
double **res_term_A;
#ifndef PUREOHM
double **hall_term_A;
#endif
//Arrays that contain the precalculated quantities neccesary to solve fluxes for B
#ifndef PUREOHM
double **hall_rflux;
double **hall_thflux;
#endif
double **res_rflux;
double **res_thflux;
//Sines precalculated at each point in the grid
double *sines;
//Minimun radius of the shell containing the magnetic field.
double rmin;
//Number of points in the radial direction both at the surface and the inner boundary for which the values
//will be solved by interpolation, and not direct calculation through the time evolution equation.
int rless;
#ifndef TOROIDAL
#ifndef SIMPLE
//Number of points used for the multipole fit outside the star
int l;
//Values of the coefficients that give the poloidal field outside the star.
double *a;
#endif
#endif
//Size of spatial steps.
double dr;
double dth;
//The value of pi
double const Pi=4*atan(1);
//value of radius and theta in a point of the grid, used multiple times
double r,th;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initial_conditions
 *  Description:  Solves the sizes of the timesteps, initializes arrays with correct sizes, and calls the functions from the namespace "initial" which gives the initial conditions for the simulation.
 * =====================================================================================
 */
	void
initial_conditions ( )
{
	//Initialize arrays with appropiate sizes
#ifndef TOROIDAL
	A=new double*[rNum];
#ifndef SIMPLE
	a=new double[l];
#endif
#endif
	B=new double*[rNum];
	for(int i=0;i<rNum;i++){
		B[i]=new double[thNum];
#ifndef TOROIDAL
		A[i]=new double[thNum];
#endif
	}

	//Get minimun radius from the value set in initial.cpp
	rmin=initial::rmin;

	//Solve size of steps in the radial and angular direction
	dr=(1.0-rmin)/(rNum-1);
	dth=Pi/(thNum-1);

	//Store the initial values of physical quantities
	for(int i=1;i<rNum-1;i++){
		r=rmin+i*dr;
		for(int j=1;j<thNum-1;j++){
			th=j*dth;
			B[i][j]=initial::B(r,th);
#ifndef TOROIDAL
			A[i][j]=initial::A(r,th);
#endif
		}
	}

	//Set boundary conditions for toroidal field
	for(int i=0;i<rNum;i++){
		B[i][0]=B[i][thNum-1]=0;
	}
	for(int j=0;j<thNum;j++){
		B[0][j]=B[rNum-1][j]=0;
	}
#ifndef TOROIDAL
	//Set boundary condition for poloidal fields
	for(int i=0;i<rNum;i++){
		A[i][0]=A[i][thNum-1]=0;
	}
	for(int j=0;j<thNum;j++){
		A[0][j]=A[rNum-1][j]=0;
	}
#endif

#ifndef TOROIDAL
#ifndef SIMPLE
	//Adjust alpha boundary conditions to fit smoothly to a poloidal field outside the star
	solve_A_boundary();
#endif
#endif
	
	return;
}		/* -----  end of function initial_conditions  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  solve_repeated_values
 *  Description:  Solves values that are constant along the entire simulation, in order to optimize the calculations.
 * =====================================================================================
 */
	void
solve_repeated_values ( )
{
	//Initialize arrays with appropiate sizes
#ifndef PUREOHM
	hall_term_A=new double*[rNum];
	hall_rflux=new double*[rNum];
	hall_thflux=new double*[rNum];
#endif
	res_term_A=new double*[rNum];
	res_rflux=new double*[rNum];
	res_thflux=new double*[rNum];
	sines=new double[thNum];
	for(int i=0;i<rNum;i++){
#ifndef PUREOHM
		hall_term_A[i]=new double[thNum];
		hall_rflux[i]=new double[thNum];
		hall_thflux[i]=new double[thNum];
#endif
		res_term_A[i]=new double[thNum];
		res_rflux[i]=new double[thNum];
		res_thflux[i]=new double[thNum];
	}
	//Solve common terms involved in the calculation of the toroidal field.
	for(int j=1;j<thNum-1;j++){
		th=j*dth;
		sines[j]=sin(th);
	}
	for(int i=0;i<rNum-1;i++){
		r=rmin+i*dr;
		for(int j=0;j<thNum-1;j++){
			th=j*dth;
#ifndef PUREOHM
			hall_term_A[i][j]=dt*sines[j]*initial::chi(r,th)/4/dr/dth;
			hall_rflux[i][j] = dt*initial::chi(r+dr/2,th)/8.0/dr/dth;
			hall_thflux[i][j]=-dt*initial::chi(r,th+dth/2)/8.0/dr/dth;
#endif
			res_term_A[i][j]=dt*thtd*initial::eta(r,th);
			res_rflux[i][j]  = dt*thtd*initial::eta(r+dr/2,th)/sin(th)/dr/dr;
			res_thflux[i][j] = dt*thtd*initial::eta(r,th+dth/2)/r/r/sin(th+dth/2.0)/dth/dth;
		}
	}
	return;
}		/* -----  end of function solve_repeated_values  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  simulate
 *  Description:  Performs the simulation. If it completes all required timesteps, returns 0, else, it returns 1.
 * =====================================================================================
 */
	int
simulate ( )
{
	//Create file where integrated quantities are logged
	io::create_integrals_file();

	//Begin simulation
	double t;
#ifndef TOROIDAL
	double Aaux[rNum][thNum];
	double gsA[rNum][thNum];
#endif
	double Baux[rNum][thNum];
	double dBr=0;
	double dBth=0;
	for(int k=0;k<=tNum;k++){
		//Log data if k is multiple of plotSteps
		if(k%plotSteps==0){
			t=k*dt;
			//Log integrated quantities
			io::log_integrals_file(t,solve_integrals());
			//Log complete profiles for A and B
			io::log_field(k);
			//No need to keep simulating if no output will be produced in next steps
			if(k+plotSteps>tNum){
				break;
			}
			std::cout << k << "/" << tNum << std::endl;
		}
#ifndef TOROIDAL
		//Update poloidal field function
		for(int i=1;i<rNum-1;i++){
			r=rmin+i*dr;
			for(int j=1;j<thNum-1;j++){
				th=j*dth;
				//Solve Grad-Shafranov operator at point
				gsA[i][j]=(A[i+1][j]+A[i-1][j]-2*A[i][j])/dr/dr+1/r/r*(A[i][j+1]+A[i][j-1]-2*A[i][j])/dth/dth-1/r/r*cos(th)/sin(th)*(A[i][j+1]-A[i][j-1])/2/dth;
				//Evolve poloidal field function at point
				Aaux[i][j]=A[i][j]+res_term_A[i][j]*gsA[i][j];
#ifndef PUREOHM
				Aaux[i][j]+= ((B[i][j+1]-B[i][j-1])*(A[i+1][j]-A[i-1][j])
							 -(B[i+1][j]-B[i-1][j])*(A[i][j+1]-A[i][j-1]))*hall_term_A[i][j];
#endif
			}
		}
#endif

		//Update toroidal field function
		for(int i=0;i<rNum-1;i++){
			for(int j=0;j<thNum-1;j++){
				//Solve radial fluxes on B/sin(th)
				if(j!=0){
					dBr=res_rflux[i][j]*(B[i+1][j]-B[i][j]);
#ifndef PUREOHM
					dBr+=hall_rflux[i][j]
						*(B[i][j]+B[i+1][j])
						*(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1]);
#ifndef TOROIDAL
					dBr+=hall_rflux[i][j]
						*(gsA[i][j]+gsA[i+1][j])
						*(A[i][j+1]+A[i+1][j+1]-A[i][j-1]-A[i+1][j-1]);
#endif
#endif
				}
				//Solve theta fluxes on B/sin(th)
				if(i!=0){
					dBth=res_thflux[i][j]*(B[i][j+1]-B[i][j]);
#ifndef PUREOHM
					dBth+=hall_thflux[i][j]
						*(B[i][j]+B[i][j+1])
						*(B[i+1][j]+B[i+1][j+1]-B[i-1][j]-B[i-1][j+1]);
#ifndef TOROIDAL
					dBth+=hall_thflux[i][j]
						*(gsA[i][j]+gsA[i][j+1])
						*(A[i+1][j]+A[i+1][j+1]-A[i-1][j]-A[i-1][j+1]);
#endif
#endif
				}
				//Add flux at point, substract it from following points
				Baux[i][j]+=(dBr+dBth)*sines[j];
				Baux[i+1][j]=B[i+1][j]-dBr*sines[j];
				Baux[i][j+1]-=dBth*sines[j+1];
			}
		}
		//Pass values from auxiliary array Baux and Aaux
		for(int i=1;i<rNum-1;i++){
			for(int j=1;j<thNum-1;j++){
				//Check for blowups, exit program if that happens
				if(isinf(Baux[i][j])
#ifndef TOROIDAL
						||isinf(Aaux[i][j])
#endif
						){
					io::report_blowup(k,i,j);
					return 1;
				}
				B[i][j]=Baux[i][j];
#ifndef TOROIDAL
				A[i][j]=Aaux[i][j];
#endif
			}
		}
		//Solve B values close to the inner shell through 3-point interpolation
		double a1,a2;
		double f1,f2;
		for(int j=1;j<thNum-1;j++){
			f1=B[rless+1][j];
			f2=B[rless+2][j];
			a1=((f1-f2)*pow(rless,2)+(4*f1-2*f2)*rless-f2+4*f1)/dr/(pow(rless,2)+3*rless+2);
			a2=-((f1-f2)*rless-f2+2*f1)/pow(dr,2)/(pow(rless,2)+3*rless+2);
			for(int n=1;n<=rless;n++){
				r=n*dr;
				B[n][j]=a1*r+a2*pow(r,2);
			}
			f1=B[rNum-1-rless-1][j];
			f2=B[rNum-1-rless-2][j];
			a1=-((f1-f2)*pow(rless,2)+(4*f1-2*f2)*rless-f2+4*f1)/dr/(pow(rless,2)+3*rless+2);
			a2=-((f1-f2)*rless-f2+2*f1)/pow(dr,2)/(pow(rless,2)+3*rless+2);
			for(int n=1;n<=rless;n++){
				r=-n*dr;
				B[rNum-1-n][j]=a1*r+a2*pow(r,2);
			}
		}
#ifndef TOROIDAL
#ifndef SIMPLE
		solve_A_boundary();
#endif
#endif
	}

	//Close file where integrated quantities are logged
	io::close_integrals_file();
	return 0;
}		/* -----  end of function simulate  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  solve_integrals
 *  Description:  Solves integrated quantities in the star:
 *  - Toroidal flux
 *  - Toroidal field energy
 *  - Poloidal field energy
 *  These are returned in an array in that order.
 * =====================================================================================
 */
	double*
solve_integrals ( )
{
	//Initialize array with adequate size, and set values to zero
#ifdef TOROIDAL
	double *integrals=new double[2];
	integrals[0]=integrals[1]=0;
#else
#ifndef SIMPLE
	//When the poloidal field is considered, the individual energy of each multipole outside the star,
	//and the sum of all the multipoles is also given.
	double *integrals=new double[4+l];
	for(int n=0;n<4+l;n++){
		integrals[n]=0;
	}
#else
	double *integrals=new double[3];
	integrals[0]=integrals[1]=integrals[2]=0;
#endif
#endif

	//Solve quantities integrated over the volume of the star
	for(int i=1;i<rNum-1;i++){
#ifndef TOROIDAL
		r=i*dr+rmin;
#endif
		for(int j=1;j<thNum-1;j++){
			//Toroidal flux
			integrals[0]+=B[i][j]/sines[j];
			//Toroidal energy
			integrals[1]+=pow(B[i][j],2)/sines[j];
#ifndef TOROIDAL
			//Poloidal energy
			integrals[2]+=pow((A[i+1][j]-A[i-1][j])/2/dr,2)/sines[j]+pow((A[i][j+1]-A[i][j-1])/2/dth,2)/r/r/sines[j];
#endif
		}
	}
	//Multiply by the area element on the r-th grid, and also by (2*Pi) and 1/(8*pi) the energies, where
	//the first factor is from integration over phi, and the second one is the factor present in the calculation
	//of the magnetic energy.
	integrals[0]=integrals[0]*dr*dth;
	integrals[1]=integrals[1]*dr*dth/4;
#ifndef TOROIDAL
	integrals[2]=integrals[2]*dr*dth/4;
#ifndef SIMPLE
	//Solve external poloidal energy using the coefficients of the expansion outside the star
	for(int n=0;n<l;n++){
		integrals[4+n]=pow(a[n],2)*(n+2)/(8*Pi);
		integrals[3]+=integrals[4+n];
	}
#endif
#endif
	return integrals;
}		/* -----  end of function solve_integrals  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  solve_A_boundary
 *  Description:  Solves the boundary values of A such that the field connects smoothly to a potential field outside the star
 * =====================================================================================
 */
#ifndef TOROIDAL
#ifndef SIMPLE
	void
solve_A_boundary ( )
{
	//Solve the required legendre polynomial on all required points
	double P[l+1][thNum];
	for(int i=0;i<=l;i++){
		for(int j=1;j<thNum;j++){
			P[i][j]=gsl_sf_legendre_Pl(i,cos(j*dth));
		}
	}
	//Solve integrated coefficients. Integrations are performed using simpsons rule, which uses a quadratic function
	//to aproximate the integrand (simpson's rule).
	double f0=0,f1=0,f2=0;
	for(int n=0;n<l;n++){
		a[n]=0;
		for(int j=1;j<thNum-2;j++){
			th=j*dth;
			f1=A[rNum-2][j]*(cos(j*dth)*P[n+1][j]-P[n][j])/sin(j*dth);
			f0=f2=0;
			if(j!=thNum-2){
				f2=A[rNum-2][j+1]*(cos(th+dth)*P[n+1][j+1]-P[n][j+1])/sin(th+dth);
			}
			if(j!=1){
				f0=A[rNum-2][j-1]*(cos(th-dth)*P[n+1][j-1]-P[n][j-1])/sin(th-dth);
			}
			a[n]+=(f0+4*f1+f2)*dth/3;
		}
		a[n]=a[n]*2*Pi*pow(1-dr,n+1)/(n+2)*sqrt((2*n+3)/4.0/Pi)*(n+1)/2;
	}
	for(int j=1;j<thNum-1;j++){
		A[rNum-1][j]=0;
	}
	//Permorm sumation
	for(int n=0;n<l;n++){
		for(int j=1;j<thNum-1;j++){
			A[rNum-1][j]+=a[n]*sqrt((2*n+3)/4.0/Pi)*(cos(j*dth)*P[n+1][j]-P[n][j]);
		}
	}
	return;
}		/* -----  end of function solve_A_boundary  ----- */
#endif
#endif

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  release_memory
 *  Description:  Destroy all used pointers to assure there are no memory leaks
 * =====================================================================================
 */
	void
release_memory ( )
{
	for(int i=0;i<rNum;i++){
		delete[] B[i];
		B[i]=NULL;
#ifndef PUREOHM
		delete[] hall_rflux[i];
		hall_rflux[i]=NULL;
		delete[] hall_thflux[i];
		hall_thflux[i]=NULL;
#endif
		delete[] res_rflux[i];
		res_rflux[i]=NULL;
		delete[] res_thflux[i];
		res_thflux[i]=NULL;
#ifndef TOROIDAL
		delete[] A[i];
		A[i]=NULL;
		delete[] res_term_A[i];
		res_term_A[i]=NULL;
#ifndef PUREOHM
		delete[] hall_term_A[i];
		hall_term_A[i]=NULL;
#endif
#endif
	}
	delete B;
	B=NULL;
#ifndef PUREOHM
	delete hall_rflux;
	hall_rflux=NULL;
	delete hall_thflux;
	hall_thflux=NULL;
#endif
	delete res_rflux;
	res_rflux=NULL;
	delete res_thflux;
	res_thflux=NULL;
#ifndef TOROIDAL
	delete A;
	A=NULL;
	delete res_term_A;
	res_term_A=NULL;
#ifndef PUREOHM
	delete hall_term_A;
	hall_term_A=NULL;
#endif
#ifndef SIMPLE
	delete[] a;
	a=NULL;
#endif
#endif
}		/* -----  end of function release_memory  ----- */

}		/* -----  end of namespace sim  ----- */
