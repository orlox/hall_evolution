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
#include <time.h>
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
double **B, **dBr, **dBth;
#ifndef TOROIDAL
//Values of the function alpha on the 2-d grid
double **A, **Aaux, **gsA;
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
//Sines and cotans precalculated at each point in the grid
double *sines, *cotans;
//Minimun radius of the shell containing the magnetic field.
double rmin;
#ifndef TOROIDAL
//Number of points used for the multipole fit outside the star
int l;
//Values of the coefficients that give the poloidal field outside the star.
double *a;
//Repeated l dependent values in the resolution of the multipole fit.
double *boundary_factors1, *boundary_factors2;
//Repeated l and th dependent values that are solved in term of combinations Legendre polynomials. This is solved in the grid, and midpoint, so two arrays are required.
double **legendre_comb;
double **legendre_comb_mid;
//Auxiliary variables used in function sim::solve_A_boundary to integrate the multipole coefficients
double f0,f1,f2;
#endif
//Size of spatial steps.
double dr,dth;
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
	A=new double*[rNum+1];
	Aaux=new double*[rNum+1];
	gsA=new double*[rNum+1];
	a=new double[l];
#endif
	B=new double*[rNum+1];
	dBr=new double*[rNum+1];
	dBth=new double*[rNum+1];
	for(int i=0;i<rNum+1;i++){
		B[i]=new double[thNum];
		dBr[i]=new double[thNum];
		dBth[i]=new double[thNum];
#ifndef TOROIDAL
		A[i]=new double[thNum];
		Aaux[i]=new double[thNum];
		gsA[i]=new double[thNum];
#endif
	}

	//Get minimun radius from the value set in initial.cpp
	rmin=initial::rmin;

	//Solve size of steps in the radial and angular direction
	dr=(1.0-rmin)/(rNum-1);
	dth=Pi/(thNum-1);

	//Store the initial values of physical quantities
	for(int i=1;i<rNum;i++){
		r=rmin+i*dr;
		for(int j=1;j<thNum-1;j++){
			th=j*dth;
#ifndef SC
			B[i][j]=initial::B(r,th);
#else
			B[i][j]=initial::B(r-dr/2,th);
#endif
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
#ifndef SC
		B[0][j]=B[rNum-1][j]=0;
		B[rNum][j]=-B[rNum-2][j];
#else
		B[rNum-1][j]=B[rNum-2][j]/3;
		B[rNum][j]=-B[rNum-1][j];
#endif
	}
#ifndef TOROIDAL
	//Set boundary condition for poloidal fields
	for(int i=0;i<rNum;i++){
		A[i][0]=A[i][thNum-1]=0;
	}
	for(int j=0;j<thNum;j++){
		A[0][j]=0;
	}
#endif

	//Solve common values repeated multiple times during the simulation
	solve_repeated_values();

#ifndef TOROIDAL
	//Adjust alpha boundary conditions to fit smoothly to a poloidal field outside the star
	solve_A_boundary();
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
	hall_term_A=new double*[rNum+1];
	hall_rflux=new double*[rNum+1];
	hall_thflux=new double*[rNum+1];
#endif
	res_term_A=new double*[rNum+1];
	res_rflux=new double*[rNum+1];
	res_thflux=new double*[rNum+1];
	sines=new double[thNum];
	cotans=new double[thNum];
	for(int i=0;i<rNum+1;i++){
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
		cotans[j]=cos(th)/sin(th);
	}
	for(int i=0;i<rNum;i++){
		r=rmin+i*dr;
		for(int j=0;j<thNum-1;j++){
			th=j*dth;
#ifndef PUREOHM
			hall_term_A[i][j]=dt*sines[j]*initial::chi(r,th)/4/dr/dth;
#ifndef SC
			hall_rflux[i][j] = dt*initial::chi(r+dr/2,th)/8.0/dr/dth;
			hall_thflux[i][j]=-dt*initial::chi(r,th+dth/2)/8.0/dr/dth;
#else
			hall_rflux[i][j] = dt*initial::chi(r,th)/8.0/dr/dth;
			hall_thflux[i][j]=-dt*initial::chi(r-dr/2,th+dth/2)/8.0/dr/dth;
#endif
#endif
			res_term_A[i][j]=dt*thtd*initial::eta(r,th);
#ifndef SC
			res_rflux[i][j]  = dt*thtd*initial::eta(r+dr/2,th)/sin(th)/dr/dr;
			res_thflux[i][j] = dt*thtd*initial::eta(r,th+dth/2)/r/r/sin(th+dth/2.0)/dth/dth;
#else
			res_rflux[i][j]  = dt*thtd*initial::eta(r,th)/sin(th)/dr/dr;
			res_thflux[i][j] = dt*thtd*initial::eta(r-dr/2,th+dth/2)/r/r/sin(th+dth/2.0)/dth/dth;
#endif
		}
	}
	//Solve common terms involved in resolution of the multipole fit outside the star
#ifndef TOROIDAL
	boundary_factors1=new double[l];
	boundary_factors2=new double[l];
	legendre_comb=new double*[l];
	legendre_comb_mid=new double*[l];
	for(int n=0;n<l;n++){
		boundary_factors1[n]=(n+1.0)/(n+2.0)*sqrt(Pi*(2*n+3))/2*dth;
		boundary_factors2[n]=-2*dr*(n+1)*sqrt((2*n+3)/(4*Pi));
		legendre_comb[n]=new double[thNum];
		legendre_comb_mid[n]=new double[thNum];
		for(int j=0;j<thNum-1;j++){
			th=j*dth;
			legendre_comb[n][j]=(cos(th)*gsl_sf_legendre_Pl(n+1,cos(th))-gsl_sf_legendre_Pl(n,cos(th)));
			legendre_comb_mid[n][j]=(cos(th+dth*0.5)*gsl_sf_legendre_Pl(n+1,cos(th+dth*0.5))-gsl_sf_legendre_Pl(n,cos(th+dth*0.5)))/sin(th+dth*0.5);
		}
	}
#endif
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
	time_t t1, t2;
	time(&t1);
	for(int k=0;k<=tNum;k++){
		//Log data if k is multiple of plotSteps
		if(k%plotSteps==0){
			io::report_progress(k);
			t=k*dt;
			//Log integrated quantities
			io::log_integrals_file(t,solve_integrals());
			//Log complete profiles for A and B
			io::log_field(k);
			//No need to keep simulating if no output will be produced in next steps
			if(k+plotSteps>tNum){
				break;
			}
		}
#ifndef TOROIDAL
		//Update poloidal field function
		solve_new_A();
#endif
		//Update toroidal field function
		solve_new_B();
		//Pass values from auxiliary array Aaux and compute change in B from fluxes in dBr and dBth
		int exit=0;
		#pragma omp parallel for collapse(2)
		for(int i=1;i<rNum;i++){
			for(int j=1;j<thNum-1;j++){
				//Check for blowups, exit program if that happens
				if(isinf(dBr[i][j])||isinf(dBth[i][j])
#ifndef TOROIDAL
						||isinf(Aaux[i][j])
#endif
						){
					io::report_blowup(k,i,j);
					exit=1;
				}
				B[i][j]+=(dBr[i][j]-dBr[i-1][j]+dBth[i][j]-dBth[i][j-1])*sines[j];
#ifndef TOROIDAL
				A[i][j]=Aaux[i][j];
#endif
			}
		}
		if(exit)
			return 1;
#ifndef TOROIDAL
		//Fix beta value just outside the star so the numerical radial  derivative at the surface corresponds
		//to solving it backwards (i.e. using only the point at the surface and the one just below). Not required
		//for purely toroidal runs.
#ifndef SC
		#pragma omp parallel for
		for(int j=0;j<thNum;j++){
			B[0][j]=B[rNum-1][j]=0;
			B[rNum-2][j]=B[rNum-3][j]/2;
			B[1][j]=B[2][j]/2;
			B[rNum][j]=-B[rNum-2][j];
		}
#else
		#pragma omp parallel for
		for(int j=0;j<thNum;j++){
			B[rNum-1][j]=B[rNum-2][j]/3;
			B[rNum][j]=-B[rNum-1][j];
		}
#endif
		//Solve boundary condition of the surface of the star for alpha
		solve_A_boundary();
#endif
	}
	time(&t2);
	std::cout << std::endl << std::endl << "time: " << difftime(t2,t1) << std::endl;

	//Close file where integrated quantities are logged
	io::close_integrals_file();
	return 0;
}		/* -----  end of function simulate  ----- */

#ifndef TOROIDAL
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  solve_new_A
 *  Description:  
 * =====================================================================================
 */
	void
solve_new_A ()
{
	#pragma omp parallel for collapse(2) private(r,th)
	for(int i=1;i<rNum;i++){
		for(int j=1;j<thNum-1;j++){
			r=rmin+i*dr;
			th=j*dth;
			//Solve Grad-Shafranov operator at point
			gsA[i][j]=(A[i+1][j]+A[i-1][j]-2*A[i][j])/dr/dr+1/r/r*(A[i][j+1]+A[i][j-1]-2*A[i][j])/dth/dth-1/r/r*cotans[j]*(A[i][j+1]-A[i][j-1])/2/dth;
			//Evolve poloidal field function at point
			Aaux[i][j]=A[i][j]+res_term_A[i][j]*gsA[i][j];
#ifndef PUREOHM
#ifndef SC
			Aaux[i][j]+= ((B[i][j+1]-B[i][j-1])*(A[i+1][j]-A[i-1][j])
						 -(B[i+1][j]-B[i-1][j])*(A[i][j+1]-A[i][j-1]))*hall_term_A[i][j];
#else
			Aaux[i][j]+= ((B[i+1][j+1]+B[i][j+1]-B[i+1][j-1]-B[i][j-1])*(A[i+1][j]-A[i-1][j])/2
						 -(B[i+1][j]-B[i][j])*(A[i][j+1]-A[i][j-1])*2)*hall_term_A[i][j];
#endif
#endif
		}
	}
	return;
}		/* -----  end of function solve_new_A  ----- */
#endif

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  solve_new_B
 *  Description:  
 * =====================================================================================
 */
	void
solve_new_B ()
{
	#pragma omp parallel for collapse(2)
	for(int i=0;i<rNum-1;i++){
		for(int j=0;j<thNum-1;j++){
			//Solve radial fluxes on B/sin(th)
			if(j!=0){
				dBr[i][j]=res_rflux[i][j]*(B[i+1][j]-B[i][j]);
#ifndef PUREOHM
				dBr[i][j]+=hall_rflux[i][j]
					*(B[i][j]+B[i+1][j])
					*(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1]);
#ifndef TOROIDAL
#ifndef SC
				dBr[i][j]+=hall_rflux[i][j]
					*(gsA[i][j]+gsA[i+1][j])
					*(A[i][j+1]+A[i+1][j+1]-A[i][j-1]-A[i+1][j-1]);
#else
				if(i!=0){
					dBr[i][j]+=hall_rflux[i][j]
						*gsA[i][j]*2
						*(A[i][j+1]-A[i][j-1])*2;
				}else{
					dBr[i][j]=0;
				}
#endif
#endif
#endif
			}
			//Solve theta fluxes on B/sin(th)
#ifndef SC
			if(i!=0){
				dBth[i][j]=res_thflux[i][j]*(B[i][j+1]-B[i][j]);
#ifndef PUREOHM
				dBth[i][j]+=hall_thflux[i][j]
					*(B[i][j]+B[i][j+1])
					*(B[i+1][j]+B[i+1][j+1]-B[i-1][j]-B[i-1][j+1]);
#ifndef TOROIDAL
				dBth[i][j]+=hall_thflux[i][j]
					*(gsA[i][j]+gsA[i][j+1])
					*(A[i+1][j]+A[i+1][j+1]-A[i-1][j]-A[i-1][j+1]);
#endif
#endif
			}
#else
			if(i==0){
				dBth[i][j]=0;
			}else if(i==1){
				dBth[i][j]=res_thflux[i][j]*(B[i][j+1]-B[i][j]);
#ifndef PUREOHM
				dBth[i][j]+=hall_thflux[i][j]
					*(B[i][j]+B[i][j+1])
					*(B[i+1][j]+B[i+1][j+1]-B[i][j]-B[i][j+1])*2;
#ifndef TOROIDAL
				dBth[i][j]+=hall_thflux[i][j]
					*(gsA[i][j]+gsA[i][j+1])
					*(A[i+1][j]+A[i+1][j+1]-A[i-1][j]-A[i-1][j+1]);
#endif
#endif
			}else{
				dBth[i][j]=res_thflux[i][j]*(B[i][j+1]-B[i][j]);
#ifndef PUREOHM
				dBth[i][j]+=hall_thflux[i][j]
					*(B[i][j]+B[i][j+1])
					*(B[i+1][j]+B[i+1][j+1]-B[i-1][j]-B[i-1][j+1]);
#ifndef TOROIDAL
				dBth[i][j]+=hall_thflux[i][j]
					*(gsA[i][j]+gsA[i][j+1]+gsA[i-1][j]+gsA[i-1][j+1])/2
					*(A[i+1][j]+A[i+1][j+1]-A[i][j]-A[i][j+1])*2;
#endif
#endif
			}

#endif
		}
	}
	return;
}		/* -----  end of function solve_new_B  ----- */

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
	//When the poloidal field is considered, the individual energy of each multipole outside the star,
	//and the sum of all the multipoles is also given.
	double *integrals=new double[4+l];
	for(int n=0;n<4+l;n++){
		integrals[n]=0;
	}
#endif

	//Solve quantities integrated over the volume of the star
	for(int i=0;i<rNum-1;i++){
#ifndef TOROIDAL
		r=(0.5+i)*dr+rmin;
#endif
		for(int j=0;j<thNum-1;j++){
			th=(0.5+j)*dth;
			//Toroidal flux
			integrals[0]+=(B[i][j]+B[i+1][j]+B[i][j+1]+B[i+1][j+1])/4/sin(th);
			//Toroidal energy
			integrals[1]+=pow((B[i][j]+B[i+1][j]+B[i][j+1]+B[i+1][j+1])/4,2)/sin(th);
#ifndef TOROIDAL
			//Poloidal energy
			integrals[2]+=pow((A[i+1][j]-A[i][j]+A[i+1][j+1]-A[i][j+1])/2/dr,2)/sin(th)+pow((A[i][j+1]-A[i][j]+A[i+1][j+1]-A[i+1][j])/2/dth,2)/r/r/sin(th);
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
	//Solve external poloidal energy using the coefficients of the expansion outside the star
	//the coefficients "a" are missing some terms to represent the ones used in the notes
	for(int n=0;n<l;n++){
		//integrals[4+n]=pow(a[n],2)*(n+2)/(8*Pi);
		integrals[4+n]=pow(a[n],2)*(n+2)/(8*Pi);
		integrals[3]+=integrals[4+n];
	}
#endif
	return integrals;
}		/* -----  end of function solve_integrals  ----- */

#ifndef TOROIDAL
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  solve_A_boundary
 *  Description:  Solves the boundary values of A such that the field connects smoothly to a potential field outside the star
 * =====================================================================================
 */
	void
solve_A_boundary ( )
{
	//Solve integrated coefficients a_m.
	#pragma omp parallel for
	for(int n=0;n<l;n++){
		a[n]=0;
		for(int j=0;j<thNum-1;j++){
			//th=(j+0.5)*dth;
			//a[n]+=(A[rNum-1][j]+A[rNum-1][j+1])/2*dth*(cos(th)*gsl_sf_legendre_Pl(n+1,cos(th))-gsl_sf_legendre_Pl(n,cos(th)))/sin(th);
			a[n]+=(A[rNum-1][j]+A[rNum-1][j+1])*legendre_comb_mid[n][j]*boundary_factors1[n];
			//a[n]+=(A[rNum-1][j]+A[rNum-1][j+1])/2*dth*legendre_comb_mid[n][j];
			//a[n]+=(A[rNum-1][j]+A[rNum-1][j+1])/2*dth*legendre_comb_mid[n][j];
		}
		//a[n]=a[n]*boundary_factors1[n];
		//a[n]=a[n]*(n+1)/(n+2)*sqrt(Pi*(2*n+3));
	}
	//Permorm sumation
	#pragma omp parallel for
	for(int j=1;j<thNum-1;j++){
		A[rNum][j]=A[rNum-2][j];
		for(int n=0;n<l;n++){
			//th=j*dth;
			//A[rNum][j]+=-2*dr*a[n]*(n+1)*sqrt((2*n+3)/(4*Pi))*(cos(th)*gsl_sf_legendre_Pl(n+1,cos(th))-gsl_sf_legendre_Pl(n,cos(th)));
			//A[rNum][j]+=-2*dr*a[n]*(n+1)*sqrt((2*n+3)/(4*Pi))*legendre_comb[n][j];
			A[rNum][j]+=boundary_factors2[n]*a[n]*legendre_comb[n][j];
		}
	}
	return;
}		/* -----  end of function solve_A_boundary  ----- */
#endif

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  release_memory
 *  Description:  Destroy all used pointers to assure there are no memory leaks. The integer
 *  argument is used to specify the context in which the memory is released (abnormal
 *  termination, or normal one).
 * =====================================================================================
 */
	void
release_memory ( int info )
{
	for(int i=0;i<rNum+1;i++){
		delete[] B[i];
		B[i]=NULL;
		delete[] dBr[i];
		dBr[i]=NULL;
		delete[] dBth[i];
		dBth[i]=NULL;
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
		delete[] Aaux[i];
		Aaux[i]=NULL;
		delete[] gsA[i];
		gsA[i]=NULL;
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
	delete dBr;
	dBr=NULL;
	delete dBth;
	dBth=NULL;
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
	delete Aaux;
	Aaux=NULL;
	delete gsA;
	gsA=NULL;
	delete res_term_A;
	res_term_A=NULL;
#ifndef PUREOHM
	delete hall_term_A;
	hall_term_A=NULL;
#endif
	delete sines;
	sines=NULL;
	delete cotans;
	cotans=NULL;
	delete[] a;
	a=NULL;
	delete[] boundary_factors1;
	boundary_factors1=NULL;
	delete[] boundary_factors2;
	boundary_factors2=NULL;
	for(int n=0;n<l;n++){
		delete[] legendre_comb[n];
		legendre_comb[n]=NULL;
	}
	delete[] legendre_comb;
	legendre_comb=NULL;
#endif
	io::report_completion(info);
}		/* -----  end of function release_memory  ----- */

}		/* -----  end of namespace sim  ----- */
