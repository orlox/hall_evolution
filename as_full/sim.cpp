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
double factor;
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
#endif
#ifdef SC
	//products of n, eta and thtd rigth at the inner boundary, required to solve superconductor boundary conditions
	double *sc_factors;
#endif
//Size of spatial steps.
double dr,dth;
//Size of timestep, set automatically by code
double dt=1;
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
		A=new double*[rNum+2];
		Aaux=new double*[rNum+2];
		gsA=new double*[rNum+2];
		a=new double[l];
	#endif
	B=new double*[rNum+2];
	dBr=new double*[rNum+2];
	dBth=new double*[rNum+2];
	for(int i=0;i<rNum+2;i++){
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
	for(int i=1;i<rNum+1;i++){
		r=rmin+(i-1)*dr;
		for(int j=1;j<thNum-1;j++){
			th=j*dth;
			B[i][j]=initial::B(r,th);
			#ifndef TOROIDAL
				A[i][j]=initial::A(r,th);
			#endif
		}
	}

	//Set boundary conditions for toroidal field at the inner and outer boundary
	for(int j=0;j<thNum;j++){
		B[rNum][j]=0;
		B[rNum+1][j]=-B[rNum-1][j];
	}
	//Set boundary conditions for toroidal and poloidal field at the axis
	for(int i=0;i<rNum+2;i++){
		B[i][0]=B[i][thNum-1]=0;
		#ifndef TOROIDAL
			A[i][0]=A[i][thNum-1]=0;
		#endif
	}

	//Solve common values repeated multiple times during the simulation
	solve_repeated_values();

	#ifdef SC
	//Solve boundary conditions for beta. This is only required if superconductor boundary conditions are used,
	//or the simulation is not purely toroidal.
	solve_B_boundary();
	#else
		#ifndef TOROIDAL
			solve_B_boundary();
		#endif
	#endif
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
		hall_term_A=new double*[rNum+2];
		hall_rflux=new double*[rNum+2];
		hall_thflux=new double*[rNum+2];
	#endif
	res_term_A=new double*[rNum+2];
	res_rflux=new double*[rNum+2];
	res_thflux=new double*[rNum+2];
	sines=new double[thNum];
	cotans=new double[thNum];
	#ifdef SC
		sc_factors=new double[thNum];
	#endif
	for(int i=0;i<rNum+2;i++){
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
		#ifdef SC
			sc_factors[j]=initial::eta(rmin,th)*initial::n(rmin,th)*thtd;
		#endif
	}
	for(int i=0;i<rNum+2;i++){
		r=rmin+(i-1)*dr;
		for(int j=0;j<thNum-1;j++){
			th=j*dth;
			#ifndef PUREOHM
				hall_term_A[i][j]=sines[j]*initial::chi(r,th)/4/dr/dth;
				hall_rflux[i][j] = initial::chi(r+dr/2,th)/8.0/dr/dth;
				hall_thflux[i][j]=-initial::chi(r,th+dth/2)/8.0/dr/dth;
			#endif
			res_term_A[i][j]=thtd*initial::eta(r,th);
			res_rflux[i][j]  = thtd*initial::eta(r+dr/2,th)/sin(th)/dr/dr;
			res_thflux[i][j] = thtd*initial::eta(r,th+dth/2)/r/r/sin(th+dth/2.0)/dth/dth;
		}
	}
	#ifndef TOROIDAL
		//Solve common terms involved in resolution of the multipole fit outside the star
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
	double t=0;
	time_t t1, t2;
	time(&t1);
	for(int k=0;k<=tNum;k++){
		/////////////////////////////////////////////////////////////Update timestep///////////////////////////////////////////////////////////
		//////////////////////////////////////////////////THIS SECTION BREAKS NON-CONSTANT eta AND n///////////////////////////////////////////
		double newdt=100000000;
		double localnewdt=100000000;
		double temp=0;
		//Solve Grad-Shafranov operator at all points
		//#pragma omp parallel for collapse(2) private(r)
		for(int i=1;i<rNum+1;i++){
			for(int j=1;j<thNum-1;j++){
				r=rmin+(i-1)*dr;
				gsA[i][j]=(A[i+1][j]+A[i-1][j]-2*A[i][j])/dr/dr+1/r/r*(A[i][j+1]+A[i][j-1]-2*A[i][j])/dth/dth-1/r/r*cotans[j]*(A[i][j+1]-A[i][j-1])/2/dth;
			}
		}
		//Obtain critical timestep from Hall Drift
		//#pragma omp parallel for private(r,th,localnewdt,temp)
		for(int i=2;i<rNum+1;i++){
			for(int j=1;j<thNum-1;j++){
				r=rmin+(i-1)*dr;
				th=j*dth;
				//Solve critical timestep
				temp=dr/sqrt(pow(gsA[i][j]/(r*sines[j]),2)+pow((B[i+1][j]-B[i-1][j])/(2*dr*r*sines[j]),2)+pow((B[i][j+1]-B[i][j-1])/(2*dth*r*r*sines[j]),2));
				if(temp<localnewdt){
					localnewdt=temp;
				}
			}
			//#pragma omp critical
			if(localnewdt<newdt){
				newdt=localnewdt;
			}
		}
		//If Ohmic critical step is smaller, use that one
		if(newdt>dr*dr/thtd){
			newdt=dr*dr/thtd;
		}
		dt=factor*newdt;
		if(dt<0.0000000001){
			dt=0.0000000001;
		}
		/////////////////////////////////////////////////////////End Update Timestep/////////////////////////////////////////////////////////////////
		//Log data if k is multiple of plotSteps
		if(k%plotSteps==0){
			std::cout<<"dt:"<<dt<<std::endl;
			io::report_progress(k,t);
			//Log integrated quantities
			io::log_integrals_file(t,solve_integrals());
			//Log complete profiles for A and B
			io::log_field(k,t);
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
		for(int i=1;i<rNum+1;i++){
			for(int j=1;j<thNum-1;j++){
				//Check for blowups, exit program if that happens
				if(isinf(dBr[i][j])||isinf(dBth[i][j])||isnan(dBr[i][j])||isnan(dBth[i][j])
						#ifndef TOROIDAL
						||isinf(Aaux[i][j])||isnan(Aaux[i][j])
						#endif
						){
					if(exit!=1){
						io::report_blowup(k,i,j);
						exit=1;
					}
				}
				#ifndef SC
					B[i][j]+=(dBr[i][j]-dBr[i-1][j]+dBth[i][j]-dBth[i][j-1])*sines[j];
				#else
					if(i==1){
						B[i][j]+=(2*dBr[i][j]+dBth[i][j]-dBth[i][j-1])*sines[j];
					}else{
						B[i][j]+=(dBr[i][j]-dBr[i-1][j]+dBth[i][j]-dBth[i][j-1])*sines[j];
					}
				#endif
				#ifndef TOROIDAL
					A[i][j]=Aaux[i][j];
				#endif
			}
		}
		if(exit){//||t>2){
			return 1;
		}
		#ifdef SC
		//Solve boundary conditions for beta. This is only required if superconductor boundary conditions are used,
		//or the simulation is not purely toroidal.
		solve_B_boundary();
		#else
			#ifndef TOROIDAL
				solve_B_boundary();
			#endif
		#endif
		#ifndef TOROIDAL
		//Solve boundary conditions for alpha
		solve_A_boundary();
		#endif
		//update simulation time
		t+=dt;
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
	for(int i=1;i<rNum+1;i++){
		for(int j=1;j<thNum-1;j++){
			r=rmin+(i-1)*dr;
			th=j*dth;
			//Evolve poloidal field function at point
			if(i!=1){
				Aaux[i][j]=A[i][j]+dt*res_term_A[i][j]*gsA[i][j];
				#ifndef PUREOHM
					Aaux[i][j]+= dt*((B[i][j+1]-B[i][j-1])*(A[i+1][j]-A[i-1][j])
								 -(B[i+1][j]-B[i-1][j])*(A[i][j+1]-A[i][j-1]))*hall_term_A[i][j];
				#endif
			}
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
	for(int i=1;i<rNum;i++){
		for(int j=0;j<thNum-1;j++){
			//Solve radial fluxes on B/sin(th)
			if(j!=0){
				dBr[i][j]=dt*res_rflux[i][j]*(B[i+1][j]-B[i][j]);
				#ifndef PUREOHM
					dBr[i][j]+=dt*hall_rflux[i][j]
						*(B[i][j]+B[i+1][j])
						*(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1]);
					#ifndef TOROIDAL
						dBr[i][j]+=dt*hall_rflux[i][j]
							*(gsA[i][j]+gsA[i+1][j])
							*(A[i][j+1]+A[i+1][j+1]-A[i][j-1]-A[i+1][j-1]);
					#endif
				#endif
			}
			//Solve theta fluxes on B/sin(th)
			if(i!=0){
				dBth[i][j]=dt*res_thflux[i][j]*(B[i][j+1]-B[i][j]);
				#ifndef PUREOHM
					dBth[i][j]+=dt*hall_thflux[i][j]
						*(B[i][j]+B[i][j+1])
						*(B[i+1][j]+B[i+1][j+1]-B[i-1][j]-B[i-1][j+1]);
					#ifndef TOROIDAL
						dBth[i][j]+=dt*hall_thflux[i][j]
							*(gsA[i][j]+gsA[i][j+1])
							*(A[i+1][j]+A[i+1][j+1]-A[i-1][j]-A[i-1][j+1]);
					#endif
				#endif
			}
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
		double *integrals=new double[3];
		integrals[0]=integrals[1]=integrals[2]=0;
	#else
		//When the poloidal field is considered, the individual energy of each multipole outside the star,
		//and the sum of all the multipoles is also given.
		double *integrals=new double[5+l];
		for(int n=0;n<5+l;n++){
			integrals[n]=0;
		}
	#endif

	//Solve quantities integrated over the volume of the star
	for(int i=1;i<rNum;i++){
		#ifndef TOROIDAL
			r=(0.5+i-1)*dr+rmin;
		#endif
		for(int j=0;j<thNum-1;j++){
			th=(0.5+j)*dth;
			//Toroidal flux
			integrals[0]+=(B[i][j]+B[i+1][j]+B[i][j+1]+B[i+1][j+1])/4/sin(th);
            //E/(dE/dt)
			integrals[1]+=pow((B[i+1][j]-B[i][j]+B[i+1][j+1]-B[i][j+1])/2/dr,2)/sin(th)+pow((B[i][j+1]-B[i][j]+B[i+1][j+1]-B[i+1][j])/2/dth,2)/r/r/sin(th);
			//Toroidal energy
			integrals[2]+=pow((B[i][j]+B[i+1][j]+B[i][j+1]+B[i+1][j+1])/4,2)/sin(th);
			#ifndef TOROIDAL
                //E/(dE/dt)
			    integrals[1]+=pow((gsA[i][j]+gsA[i+1][j]+gsA[i][j+1]+gsA[i+1][j+1])/4,2)/sin(th);
				//Poloidal energy
				integrals[3]+=pow((A[i+1][j]-A[i][j]+A[i+1][j+1]-A[i][j+1])/2/dr,2)/sin(th)+pow((A[i][j+1]-A[i][j]+A[i+1][j+1]-A[i+1][j])/2/dth,2)/r/r/sin(th);
			#endif
		}
	}
	//Multiply by the area element on the r-th grid, and also by (2*Pi) and 1/(8*pi) the energies, where
	//the first factor is from integration over phi, and the second one is the factor present in the calculation
	//of the magnetic energy.
	integrals[0]=integrals[0]*dr*dth;
	integrals[1]=integrals[1]*dr*dth/2*thtd;
	integrals[2]=integrals[2]*dr*dth/4;
	#ifndef TOROIDAL
		integrals[3]=integrals[3]*dr*dth/4;
		//Solve external poloidal energy using the coefficients of the expansion outside the star
		//the coefficients "a" are missing some terms to represent the ones used in the notes
		for(int n=0;n<l;n++){
			//integrals[4+n]=pow(a[n],2)*(n+2)/(8*Pi);
			integrals[5+n]=pow(a[n],2)*(n+2)/(8*Pi);
			integrals[4]+=integrals[5+n];
		}
	#endif
    integrals[1]=(integrals[2]+integrals[3]+integrals[4])/integrals[1];
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
	//Solve integrated coefficients a_l.
	#pragma omp parallel for
	for(int n=0;n<l;n++){
		a[n]=0;
		for(int j=0;j<thNum-1;j++){
			a[n]+=(A[rNum][j]+A[rNum][j+1])*legendre_comb_mid[n][j]*boundary_factors1[n];
		}
	}
	//Fix A value at the inner boundary to zero and solve value of inner point to set superconductor boundary condition.
    //If not using SC boundary conditions, set A value such that the second radial derivative of A is zero at the boundary,
    //which is equivalent to gsA=0 there.
	//Then, perform summation to impose boundary condition on the surface.
	#pragma omp parallel for
	for(int j=1;j<thNum-1;j++){
		A[1][j]=0;
		#ifdef SC
			#ifndef PUREOHM
				A[0][j]=-A[2][j]*(sc_factors[j]*thtd/dr+(B[1][j+1]-B[1][j-1])/(pow(rmin,2)*sines[j]*4*dth))/(sc_factors[j]*thtd/dr-(B[1][j+1]-B[1][j-1])/(pow(rmin,2)*sines[j]*4*dth));
			#endif
        #else
            A[0][j]=-A[2][j];
		#endif
		A[rNum+1][j]=A[rNum-1][j];
		for(int n=0;n<l;n++){
			A[rNum+1][j]+=boundary_factors2[n]*a[n]*legendre_comb[n][j];
		}
	}
	return;
}		/* -----  end of function solve_A_boundary  ----- */
#endif

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  solve_B_boundary
 *  Description:  Apply boundary conditions for beta
 * =====================================================================================
 */
	void
solve_B_boundary ( )
{
	//Fix beta value just outside the star so the numerical radial  derivative at the surface corresponds
	//to solving it backwards (i.e. using only the point at the surface and the one just below). Also, Apply linear
	//interpolation of beta for the point just below the surface and just above the inner boundary (this last one
	//only in the case of zero boundary conditions at the center).
	//Interpolation at both the inner and outer boundary is unneccesary for purely toroidal runs, while on superconductor
	//runs interpolation in the inner boundary is unnecesary.
	#pragma omp parallel for
	for(int j=1;j<thNum-1;j++){
		#ifndef SC
			B[1][j]=0;
			B[2][j]=B[3][j]/2;
		#else
			#ifndef PUREOHM
				B[0][j]=2*dr/sc_factors[j]*(B[1][j]*(B[1][j+1]-B[1][j-1])/(2*dth*pow(rmin,2)*sines[j]))+B[2][j];
			#else
				B[0][j]=B[2][j];
			#endif
		#endif
		B[rNum][j]=0;
		//B[rNum-1][j]=B[rNum-2][j]/2;
		B[rNum+1][j]=-B[rNum-1][j];
	}
	return;
}		/* -----  end of function solve_B_boundary  ----- */

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
	for(int i=0;i<rNum+2;i++){
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
