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
//values given as cli arguments.
int rNum;
int thNum;
double dt;
int tNum;
int plotSteps;
double thtd;
//physical values of functions defining the magnetic field.
double **B;
#ifndef TOROIDAL
double **A;
#endif
//physical values that describe the structure of the star.
double **eta;
#ifndef PUREOHM
double **chi;
#endif
//arrays that contains the precalculated quantities neccesary to solve fluxes for B
#ifndef PUREOHM
double **hall_rflux;
double **hall_thflux;
#endif
double **res_rflux;
double **res_thflux;
//sines precalculated at each point in the grid
double *sines;
//minimun radius of the shell containing the magnetic field.
double rmin;
//size of spatial steps.
double dr;
double dth;

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
#endif
	B=new double*[rNum];
	eta=new double*[rNum];
#ifndef PUREOHM
	chi=new double*[rNum];
	hall_rflux=new double*[rNum];
	hall_thflux=new double*[rNum];
#endif
	res_rflux=new double*[rNum];
	res_thflux=new double*[rNum];
	sines=new double[thNum];
	for(int i=0;i<rNum;i++){
#ifndef TOROIDAL
		A[i]=new double[thNum];
#endif
		B[i]=new double[thNum];
		eta[i]=new double[thNum];
#ifndef PUREOHM
		chi[i]=new double[thNum];
		hall_rflux[i]=new double[thNum];
		hall_thflux[i]=new double[thNum];
#endif
		res_rflux[i]=new double[thNum];
		res_thflux[i]=new double[thNum];
	}

	//get minimun radius from the value set in initial.cpp
	rmin=initial::rmin;

	//solve size of steps in the radial and angular direction
	dr=(1.0-rmin)/(rNum-1);
	double const Pi=4*atan(1);
	dth=Pi/(thNum-1);

	//store the initial values of physical quantities
	double r,th;
	for(int i=1;i<rNum-1;i++){
		r=rmin+i*dr;
		for(int j=1;j<thNum-1;j++){
			th=j*dth;
#ifndef TOROIDAL
			A[i][j]=initial::A(r,th);
#endif
			B[i][j]=initial::B(r,th);
			eta[i][j]=initial::eta(r,th);
#ifndef PUREOHM
			chi[i][j]=initial::chi(r,th);
#endif
		}
	}

	//set boundary conditions
	for(int i=0;i<rNum;i++){
#ifndef TOROIDAL
		A[i][0]=A[i][thNum-1]=0;
#endif
		B[i][0]=B[i][thNum-1]=0;
	}
	for(int j=0;j<thNum;j++){
#ifndef TOROIDAL
		A[0][j]=A[rNum-1][j]=0;
#endif
		B[0][j]=B[rNum-1][j]=0;
	}
#ifndef TOROIDAL
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
	//Solve common terms involved in the calculation of the toroidal field.
	double r,th;
	for(int i=0;i<rNum-1;i++){
		r=rmin+i*dr;
		for(int j=0;j<thNum-1;j++){
			th=j*dth;
#ifndef PUREOHM
			hall_rflux[i][j] = dt*initial::chi(r+dr/2,th)/8.0/dr/dth;
			hall_thflux[i][j]=-dt*initial::chi(r,th+dth/2)/8.0/dr/dth;
#endif
			res_rflux[i][j]  = dt*thtd*initial::eta(r+dr/2,th)/sin(th)/dr/dr;
			res_thflux[i][j] = dt*thtd*initial::eta(r,th+dth/2)/r/r/sin(th+dth/2.0)/dth/dth;
		}
	}
	for(int j=1;j<thNum-1;j++){
		th=j*dth;
		sines[j]=sin(th);
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
	for(int k=0;k<=tNum;k++){
		//log data if k is multiple of plotSteps
		if(k%plotSteps==0){
			t=k*dt;
			//log integrated quantities
			io::log_integrals_file(t,solve_integrals());
			//log complete profiles for A and B
			io::log_field(k);
			//No need to keep simulating if no output will be produced in next steps
			if(k+plotSteps>tNum){
				break;
			}
			std::cout << k << "/" << tNum << std::endl;
		}
#ifndef TOROIDAL
		//Solve grad shafranov operator acting on A on all grid (except boundaries)
		for(int i=1;i<rNum-1;i++){
			double r=rmin+i*dr;
			for(int j=1;j<thNum-1;j++){
				double th=j*dth;
				gsA[i][j]=(A[i+1][j]+A[i-1][j]-2*A[i][j])/dr/dr+1/r/r*(A[i][j+1]+A[i][j-1]-2*A[i][j])/dth/dth-1/r/r*cos(th)/sin(th)*(A[i][j+1]-A[i][j-1])/2/dth;
			}
		}

		//update poloidal field function
		for(int i=1;i<rNum-1;i++){
			for(int j=1;j<thNum-1;j++){
#ifndef PUREOHM
//				//cut timestep by half if displacement is too large
//				double dispr=dt*sines[j]*chi[i][j]*(B[i][j+1]-B[i][j-1])/2/dth;
//				double dispth=-dt*sines[j]*chi[i][j]*(B[i+1][j]-B[i-1][j])/2/dr;
//				double r=rmin+i*dr+dispr;
//				double th=j*dth-dispth;
//				//solve 4 point grid to interpolate (or perhaps extrapolate)
//				int imoved=(r-rmin)/dr;
//				int jmoved=th/dth;
//				if(imoved==0){
//					imoved=1;
//				}
//				if(jmoved==0){
//					jmoved=1;
//				}
//				if(imoved==rNum-2){
//					imoved=rNum-3;
//				}
//				if(jmoved==thNum-2){
//					jmoved=thNum-3;
//				}
//				r-=(rmin+imoved*dr);
//				th-=(jmoved*dth);
//				//using the four points, approximate the grad-shafranov of alpha
//				double a=(gsA[imoved+1][jmoved]-gsA[imoved][jmoved])/dr;
//				double b=(gsA[imoved][jmoved+1]-gsA[imoved][jmoved])/dth;
//				double c=(gsA[imoved][jmoved]+gsA[imoved+1][jmoved+1]
//						-gsA[imoved+1][jmoved]-gsA[imoved][jmoved+1])/dr/dth;
//
//				double a2=(A[imoved+1][jmoved]-A[imoved][jmoved])/dr;
//				double b2=(A[imoved][jmoved+1]-A[imoved][jmoved])/dth;
//				double c2=(A[imoved][jmoved]+A[imoved+1][jmoved+1]
//						-A[imoved+1][jmoved]-A[imoved][jmoved+1])/dr/dth;
//				Aaux[i][j]=(A[imoved][jmoved]+a2*r+b2*th+c2*r*th)
//					+dt*thtd*eta[i][j]*(gsA[imoved][jmoved]+a*r+b*th+c*r*th);
				Aaux[i][j]=A[i][j]
					+dt*sines[j]*chi[i][j]*(B[i][j+1]-B[i][j-1])/2/dth*(A[i+1][j]-A[i-1][j])/2/dr
					-dt*sines[j]*chi[i][j]*(B[i+1][j]-B[i-1][j])/2/dr*(A[i][j+1]-A[i][j-1])/2/dth
					+dt*thtd*eta[i][j]*gsA[i][j];
#else
				Aaux[i][j]=A[i][j]+dt*thtd*eta[i][j]*gsA[i][j];
#endif
			}
		}
#endif

		//Update toroidal field function
		for(int i=0;i<rNum-1;i++){
			for(int j=0;j<thNum-1;j++){
				if(i==0&&j==0)
					continue;
				double dBr=0;
				double dBth=0;
				if(j!=0){
#ifndef PUREOHM
					dBr+=hall_rflux[i][j]
						*(B[i+1][j]+B[i][j])
						*(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1]);
#endif
					dBr+=res_rflux[i][j]*(B[i+1][j]-B[i][j]);
#ifndef TOROIDAL
#ifndef PUREOHM
					if(i!=0){
						dBr+=dt*initial::chi(rmin+i*dr+dr/2,j*dth)
							*(gsA[i][j]+gsA[i+1][j])
							*(A[i][j+1]+A[i+1][j+1]-A[i][j-1]-A[i+1][j-1])
							/8/dr/dth;
					}else{
						dBr+=dt*initial::chi(rmin+i*dr+dr/2,j*dth)
							*(3*gsA[i+1][j]-gsA[i+2][j])
							*(A[i][j+1]+A[i+1][j+1]-A[i][j-1]-A[i+1][j-1])
							/8/dr/dth;
					}
#endif
#endif
				}
				if(i!=0){
#ifndef PUREOHM
					dBth+=hall_thflux[i][j]
						*(B[i][j+1]+B[i][j])
						*(B[i+1][j]+B[i+1][j+1]-B[i-1][j]-B[i-1][j+1]);
#endif
					dBth+=res_thflux[i][j]*(B[i][j+1]-B[i][j]);
#ifndef TOROIDAL
#ifndef PUREOHM
					if(j!=0){
						dBth+=-dt*initial::chi(rmin+i*dr,j*dth+dth/2)
							*(gsA[i][j]+gsA[i][j+1])
							*(A[i+1][j]+A[i+1][j+1]-A[i-1][j]-A[i-1][j+1])
							/8/dr/dth;
					}else{
						dBth+=-dt*initial::chi(rmin+i*dr,j*dth+dth/2)
							*(3*gsA[i][j+1]-gsA[i][j+2])
							*(A[i+1][j]+A[i+1][j+1]-A[i-1][j]-A[i-1][j+1])
							/8/dr/dth;
					}
#endif
#endif
				}
				Baux[i][j]+=(dBr+dBth)*sines[j];
				Baux[i+1][j]=B[i+1][j]-dBr*sines[j];
				Baux[i][j+1]-=dBth*sines[j+1];
			}
		}
#ifndef TOROIDAL
		solve_A_boundary();
#endif
		for(int i=1;i<rNum;i++){
			for(int j=1;j<thNum-1;j++){
				//check for blowups, exit program if that happens
				if(isinf(Baux[i][j])
#ifndef TOROIDAL
						||isinf(Aaux[i][j])
#endif
						){
					io::report_blowup(k,i,j);
					return 1;
				}
				if(i!=rNum-1)
					B[i][j]=Baux[i][j];
#ifndef TOROIDAL
				A[i][j]=Aaux[i][j];
#endif
			}
		}
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
#ifdef TOROIDAL
	double *integrals=new double[2];
	integrals[0]=integrals[1]=0;
#else
	double *integrals=new double[3];
	integrals[0]=integrals[1]=integrals[2]=0;
#endif
	double r;
	for(int i=1;i<rNum-1;i++){
		r=i*dr+rmin;
		for(int j=1;j<thNum-1;j++){
			//Toroidal flux
			integrals[0]+=B[i][j]/sines[j];
			//Toroidal energy
			integrals[1]+=pow(B[i][j],2)/r;
#ifndef TOROIDAL
			//Poloidal energy
			integrals[2]+=pow((A[i+1][j]-A[i-1][j])/2/dr,2)/sines[j]+pow((A[i][j+1]-A[i][j-1])/2/dth,2)/r/r/sines[j];
#endif
		}
	}
	integrals[0]=integrals[0]*dr*dth;
	integrals[1]=integrals[1]*dr*dth;
#ifndef TOROIDAL
	integrals[2]=integrals[2]*dr*dth;
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
	void
solve_A_boundary ( )
{
	//maximun l for the spherical harmonics expansion
	int l=10;
	//solve the required legendre polynomial on all required points
	double P[l+1][thNum];
	for(int i=0;i<=l;i++){
		for(int j=1;j<thNum-1;j++){
			P[i][j]=gsl_sf_legendre_Pl(i,cos(j*dth));
		}
	}
	//solve integrated coefficients
	double a[l];
	for(int i=0;i<l;i++){
		a[i]=0;
		for(int j=1;j<thNum-1;j++){
			a[i]+=(A[rNum-2][j+1]-A[rNum-2][j-1])/2*P[i+1][j]*dr*(i+1)*(2*i+3)/(i+1);
		}
	}
	//set boundary to A[rNum-3][j] to start summation
	for(int j=1;j<thNum-1;j++){
		A[rNum-1][j]=A[rNum-3][j];
	}
	//permorm sumation
	for(int i=0;i<l;i++){
		for(int j=1;j<thNum-1;j++){
			A[rNum-1][j]+=a[i]*(cos(j*dth)*P[i+1][j]-P[i][j]);
		}
	}
	return;
}		/* -----  end of function solve_A_boundary  ----- */
#endif

}		/* -----  end of namespace sim  ----- */
