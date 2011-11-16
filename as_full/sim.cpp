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
	A=new double*[2*rNum-1];
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
#ifndef TOROIDAL
	for(int i=0;i<2*rNum-1;i++){
		A[i]=new double[2*thNum-1];
	}
#endif

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
			B[i][j]=initial::B(r,th);
			eta[i][j]=initial::eta(r,th);
#ifndef PUREOHM
			chi[i][j]=initial::chi(r,th);
#endif
		}
	}
#ifndef TOROIDAL
	for(int i=1;i<2*rNum-2;i++){
		r=rmin+i*dr/2;
		for(int j=1;j<2*thNum-2;j++){
			th=j*dth/2;
			A[i][j]=initial::A(r,th);
		}
	}
#endif

	//set boundary conditions for toroidal field
	for(int i=0;i<rNum;i++){
		B[i][0]=B[i][thNum-1]=0;
	}
	for(int j=0;j<thNum;j++){
		B[0][j]=B[rNum-1][j]=0;
	}
#ifndef TOROIDAL
	for(int i=0;i<2*rNum-1;i++){
		A[i][0]=A[i][2*thNum-2]=0;
	}
	for(int j=0;j<2*thNum-1;j++){
		A[0][j]=A[2*rNum-2][j]=0;
	}
#endif

#ifndef TOROIDAL
#ifndef SIMPLE
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
	for(int i=0;i<rNum;i++){
		gsA[i][0]=gsA[i][thNum-1]=0;
	}
	for(int j=0;j<thNum;j++){
		gsA[0][j]=gsA[rNum-1][j]=0;
	}
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
				gsA[i][j]=(105.0/2.0*(-pow(r,2)+pow(r,4))*pow(sin(th),2));
			}
		}

		//update poloidal field function
		for(int i=1;i<rNum-1;i++){
			//double r=rmin+i*dr;
			for(int j=1;j<thNum-1;j++){
				//double th=i*dth;
				Aaux[i][j]=A[i][j];
#ifndef PUREOHM
				Aaux[i][j]=
					 dt*sines[j]*chi[i][j]*(B[i][j+1]-B[i][j-1])/2/dth*(A[i+1][j]-A[i-1][j])/2/dr
					-dt*sines[j]*chi[i][j]*(B[i+1][j]-B[i-1][j])/2/dr*(A[i][j+1]-A[i][j-1])/2/dth;
#endif
				Aaux[i][j]+=A[i][j]+dt*thtd*eta[i][j]*gsA[i][j];
			}
		}
#endif

		//Update toroidal field function
		for(int i=10;i<rNum-11;i++){
			double r=rmin+i*dr;
			for(int j=3;j<thNum-4;j++){
			double th=j*dth;
				if(i==0&&j==0)
					continue;
				double dBr=0;
				double dBth=0;
//				double dB1=0;
//				double dB2=0;
				if(j!=0){
					dBr+=thtd*initial::eta(r+dr/2,th)/sin(th)
						*(B[i+1][j]-B[i][j])/dr;

					dBr+=initial::chi(r+dr/2,th)
						*(B[i][j]+B[i+1][j])/2
						*(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1])/4/dth;

					dBr+=initial::chi(r+dr/2,th)
						*(gsA[i][j]+gsA[i+1][j])/2
						*(A[i][j+1]+A[i+1][j+1]-A[i][j-1]-A[i+1][j-1])/4/dth;

					dBr=dBr/dr*dt;
				}
				if(i!=0){
					dBth+=thtd*initial::eta(r,th+dth/2)/sin(th+dth/2)/pow(r,2)
						*(B[i][j+1]-B[i][j])/dth;

					dBth+=-initial::chi(r,th+dth/2)
						*(B[i][j]+B[i][j+1])/2
						*(B[i+1][j]+B[i+1][j+1]-B[i-1][j]-B[i-1][j+1])/4/dr;

					dBth+=-initial::chi(r,th+dth/2)
						*(gsA[i][j]+gsA[i][j+1])/2
						*(A[i+1][j]+A[i+1][j+1]-A[i-1][j]-A[i-1][j+1])/4/dr;

					dBth=dBth/dth*dt;
				}
//
//				if(i!=0&&j!=0){
//					double a2,a3,a4,a5,a6,a7,a8,a9;
//					double f1=A[i-1][j+1];
//					double f2=A[i][j+1];
//					double f3=A[i+1][j+1];
//					double f4=A[i-1][j];
//					double f5=A[i][j];
//					double f6=A[i+1][j];
//					double f7=A[i-1][j-1];
//					double f8=A[i][j-1];
//					double f9=A[i+1][j-1];
//					a2=-(f4-f6)/(2*dr);
//					a3=(f6-2*f5+f4)/(2*dr*dr);
//					a4=(f2-f8)/(2*dth);
//					a5=(f8-2*f5+f2)/(2*dth*dth);
//					a6=-(f9-f7-f3+f1)/(4*dr*dth);
//					a7=(-f9+2*f8-f7+f3-2*f2+f1)/(4*dr*dr*dth);
//					a8=-(-f9+f7+2*f6-2*f4-f3+f1)/(4*dr*dth*dth);
//					a9=(f9-2*f8+f7-2*f6+4*f5-2*f4+f3-2*f2+f1)/(4*dr*dr*dth*dth);
//					int caca=8;
//					dr=dr/caca;
//					dth=dth/caca;
//					int st=8;
//					double thaux=-dth/2+dth/2/st;
//					for(int n=0;n<st;n++){
//						dB1+=initial::chi(r+dr/2,th+thaux)
//							*(
//								2*a3+2*a7*thaux+2*a9*pow(thaux,2)
//								+1/pow(r+dr/2,2)*(2*a5+2*a8*(dr/2)+2*a9*pow(dr/2,2))
//								-1/pow(r+dr/2,2)*cos(th+thaux)/sin(th+thaux)*(a4+2*a5*thaux+a6*(dr/2)+a7*pow(dr/2,2)+2*a8*(dr/2)*thaux+2*a9*pow(dr/2,2)*thaux)
//							 )
//							*(a4+2*a5*thaux+a6*(dr/2)+a7*pow(dr/2,2)+2*a8*(dr/2)*thaux+2*a9*pow(dr/2,2)*thaux)
//							/dr/st;
//						dB1-=initial::chi(r-dr/2,th+thaux)
//							*(
//								2*a3+2*a7*thaux+2*a9*pow(thaux,2)
//								+1/pow(r-dr/2,2)*(2*a5+2*a8*(-dr/2)+2*a9*pow(-dr/2,2))
//								-1/pow(r-dr/2,2)*cos(th+thaux)/sin(th+thaux)*(a4+2*a5*thaux+a6*(-dr/2)+a7*pow(-dr/2,2)+2*a8*(-dr/2)*thaux+2*a9*pow(-dr/2,2)*thaux)
//							 )
//							*(a4+2*a5*thaux+a6*(-dr/2)+a7*pow(-dr/2,2)+2*a8*(-dr/2)*thaux+2*a9*pow(-dr/2,2)*thaux)
//							/dr/st;
//
//						thaux+=dth/st;
//					}
//					double raux=-dr/2+dr/2/st;
//					for(int n=0;n<st;n++){
//						dB2+=-initial::chi(r+raux,th+dth/2)*
//							(
//								2*a3+2*a7*(dth/2)+2*a9*pow(dth/2,2)
//								+1/pow(r+raux,2)*(2*a5+2*a8*raux+2*a9*pow(raux,2))
//								-1/pow(r+raux,2)*cos(th+dth/2)/sin(th+dth/2)*(a4+2*a5*(dth/2)+a6*raux+a7*pow(raux,2)+2*a8*raux*(dth/2)+2*a9*pow(raux,2)*(dth/2))
//							 )
//							*(a2+2*a3*raux+a6*(dth/2)+2*a7*raux*(dth/2)+a8*pow(dth/2,2)+2*a9*raux*pow(dth/2,2))
//							/dth/st;
//						dB2-=-initial::chi(r+raux,th-dth/2)*
//							(
//								2*a3+2*a7*(-dth/2)+2*a9*pow(-dth/2,2)
//								+1/pow(r+raux,2)*(2*a5+2*a8*raux+2*a9*pow(raux,2))
//								-1/pow(r+raux,2)*cos(th-dth/2)/sin(th-dth/2)*(a4+2*a5*(-dth/2)+a6*raux+a7*pow(raux,2)+2*a8*raux*(-dth/2)+2*a9*pow(raux,2)*(-dth/2))
//							 )
//							*(a2+2*a3*raux+a6*(-dth/2)+2*a7*raux*(-dth/2)+a8*pow(-dth/2,2)+2*a9*raux*pow(-dth/2,2))
//							/dth/st;
//
//						raux+=dr/st;
//					}
//					dr=caca*dr;
//					dth=caca*dth;
//					
//				}

//				Baux[i][j]+=(dBr+dBth)*sines[j];
//				Baux[i+1][j]=B[i+1][j]-dBr*sines[j];
//				Baux[i][j+1]-=dBth*sines[j+1];
//
//				Aaux[i][j]=dBr*sines[j];
//				Baux[i][j]=dBth*sines[j];


				Baux[i][j]+=(dBr)*sines[j];
				Baux[i+1][j]=B[i+1][j]-dBr*sines[j];

				Aaux[i][j]+=(dBth)*sines[j];
				Aaux[i][j+1]-=dBth*sines[j+1];
			}
		}
		//pass values from auxiliary array Baux
		for(int i=1;i<rNum-1;i++){
			for(int j=1;j<thNum-1;j++){
				//check for blowups, exit program if that happens
				if(isinf(Baux[i][j])){
					io::report_blowup(k,i,j);
					return 1;
				}
				B[i][j]=Baux[i][j];
			}
		}
		//pass values from auxiliary array Aaux
		for(int i=1;i<2*rNum-2;i++){
			for(int j=1;j<2*thNum-2;j++){
				//check for blowups, exit program if that happens
				if(isinf(Aaux[i][j])){
					io::report_blowup(k,i,j);
					return 1;
				}
				A[i][j]=Aaux[i][j];
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
	//initialize variables
#ifdef TOROIDAL
	double *integrals=new double[2];
	integrals[0]=integrals[1]=0;
#else
	double *integrals=new double[3];
	integrals[0]=integrals[1]=integrals[2]=0;
#endif
	double r,th;

	//solve toroidal quantities
	for(int i=1;i<rNum-1;i++){
		r=i*dr+rmin;
		for(int j=1;j<thNum-1;j++){
			//Toroidal flux
			integrals[0]+=B[i][j]/sines[j];
			//Toroidal energy
			integrals[1]+=pow(B[i][j],2)/r;
		}
	}
	integrals[0]=integrals[0]*dr*dth;
	integrals[1]=integrals[1]*dr*dth;

	//solve poloidal quantities
#ifndef TOROIDAL
	for(int i=1;i<2*rNum-2;i++){
		r=rmin+i*dr/2;
		for(int j=1;j<2*thNum-2;j++){
			th=j*dth/2;
			//Poloidal energy
			integrals[2]+=pow((A[i+1][j]-A[i-1][j])/4/dr,2)/sines[j]+pow((A[i][j+1]-A[i][j-1])/4/dth,2)/r/r/sin(th);
		}
	}
	integrals[2]=integrals[2]*dr*dth/4;
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
	//maximum l for the spherical harmonics expansion
	int l=5;
	//solve the required legendre polynomial on all required points
	double P[l+1][2*thNum-1];
	for(int i=0;i<=l;i++){
		for(int j=1;j<2*thNum-2;j++){
			P[i][j]=gsl_sf_legendre_Pl(i,cos(j*dth/2));
		}
	}
	//solve integrated coefficients
	double a[l];
	double f0=0,f1=0,f2=0;
	for(int i=0;i<l;i++){
		a[i]=0;
<<<<<<< HEAD
		for(int j=1;j<2*thNum-2;j++){
			a[i]+=(A[2*rNum-3][j+1]-A[2*rNum-3][j-1])*P[i+1][j]*dr/2*(i+1)*(2*i+3)/(i+2)/4*pow(1-dr/2,i+1);
=======
		for(int j=1;j<thNum-2;j++){
			double th=j*dth;
			f1=A[rNum-2][j]*(cos(j*dth)*P[i+1][j]-P[i][j])/sin(j*dth);
			f0=f2=0;
			if(j!=thNum-2){
				f2=A[rNum-2][j+1]*(cos(th+dth)*P[i+1][j+1]-P[i][j+1])/sin(th+dth);
			}
			if(j!=1){
				f0=A[rNum-2][j-1]*(cos(th-dth)*P[i+1][j-1]-P[i][j-1])/sin(th-dth);
			}
			a[i]+=(f0+4*f1+f2)*dth/3;
>>>>>>> Use hybrid approach in evolution of B
		}
		a[i]=a[i]*2*(4*atan(1))*pow(1-dr,i+1)/(i+2)*sqrt((2*i+3)/4.0/4.0/atan(1))*(i+1)/2;
	}
<<<<<<< HEAD
	//set boundary to A[rNum-2][j] to start summation
	for(int j=1;j<2*thNum-2;j++){
		A[2*rNum-2][j]=A[2*rNum-3][j];
	}
	//permorm sumation
	for(int i=0;i<l;i++){
		for(int j=1;j<2*thNum-2;j++){
			A[2*rNum-2][j]+=a[i]*(cos(j*dth/2)*P[i+1][j]-P[i][j])/pow(1-dr/4,i+2);
=======
	for(int j=1;j<thNum-1;j++){
		A[rNum-1][j]=0;
	}
	//permorm sumation
	for(int i=0;i<l;i++){
		for(int j=1;j<thNum-1;j++){
			A[rNum-1][j]+=a[i]*sqrt((2*i+3)/4.0/4.0/atan(1))*(cos(j*dth)*P[i+1][j]-P[i][j]);
>>>>>>> Use hybrid approach in evolution of B
		}
	}
	return;
}		/* -----  end of function solve_A_boundary  ----- */
#endif
#endif

}		/* -----  end of namespace sim  ----- */
