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
				Aaux[i][j]=A[i][j];
			}
		}
#endif

		//Update toroidal field function
		int st=10;
		double a1=0,a2=0,a3=0,a4=0,a5=0,a6=0;
		double b1=0,b2=0,b3=0,b4=0,b5=0,b6=0;
		double c1=0,c2=0,c3=0,c4=0,c5=0,c6=0;
		for(int i=0;i<rNum-1;i++){
			double r=rmin+i*dr;
			for(int j=0;j<thNum-1;j++){
			double th=j*dth;
				if(i==0&&j==0)
					continue;
				double dBr=0;
				double dBth=0;
				if(j!=0){
					double raux=r+dr/2;
					double thaux=th-dth/2+dth/2/st;
					a1=(B[i][j]+B[i+1][j])/2;
					a2=(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1])/4/dth;
					a3=(B[i+1][j]-B[i][j])/dr;
					a4=(B[i][j+1]+B[i+1][j+1]-2*B[i][j]-2*B[i+1][j]+B[i][j-1])+B[i+1][j-1]/4/dth/dth;
					a5=(B[i+1][j+1]+B[i][j-1]-B[i][j+1]-B[i+1][j-1])/2/dr/dth;
					a6=(-B[i][j+1]+B[i+1][j+1]+2*B[i][j]-2*B[i+1][j]-B[i][j-1]+B[i+1][j-1])/2/dr/dth/dth;

					b1=(A[i][j]+A[i+1][j])/2;
					b2=(A[i][j+1]+A[i+1][j+1]-A[i][j-1]-A[i+1][j-1])/4/dth;
					b3=(A[i+1][j]-A[i][j])/dr;
					b4=(A[i][j+1]+A[i+1][j+1]-2*A[i][j]-2*A[i+1][j]+A[i][j-1])+A[i+1][j-1]/4/dth/dth;
					b5=(A[i+1][j+1]+A[i][j-1]-A[i][j+1]-A[i+1][j-1])/2/dr/dth;
					b6=(-A[i][j+1]+A[i+1][j+1]+2*A[i][j]-2*A[i+1][j]-A[i][j-1]+A[i+1][j-1])/2/dr/dth/dth;

					c1=(gsA[i][j]+gsA[i+1][j])/2;
					c2=(gsA[i][j+1]+gsA[i+1][j+1]-gsA[i][j-1]-gsA[i+1][j-1])/4/dth;
					c3=(gsA[i+1][j]-gsA[i][j])/dr;
					c4=(gsA[i][j+1]+gsA[i+1][j+1]-2*gsA[i][j]-2*gsA[i+1][j]+gsA[i][j-1])+gsA[i+1][j-1]/4/dth/dth;
					c5=(gsA[i+1][j+1]+gsA[i][j-1]-gsA[i][j+1]-gsA[i+1][j-1])/2/dr/dth;
					c6=(-gsA[i][j+1]+gsA[i+1][j+1]+2*gsA[i][j]-2*gsA[i+1][j]-gsA[i][j-1]+gsA[i+1][j-1])/2/dr/dth/dth;
					for(int n=0;n<st;n++){
						dBr+=dt*thtd*initial::eta(raux,thaux)/sin(thaux)
							*(a3+a5*(thaux-th)+a6*pow(thaux-th,2))/st;
						dBr+=dt*initial::chi(raux,thaux)
							*(a1+a2*(thaux-th)+a3*(raux-r)+a4*pow(thaux-th,2)+a5*(raux-r)*(thaux-th)+a6*(raux-r)*pow(thaux-th,2))
							*(a2+2*a4*(thaux-th))/st;
//						dBr+=dt*initial::chi(raux,thaux)
//							*(c1+c2*(thaux-th)+c3*(raux-r)+c4*pow(thaux-th,2)+c5*(raux-r)*(thaux-th)+c6*(raux-r)*pow(thaux-th,2))
//							*(b2+2*b4*(thaux-th))/st;
						thaux+=dth/st;
					}
					dBr+=dt*initial::chi(r+dr/2,th)*(gsA[i+1][j]+gsA[i][j])/2*b2;				
//					if(i==rNum-2){
//						std::cout<<initial::chi(r+dr/2,th)*(gsA[i+1][j]+gsA[i][j])/2*b2<<" "<<i<<std::endl;
//					}
//					if(i==rNum-3){
//						std::cout<<initial::chi(r+dr/2,th)*(gsA[i+1][j]+gsA[i][j])/2*b2<<" "<<i<<std::endl;
//					}
					if(i==rNum-2){
						std::cout<<dBr<<" "<<i<<std::endl;
					}
					if(i==rNum-3){
						std::cout<<dBr<<" "<<i<<std::endl;
					}
					dBr=dBr/dr;
				}
				if(i!=0){
					double raux=r-dr/2+dr/2/st;
					double thaux=th+dth/2;
					a1=(B[i][j]+B[i][j+1])/2;
					a2=(B[i+1][j]+B[i+1][j+1]-B[i-1][j]-B[i-1][j+1])/4/dr;
					a3=(B[i][j+1]-B[i][j])/dth;
					a4=(B[i+1][j]+B[i+1][j+1]-2*B[i][j]-2*B[i][j+1]+B[i-1][j])+B[i-1][j+1]/4/dr/dr;
					a5=(B[i+1][j+1]+B[i-1][j]-B[i+1][j]-B[i-1][j+1])/2/dr/dth;
					a6=(-B[i+1][j]+B[i+1][j+1]+2*B[i][j]-2*B[i][j+1]-B[i-1][j]+B[i-1][j+1])/2/dr/dth/dth;

					b1=(A[i][j]+A[i][j+1])/2;
					b2=(A[i+1][j]+A[i+1][j+1]-A[i-1][j]-A[i-1][j+1])/4/dr;
					b3=(A[i][j+1]-A[i][j])/dth;
					b4=(A[i+1][j]+A[i+1][j+1]-2*A[i][j]-2*A[i][j+1]+A[i-1][j])+A[i-1][j+1]/4/dr/dr;
					b5=(A[i+1][j+1]+A[i-1][j]-A[i+1][j]-A[i-1][j+1])/2/dr/dth;
					b6=(-A[i+1][j]+A[i+1][j+1]+2*A[i][j]-2*A[i][j+1]-A[i-1][j]+A[i-1][j+1])/2/dr/dth/dth;

					c1=(gsA[i][j]+gsA[i][j+1])/2;
					c2=(gsA[i+1][j]+gsA[i+1][j+1]-gsA[i-1][j]-gsA[i-1][j+1])/4/dr;
					c3=(gsA[i][j+1]-gsA[i][j])/dth;
					c4=(gsA[i+1][j]+gsA[i+1][j+1]-2*gsA[i][j]-2*gsA[i][j+1]+gsA[i-1][j])+gsA[i-1][j+1]/4/dr/dr;
					c5=(gsA[i+1][j+1]+gsA[i-1][j]-gsA[i+1][j]-gsA[i-1][j+1])/2/dr/dth;
					c6=(-gsA[i+1][j]+gsA[i+1][j+1]+2*gsA[i][j]-2*gsA[i][j+1]-gsA[i-1][j]+gsA[i-1][j+1])/2/dr/dth/dth;

					for(int n=0;n<st;n++){
						dBth+=dt*thtd*initial::eta(raux,thaux)/sin(thaux)/raux/raux
							*(a3+a5*(raux-r)+a6*pow(raux-r,2))/st;
						dBth+=-dt*initial::chi(raux,thaux)
							*(a1+a2*(raux-r)+a3*(thaux-th)+a4*pow(raux-r,2)+a5*(thaux-th)*(raux-r)+a6*(thaux-th)*pow(raux-r,2))
							*(a2+2*a4*(raux-r))/st;
//						dBth+=-dt*initial::chi(raux,thaux)
//							*(c1+c2*(raux-r)+c3*(thaux-th)+c4*pow(raux-r,2)+c5*(thaux-th)*(raux-r)+c6*(thaux-th)*pow(raux-r,2))
//							*(b2+2*b4*(raux-r))/st;
						raux+=dr/st;
					}
					dBth=dBth/dth;
				}
				Baux[i][j]+=(dBr+dBth)*sines[j];
				Baux[i+1][j]=B[i+1][j]-dBr*sines[j];
				Baux[i][j+1]-=dBth*sines[j+1];
				Aaux[i][j]=dBr;
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
