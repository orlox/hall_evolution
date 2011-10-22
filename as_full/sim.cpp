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
#include "math.h"
#include "sim.h"
#include "io.h"
#include "initial.h"
namespace sim{
//values given as cli arguments.
int rNum;
int thNum;
double dt;
int tNum;
int plotSteps;
int thtd;
//physical values of functions defining the magnetic field.
double **B;
double **A;
//physical values that describe the structure of the star.
double **chi;
double **eta;
//arrays that contains the precalculated quantities neccesary to solve fluxes for B
double **hall_rflux;
double **hall_thflux;
double **res_rflux;
double **res_thflux;
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
	A=new double*[rNum];
	B=new double*[rNum];
	chi=new double*[rNum];
	eta=new double*[rNum];
	hall_rflux=new double*[rNum];
	hall_thflux=new double*[rNum];
	res_rflux=new double*[rNum];
	res_thflux=new double*[rNum];
	for(int i=0;i<rNum;i++){
		A[i]=new double[thNum];
		B[i]=new double[thNum];
		chi[i]=new double[thNum];
		eta[i]=new double[thNum];
		hall_rflux[i]=new double[thNum];
		hall_thflux[i]=new double[thNum];
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
		for(int j=1;j<rNum-1;j++){
			th=j*dth;
			A[i][j]=initial::A(r,th);
			B[i][j]=initial::B(r,th);
			chi[i][j]=initial::chi(r,th);
			eta[i][j]=initial::eta(r,th);
		}
	}

	//set boundary conditions
	for(int i=0;i<rNum;i++){
		A[i][0]=A[i][thNum-1]=0;
		B[i][0]=A[i][thNum-1]=0;
	}
	for(int j=0;j<thNum;j++){
		A[0][j]=A[rNum-1][thNum-1]=0;
		B[0][j]=A[rNum-1][thNum-1]=0;
	}
	
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
	for(int i=1;i<rNum-1;i++){
		r=rmin+i*dr;
		for(int j=0;j<rNum;j++){
			th=j*dth;
			hall_rflux[i][j] = dt*initial::chi(r+dr/2,th)/8.0/dr/dth;
			hall_thflux[i][j]=-dt*initial::chi(r,th+dth/2)/8.0/dr/dth;
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
	return 0;
}		/* -----  end of function simulate  ----- */

}		/* -----  end of namespace sim  ----- */
