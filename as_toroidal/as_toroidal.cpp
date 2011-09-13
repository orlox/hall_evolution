/*
 * =====================================================================================
 *
 *       Filename:  as_toroidal.cpp
 *
 *    Description:  This code follows the evolution of an axisymmetric purely toroidal
 *    MF on a NS due to the hall effect.
 *
 *        Version:  1.0
 *        Created:  05/09/11 19:18:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Pablo Marchant
 *        Company:  Pontificia Universidad Católica de Chile
 *
 * =====================================================================================
 */

#include	<iostream>
#include	<math.h>
//Number of steps in r and theta
#define rNum 100
#define thNum 100
//Ratio of hall to dissipation timescales
#define thtd 0.01
//define timestep and number of timesteps in simulation
#define dt 0.001
#define tNum 1000000
#define plotSteps 1000

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Bi
 *  Description:  Gives the initial conditions for B at a given r and theta
 * =====================================================================================
 */
	double
Bi ( float r, float th )
{
	return pow(r,2)*pow(1-r,2)*pow(sin(th),2)*16;
}		/* -----  end of function Bi  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  chi
 *  Description:  gives the value of chi at certain radius
 * =====================================================================================
 */
	float
chiValue ( float r )
{
	return 1.0/(1.0-r*r);
}		/* -----  end of function chi  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  dchi
 *  Description:  returns the derivative of xi at a certain radius
 * =====================================================================================
 */
	float
dchiValue ( float r )
{
	return 2*r/pow(1-r*r,2);
}		/* -----  end of function dchi  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  Sets the initial conditions, runs the simulation loop, and logs
 *  the results.
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	double const Pi=4*atan(1);
	//define size of steps
	float dr=1.0/rNum;
	float dth=Pi/thNum;

	//Create arrays for B, chi, and dchi, and the needed sines and cosines
	float B[rNum][thNum];
	float Baux[rNum][thNum];
	//float chi[rNum][thNum];
	float dchi[rNum][thNum];
	float cosines[thNum];
	float sines[thNum];

	//variables used to store the coordinate values at a given point
	float r,th;

	//set all arrays to their initial values
	for(int i=1;i<rNum-1;i++){
		r=i*dr;
		for(int j=1;j<thNum-1;j++){
			th=j*dth;
			B[i][j]=Bi(r,th);
			//chi[i][j]=chiValue(r);
			dchi[i][j]=dchiValue(r);
		}
	}
	for(int j=1;j<thNum-1;j++){
		th=j*dth;
		cosines[j]=cos(th);
		sines[j]=sin(th);
	}

	//set boundary values for B
	for(int i=0;i<rNum;i++)
		B[i][0]=B[i][thNum-1]=0;
	for(int j=0;j<thNum;j++)
		B[rNum-1][j]=0;

	//simulate, use dBdt to store temporal derivative to simplify code
	float dBdt;
	for(int k=0;k<tNum;k++){
		for(int i=1;i<rNum-1;i++){
			r=i*dr;
			for(int j=1;j<thNum-1;j++){
				//add hall contribution
				dBdt=B[i][j]*sines[j]/dth/2*(B[i][j+1]-B[i][j-1])*dchi[i][j];
				//add resistive contribution
				dBdt+=thtd*(
						r*pow(sines[j],2)/dr*(B[i+1][j]-B[i-i][j])
						+pow(r*sines[j],2)/dr*(B[i+1][j]+B[i-1][j]-2*B[i][j])
						+cosines[j]*sines[j]/dth/2*(B[i][j+1]-B[i][j-1])
						+pow(sines[j],2)/dth*(B[i][j+1]+B[i][j-1]-2*B[i][j])
						);
				//update value in auxiliary variable
				Baux[i][j]=B[i][j]+dt*dBdt;
			}
		}
		//update array with auxiliary values

		for(int i=1;i<rNum-1;i++){
			for(int j=1;j<thNum-1;j++){
				if(isinf(Baux[i][j])){
					std::cout << "Blew up at step " << k << "in place " << i << "," << j << std::endl;
					return 1;
				}
				B[i][j]=Baux[i][j];
			}
		}
		if(k%plotSteps==0)
			std::cout << k << "\t" << B[82][19] << "\t" << B[83][19]<< "\t" << B[84][19]<< std::endl;
	}

	return 0;
}				/* ----------  end of function main  ---------- */
