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
#include	<fstream>
#include	<math.h>
//Number of steps in r and theta
#define rNum 50
#define thNum 20
//value of the minimal r inside which there is no MF
#define rmin 0.5
//Ratio of hall to dissipation timescales
#define thtd 0.01
//define timestep and number of timesteps in simulation
#define dt 0.00001
#define tNum 1000000
#define plotSteps 1000

#define conservative

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Bi
 *  Description:  Gives the initial conditions for B at a given r and theta
 * =====================================================================================
 */
	double
Bi ( double r, double th )
{
	return pow(r-rmin,3)*pow(1-r,3)*pow(sin(th),3);
}		/* -----  end of function Bi  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  chi
 *  Description:  gives the value of chi at certain radius
 * =====================================================================================
 */
	double
chiValue ( double r )
{
	return 1.0/(1.0-r*r);
}		/* -----  end of function chi  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  dchi
 *  Description:  returns the derivative of xi at a certain radius
 * =====================================================================================
 */
	double
dchiValue ( double r )
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
	double dr=(1.0-rmin)/rNum;
	double dth=Pi/thNum;

	//create array to store results in each timestep
	double results[rNum][tNum/plotSteps+1];
	double results2[rNum][tNum/plotSteps+1];
	double integrals[tNum/plotSteps+1];
	double integrals2[tNum/plotSteps+1];

	//Create arrays for B, chi, and dchi, and the needed sines and cosines
	double B[rNum][thNum];
	double Baux[rNum][thNum];
#ifdef conservative
	double chi[rNum];
#else
	double dchi[rNum];
	double cosines[thNum];
#endif
	double sines[thNum];

	//variables used to store the coordinate values at a given point
	double r,th;

	//set all arrays to their initial values
	for(int i=1;i<rNum-1;i++){
		r=i*dr+rmin;
#ifdef conservative
		//The values of chi are stored at midpoints, as needed by the
		//conservative scheme
		chi[i]=chiValue(r+dr/2);
#else
		dchi[i]=dchiValue(r);
#endif
		for(int j=1;j<thNum-1;j++){
			th=j*dth;
			B[i][j]=Bi(r,th);
		}
	}
	for(int j=1;j<thNum-1;j++){
		th=j*dth;
#ifndef conservative
		cosines[j]=cos(th);
#endif
		sines[j]=sin(th);
	}
	//solve all the commonly used constant terms in the finite difference scheme
#ifndef conservative
	double A1[rNum][thNum];
	double A2[rNum][thNum];
	double A3[rNum][thNum];
	double A4[rNum][thNum];
	double A5[rNum][thNum];
	for(int i=1;i<rNum-1;i++){
		r=i*dr+rmin;
		for(int j=1;j<thNum-1;j++){
			A1[i][j]=sines[j]/dth/2*dchi[i][j];
			A2[i][j]=r*pow(sines[j],2)/dr*thtd;
			A3[i][j]=pow(r*sines[j],2)/dr*thtd;
			A4[i][j]=cosines[j]*sines[j]/dth/2*thtd;
			A5[i][j]=pow(sines[j],2)/dth*thtd;
		}
	}
#endif

	//set boundary values for B
	for(int i=0;i<rNum;i++)
		B[i][0]=B[i][thNum-1]=0;
	for(int j=0;j<thNum;j++)
		B[rNum-1][j]=0;

	//simulate, use dBdt to store temporal derivative to simplify code
	double dBdt;
	int plotStep=0;
	for(int k=0;k<tNum;k++){
		for(int i=1;i<rNum-1;i++){
			r=i*dr+rmin;
			for(int j=1;j<thNum-1;j++){
#ifdef conservative
				dBdt=0;
				//add hall contribution
				dBdt+=sines[j]/dr*(
						(B[i+1][j]+B[i][j])/2.0
						*(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1])/4.0/dth
						*chi[i]
						);
				dBdt+=-sines[j]/dr*(
						(B[i][j]+B[i-1][j])/2.0
						*(B[i][j+1]+B[i-1][j+1]-B[i][j-1]-B[i-1][j-1])/4.0/dth
						*chi[i-1]
						);
				dBdt+=sines[j]/dth*(
						(B[i][j+1]+B[i][j])/2.0
						*(B[i+1][j+1]+B[i+1][j]-B[i-1][j]-B[i-1][j+1])/4.0/dth
						*chi[i]
						);
				dBdt+=-sines[j]/dth*(
						(B[i][j]+B[i][j-1])/2.0
						*(B[i+1][j-1]+B[i+1][j]-B[i-1][j]-B[i-1][j-1])/4.0/dth
						*chi[i]
						);
				//add resistive contribution
				dBdt+=thtd/dr*(
						(B[i+1][j]-B[i][j])/dr
						);
				dBdt+=-thtd/dr*(
						(B[i][j]-B[i-1][j])/dr
						);
				dBdt+=sines[j]*thtd/dth/r/r*(
						1/sin(dth*(2*j+1))
						*(B[i][j+1]-B[i][j])/dth
						);
				dBdt+=-sines[j]*thtd/dth/r/r*(
						1/sin(dth*(2*j-1))
						*(B[i][j]-B[i][j-1])/dth
						);
#else
				//add hall contribution
//				dBdt=B[i][j]*sines[j]/dth/2*(B[i][j+1]-B[i][j-1])*dchi[i][j];
				dBdt=B[i][j]*(B[i][j+1]-B[i][j-1])*A1[i][j];
				//add resistive contribution
//				dBdt+=thtd*(
//						r*pow(sines[j],2)/dr*(B[i+1][j]-B[i-i][j])
//						+pow(r*sines[j],2)/dr*(B[i+1][j]+B[i-1][j]-2*B[i][j])
//						+cosines[j]*sines[j]/dth/2*(B[i][j+1]-B[i][j-1])
//						+pow(sines[j],2)/dth*(B[i][j+1]+B[i][j-1]-2*B[i][j])
//						);
				dBdt+=	A2[i][j]*(B[i+1][j]-B[i-i][j])
						+A3[i][j]*(B[i+1][j]+B[i-1][j]-2*B[i][j])
						+A4[i][j]*(B[i][j+1]-B[i][j-1])
						+A5[i][j]*(B[i][j+1]+B[i][j-1]-2*B[i][j]);
#endif
				//update value in auxiliary variable
				Baux[i][j]=B[i][j]+dt*dBdt;
			}
		}
		//update array with auxiliary values, and solve supposedly conserved quantities
		double integral=0;
		double integral2=0;
		for(int i=1;i<rNum-1;i++){
			for(int j=1;j<thNum-1;j++){
				if(isinf(Baux[i][j])){
					std::cout << "Blew up at step " << k << "in place " << i << "," << j << std::endl;
					return 1;
				}
				B[i][j]=Baux[i][j];
				integral+=B[i][j]*sines[j];
				integral2+=B[i][j]*i*dr;
			}
		}
		//store results to be logged
		if(k%plotSteps==0){
			std::cout << k << "/" << tNum << std::endl;
			for(int i=0;i<rNum;i++){
				results[i][plotStep]=B[i][thNum/2];
			}
			for(int j=0;j<thNum;j++){
				results2[j][plotStep]=B[rNum/2][j];
			}
			integrals[plotStep]=integral;
			integrals2[plotStep]=integral2;
			plotStep++;
		}
	}

	//log results to file
	//resultsR stores the results along a line with theta=Pi/2
	std::ofstream myfile;
	myfile.open("resultsR.dat");
	for(int i=0;i<rNum;i++){
		r=i*dr+rmin;
		myfile << r << " ";
		for(int k=0;k<tNum/plotSteps;k++){
			myfile << results[i][k] << " ";
		}
		myfile << std::endl;
	}
	myfile.close();

	//resultsTH stores the results along a line with r=R/2
	std::ofstream myfile2;
	myfile2.open("resultsTH.dat");
	for(int j=0;j<thNum;j++){
		th=j*dth;
		myfile2 << th << " ";
		for(int k=0;k<tNum/plotSteps;k++){
			myfile2 << results2[j][k] << " ";
		}
		myfile2 << std::endl;
	}
	myfile2.close();

	//resultsINT stores the values at each timestep of the supposedly constant integrals
	std::ofstream myfile3;
	myfile3.open("resultsINT.dat");
	double t;
	for(int k=0;k<tNum/plotSteps;k++){
		t=k*dt;
		myfile3 << t << " ";
		myfile3 << integrals[k] << " ";
		myfile3 << integrals2[k] << " ";
		myfile3 << std::endl;
	}
	myfile3.close();

	return 0;
}				/* ----------  end of function main  ---------- */
