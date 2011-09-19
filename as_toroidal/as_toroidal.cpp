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
#include	<sstream>
#include	<fstream>
#include	<math.h>
#include	<time.h>
#include	<string>
#include	<stdlib.h>

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
	//value of the function used to define the region where the toroidal field is present
	double alpha=(35.0/8.0*pow(r,2)-21.0/4.0*pow(r,4)+15.0/8.0*pow(r,6))*pow(sin(th),2);
	if(alpha<=1)
		return 0;
	else
		return 40.0*pow(alpha-1,2);
}		/* -----  end of function Bi  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  chi
 *  Description:  gives the value of chi at certain point in the star
 * =====================================================================================
 */
	double
chiValue ( double r, double th )
{
	return 0.25/(1.0-r*r)/pow(r*sin(th),2);
}		/* -----  end of function chi  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  dchi
 *  Description:  returns the derivative of xi at a certain point in the star
 * =====================================================================================
 */
	double
dchiValue ( double r, double th )
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
	//Number of steps in r and theta
	int rNum=250;
	int thNum=90;
	//value of the minimal r inside which there is no MF
	double rmin=0.0;
	//size of the theta region, used for practical reasons to define well the zone
	//that contains the MF without unnecesary points
	double thsize=Pi;
	//Ratio of hall to dissipation timescales
	double thtd=0.01;
	//define timestep and number of timesteps in simulation
	double dt=0.0000001;
	int tNum=10000000;
	int plotSteps=100;

	//define size of steps
	double dr=(1.0-rmin)/rNum;
	double dth=thsize/thNum;

	//Create arrays for B, chi, and dchi, and the needed sines and cosines
	double B[rNum][thNum];
	double Baux[rNum][thNum];
#ifdef conservative
	//The values of chi are stored at midpoints, as needed by the conservative scheme
	//chiR stores the values at midpoints in the radial coordinate, while
	//chiTH stores the valus at midpoints in the polar coordinate
	double chiR[rNum][thNum];
	double chiTH[rNum][thNum];
#else
	double dchi[rNum][thNum];
	double cosines[thNum];
#endif
	double sines[thNum];

	//variables used to store the coordinate values at a given point, and the time
	double r,th,t;

	//set all arrays to their initial values
	for(int i=1;i<rNum-1;i++){
		r=i*dr+rmin;
		for(int j=1;j<thNum-1;j++){
			th=(Pi-thsize)/2+j*dth;
			B[i][j]=Bi(r,th);
#ifdef conservative
			chiR[i][j]=chiValue(r+dr/2,th);
			chiTH[i][j]=chiValue(r,th+dth/2);
#else
			dchi[i][j]=dchiValue(r,th);
#endif
		}
	}
	for(int j=1;j<thNum-1;j++){
		th=(Pi-thsize)/2+j*dth;
#ifndef conservative
		cosines[j]=cos(th);
#endif
		sines[j]=sin(th);
	}
	//set boundary values for B
	for(int i=0;i<rNum;i++)
		B[i][0]=B[i][thNum-1]=0;
	for(int j=0;j<thNum;j++)
		B[rNum-1][j]=0;

	//create directory for log files
	std::stringstream timeStream;
	timeStream << time(0);
	std::string mkdir="mkdir results_"+timeStream.str();
	system(mkdir.c_str());
	//resultsINT stores the values at each timestep of the supposedly constant integrals
	std::ofstream resultsINT;
	std::string filename="results_"+timeStream.str()+"/resultsINT.dat";
	resultsINT.open(filename.c_str());
	
	//simulate, use dBdt to store temporal derivative to simplify code
	double dBdt;
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
						*chiR[i][j]
						);
				dBdt+=-sines[j]/dr*(
						(B[i][j]+B[i-1][j])/2.0
						*(B[i][j+1]+B[i-1][j+1]-B[i][j-1]-B[i-1][j-1])/4.0/dth
						*chiR[i-1][j]
						);
				dBdt+=sines[j]/dth*(
						(B[i][j+1]+B[i][j])/2.0
						*(B[i+1][j+1]+B[i+1][j]-B[i-1][j]-B[i-1][j+1])/4.0/dth
						*chiTH[i][j]
						);
				dBdt+=-sines[j]/dth*(
						(B[i][j]+B[i][j-1])/2.0
						*(B[i+1][j-1]+B[i+1][j]-B[i-1][j]-B[i-1][j-1])/4.0/dth
						*chiTH[i][j-1]
						);
				//add resistive contribution
				dBdt+=thtd/dr*(
						(B[i+1][j]-B[i][j])/dr
						);
				dBdt+=-thtd/dr*(
						(B[i][j]-B[i-1][j])/dr
						);
				dBdt+=sines[j]*thtd/dth/r/r*(
						1/sin((Pi-thsize)/2+dth*(2*j+1.0)/2.0)
						*(B[i][j+1]-B[i][j])/dth
						);
				dBdt+=-sines[j]*thtd/dth/r/r*(
						1/sin((Pi-thsize)/2+dth*(2*j-1.0)/2.0)
						*(B[i][j]-B[i][j-1])/dth
						);
#else
				//add hall contribution
				dBdt=dchi[i]*B[i][j]*(B[i][j+1]-B[i][j-1])/dth/2.0;
				//add resistive contribution
				dBdt+=thtd*(
						(B[i+1][j]+B[i-1][j]-B[i][j])/pow(dr,2)
						+cosines[j]/sines[j]/pow(r,2)*(B[i][j+1]-B[i][j-1])/dth/2.0
						+1/pow(r,2)*(B[i][j+1]+B[i][j-1]-B[i][j])/pow(dth,2)
						);
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
				integral+=B[i][j]/sines[j];
				integral2+=B[i][j]*i*dr;
			}
		}
		//log results for the integrals
		if(k%plotSteps==0){
			std::cout << k << "/" << tNum << std::endl;
			t=k*dt;
			resultsINT << t << " " << integral << " " << integral2 << std::endl;
			//log results for this timestep
			std::stringstream fnStream;
			fnStream << "results_" << timeStream.str() << "/data_" << k;
			filename=fnStream.str();
			std::ofstream results;
			results.open(filename.c_str());
			for(int i=0;i<rNum;i++){
				r=i*dr+rmin;
				for(int j=0;j<thNum;j++){
					th=(Pi-thsize)/2+j*dth;
					results << r << " " << th << " " << B[i][j] << std::endl;
				}
				results << std::endl;
			}
		}
	}
	//close integrals log file
	resultsINT.close();

	return 0;
}				/* ----------  end of function main  ---------- */
