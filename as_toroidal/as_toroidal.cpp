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

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  chi
 *  Description:  gives the value of chi at certain point in the star
 * =====================================================================================
 */
	double
chiValue ( double r, double th )
{
	return 1/pow(r*sin(th),2);//0.25/(1.0-r*r)/pow(r*sin(th),2);
}		/* -----  end of function chi  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Bi
 *  Description:  Gives the initial conditions for B at a given r and theta, the minimun radius is part of the initial condition, so it's included in this block.
 * =====================================================================================
 */
//double rmin=0.5;
double rmin=0;

	double
Bi ( double r, double th )
{
	//This initial condition represent the fundamental eigenmode of pure ohm
	//evolution. k and b depend on 
	//double k=12.72135634741802;
	//double b=-4.205800731231387;
	//return ((sin(k*r)/(k*r)-cos(k*r))/(k*r)+(b*(-sin(k*r)-cos(k*r)/(k*r)))/(k*r))*sin(th)*r*sin(th);
	double k=4.493328387732195;
	return ((sin(k*r)/(k*r)-cos(k*r))/(k*r))*sin(th)*r*sin(th);
}		/* -----  end of function Bi  ----- */

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
	//check if number of arguments is correct
	if(argc!=7){
		std::cout<<"ERROR: Not enough arguments provided"<<std::endl;
		std::cout<<"usage is:"<<std::endl<<std::endl;
		std::cout<<"as_toroidal <rNum> <thNum> <dt> <tNum> <tNumPlot> <thtd>"<<std::endl<<std::endl;
		return 0;
	}
	//########################READ ARGUMENTS FROM ARGV######################
	//NO CHECK ARE DONE ON THE INPUT FROM TERMINAL!!
	//Number of steps in r and theta
	int rNum=atoi(argv[1]);
	int thNum=atoi(argv[2]);
	//define timestep and number of timesteps in simulation
	double dt=atof(argv[3]);
	int tNum=atoi(argv[4]);
	int plotSteps=atoi(argv[5]);
	//Ratio of hall to dissipation timescales
	double thtd=atof(argv[6]);
	//######################################################################

	//define size of steps
	double dr=(1.0-rmin)/(rNum-1);
	double const Pi=4*atan(1);
	double dth=Pi/(thNum-1);

	//Create arrays for B, chi, and dchi, and the needed sines and cosines
	double B[rNum][thNum];
	double Baux[rNum][thNum];
	double sines[thNum];

	//variables used to store the coordinate values at a given point, and the time
	double r,th,t;

	//set initial values
	for(int i=1;i<rNum-1;i++){
		r=i*dr+rmin;
		for(int j=1;j<thNum-1;j++){
			th=j*dth;
			B[i][j]=Bi(r,th);
		}
	}
	//set boundary values for B
	for(int i=0;i<rNum;i++)
		B[i][0]=B[i][thNum-1]=0;
	for(int j=0;j<thNum;j++)
		B[0][j]=B[rNum-1][j]=0;
	//solve sines at each point in the grid
	for(int j=1;j<thNum-1;j++){
		th=j*dth;
		sines[j]=sin(th);
	}
	//Solve repeated values in the calculations, includes multiplication by dt
	float hall_rflux_plus[rNum][thNum];
	float hall_thflux_plus[rNum][thNum];
	float res_rflux=dt*thtd/dr/dr;
	float res_thflux_plus[rNum][thNum];
	for(int i=0;i<rNum-1;i++){
		r=i*dr+rmin;
		for(int j=0;j<thNum-1;j++){
			th=j*dth;
			hall_rflux_plus[i][j]=dt*chiValue(r+dr/2,th)/8.0/dr/dth;
			hall_thflux_plus[i][j]=dt*chiValue(r,th+dth/2)/8.0/dr/dth;
			res_thflux_plus[i][j]=dt*thtd/dth/dth/r/r/sin(dth*(2*j+1.0)/2.0);
		}
	}

	//create directory for log files
	std::stringstream timeStream;
	timeStream << time(0);
	std::string mkdir="mkdir results_"+timeStream.str();
	system(mkdir.c_str());
	//resultsINT stores the values at each timestep of the supposedly constant integrals
	std::ofstream resultsINT;
	std::string filename="results_"+timeStream.str()+"/resultsINT.dat";
	resultsINT.open(filename.c_str());

	//Display info
	std::cout<<std::endl;
	std::cout<< "#################################################"<< std::endl;
	std::cout<< "##########Axisymmetric Hall Evolution############"<< std::endl;
	std::cout<< "#################################################"<< std::endl;
	std::cout<< "PARAMETERS CHOSEN:"<< std::endl;
	std::cout<< "-Number of Radial steps: "<< rNum << std::endl;
	std::cout<< "-Number of Angular steps: "<< thNum << std::endl;
	std::cout<< "-Size of timestep: "<< dt << std::endl;
	std::cout<< "-Number of Time steps: "<< tNum << std::endl;
	std::cout<< "-Save output every "<< plotSteps << " steps" << std::endl;
	std::cout<< "-Minimun radius: "<< rmin << std::endl;
	std::cout<< "-Ratio of hall to dissipative timescales: "<< thtd << std::endl;
	std::cout<<std::endl;
	std::cout<< "Results stored in folder results_" << timeStream.str() << std::endl;
	std::cout<<std::endl;

	//_params.dat stores the simulation parameters
	std::ofstream params;
	std::string filename2="results_"+timeStream.str()+"/params.dat";
	params.open(filename2.c_str());
	params<< "rNum:"<< rNum << std::endl;
	params<< "thNum:"<< thNum << std::endl;
	params<< "dt:"<< dt << std::endl;
	params<< "tNum:"<< tNum << std::endl;
	params<< "plotSteps:"<< plotSteps << std::endl;
	params<< "rmin:"<< rmin << std::endl;
	params<< "thtd:"<< thtd << std::endl;
	params.close();
	
	//simulate, use dBr and dBth to store changes due to flux
	//in the r and theta directions
	//flux and energy are the supposedly conserved quantities
	double dBr;
	double dBth;
	double flux=0;
	double energy=0;
	std::cout<< "Starting Simulation!" << std::endl;
	for(int k=0;k<=tNum;k++){
		//log results for the integrals and the conserved values, if it corresponds
		if(k%plotSteps==0){
			flux=0;
			energy=0;
			//Solve conserved quantities
			for(int i=1;i<rNum-1;i++){
				r=i*dr+rmin;
				for(int j=1;j<thNum-1;j++){
					flux+=B[i][j]/sines[j];
					energy+=pow(B[i][j],2)/r;
				}
			}
			flux=flux*dr*dth;
			energy=energy*dr*dth;
			//log conserved quantities to file
			t=k*dt;
			resultsINT << t << " " << flux << " " << energy << std::endl;
			//log the values of B for this timestep
			std::stringstream fnStream;
			fnStream << "results_" << timeStream.str() << "/data_" << k;
			filename=fnStream.str();
			std::ofstream results;
			results.open(filename.c_str());
			for(int i=0;i<rNum;i++){
				r=i*dr+rmin;
				for(int j=0;j<thNum;j++){
					results << B[i][j] << " ";
				}
				results << std::endl;
			}
			//put messege to std::cout to inform progress
			std::cout << k << "/" << tNum << std::endl;
		}
		//exit program if timesteps are done
		if(k==tNum){
			break;
		}
		for(int i=0;i<rNum-1;i++){
			for(int j=0;j<thNum-1;j++){
				dBr=0;
				dBth=0;
				//Hall contribution, radial flux
				//dBr+=dt*sines[j]/dr*(
				//		(B[i+1][j]+B[i][j])/2.0
				//		*(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1])/4.0/dth
				//		*chiR[i][j]
				//		);
				dBr+=hall_rflux_plus[i][j]
						*(B[i+1][j]+B[i][j])
						*(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1]);
				//Hall contribution, theta flux
				//dBth+=sines[j]/dth*(
				//		-(B[i][j+1]+B[i][j])/2.0
				//		*(B[i+1][j+1]+B[i+1][j]-B[i-1][j]-B[i-1][j+1])/4.0/dr
				//		*chiTH[i][j]
				//		);
				dBth+=-hall_thflux_plus[i][j]
						*(B[i][j+1]+B[i][j])
						*(B[i+1][j]+B[i+1][j+1]-B[i-1][j]-B[i-1][j+1]);
				//Resistive contribution, radial flux
				//dB+=thtd/dr*(
				//		(B[i+1][j]-B[i][j])/dr
				//		);
				dBr+=res_rflux*(B[i+1][j]-B[i][j])/sines[j];
				//Resistive contribution, theta flux
				//dB+=sines[j]*thtd/dth/r/r*(
				//		1/sin(dth*(2*j+1.0)/2.0)
				//		*(B[i][j+1]-B[i][j])/dth
				//		);
				dBth+=res_thflux_plus[i][j]*(B[i][j+1]-B[i][j]);
				//Update value in auxiliary array, updating also points forward
				//in the grid. Values in the boundaries are irrelevant, and not
				//passed from the auxiliary to the definite array.
				Baux[i][j]+=(dBr+dBth)*sines[j];
				Baux[i+1][j]=B[i+1][j]-dBr*sines[j];
				Baux[i][j+1]-=dBth*sines[j+1];
			}
		}
		//update array with auxiliary values
		for(int i=1;i<rNum-1;i++){
			for(int j=1;j<thNum-1;j++){
				//check for blowups, exit program if that happens
				if(isinf(Baux[i][j])){
					std::cout << "Blew up at step " << k << "in place " << i << "," << j << std::endl;
					return 1;
				}
				B[i][j]=Baux[i][j];
			}
		}
	}
	//close integrals log file
	resultsINT.close();

	return 0;
}				/* ----------  end of function main  ---------- */
