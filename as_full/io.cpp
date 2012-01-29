/*
 * =====================================================================================
 *
 *       Filename:  io.cpp
 *
 *    Description:  This file has all functions pertaining input and output. This includes reading cli parameters, and producing all the required output files.
 *
 *        Version:  1.0
 *        Created:  10/21/2011 08:54:49 PM
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
#include	<string>
#include	<time.h>
#include	<stdlib.h>
#include	<signal.h>
#include	"string.h"
#include	"io.h"
#include	"sim.h"

using namespace std;
namespace io{
//std::stringstream used to store the timestamp
stringstream timeStream;
//std::ofstream which points to the file where integrated quantities are stored
ofstream integrals_file;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  read_args
 *  Description:  Reads the cli arguments, checks them for integrity, and stores them for the simulation.
 *
 *  Exit code 0 implies success, 1 implies failure
 * =====================================================================================
 */
	int
read_args ( int argc, char *argv[] )
{
	//Check if correct number of arguments is provided
#ifdef TOROIDAL
	if(argc!=8)
#else
	if(argc!=9)
#endif
	{
		cerr<<"ERROR: Not enough arguments provided"<<endl;
		cerr<<"usage is:"<<endl;
		cerr<<"as_full_run <rNum(int)> <thNum(int)> <dt(float)> <tNum(int)> <tNumPlot(int)> <th/td(float)> <rless(int)>";
#ifndef TOROIDAL
		cerr<<" <l(int)>";
#endif
		cerr<<endl;
		return 1;
	}
	//Check integrity of all parameters
	if(
			!io::is_integer(argv[1])||
			!io::is_integer(argv[2])||
			!io::is_float(argv[3])||
			!io::is_integer(argv[4])||
			!io::is_integer(argv[5])||
			!io::is_float(argv[6])||
			!io::is_integer(argv[7])
#ifndef TOROIDAL
			||!io::is_integer(argv[8])
#endif
			)
	{
		cerr<<"ERROR: Arguments not in the right format"<<endl;
		cerr<<"usage is:"<<endl;
		cerr<<"as_full_run <rNum(int)> <thNum(int)> <dt(float)> <tNum(int)> <tNumPlot(int)> <th/td(float)> <rless(int)>";
#ifndef TOROIDAL
		cerr<<" <l(int)>";
#endif
		return 1;
	}
	//Store values for simulation
	sim::rNum=atoi(argv[1]);
	sim::thNum=atoi(argv[2]);
	sim::dt=atof(argv[3]);
	sim::tNum=atoi(argv[4]);
	sim::plotSteps=atoi(argv[5]);
	sim::thtd=atof(argv[6]);
	sim::rless=atoi(argv[7]);
#ifndef TOROIDAL
	sim::l=atoi(argv[8]);
#endif
	return 0;
}		/* -----  end of function read_args  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  is_integer
 *  Description:  Check if given string is an int. Returns 1 if it is, else, returns 0.
 * =====================================================================================
 */
	int
is_integer ( char *value )
{
	//Check that each individual character is an integer and return 0 if one is not
	for (unsigned int i=0;i<strlen(value);i++)
	{
		if (!isdigit(value[i]))
			return 0;
	}
	return 1;
}		/* -----  end of function is_integer  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  is_float
 *  Description:  Check if given string is a float. Returns 1 if it is, else, returns 0.
 * =====================================================================================
 */
	int
is_float ( char *value )
{
	//Check that each individual character is an integer and return 0 if one is not.
	//Also, it accepts a single point.
	int points=0;
	for (unsigned int i=0;i<strlen(value);i++)
	{
		if(points==0&&value[i]=='.'){
			points=1;
			continue;
		}
		if (!isdigit(value[i]))
			return 0;
	}
	return 1;
}		/* -----  end of function is_float  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_folder
 *  Description:  Creates folder using as a distinct id the current timestamp. It also creates the file params.dat which stores the simulation parameters,
 *  and appends a line with info to the file results_summary
 * =====================================================================================
 */
	void
create_folder ( )
{
	//Store timestamp, will be needed also later to access folder
	timeStream << time(0);
	//Create directory. THIS WONT WORK ON WINDOWS!!
	string mkdir="mkdir results_"+timeStream.str();
	system(mkdir.c_str());
	//Create params.dat file, where simulation parameters are logged
	ofstream params;
	string filename="results_"+timeStream.str()+"/params.dat";
	params.open(filename.c_str());
	params << "rNum:"<< sim::rNum << std::endl;
	params << "thNum:"<< sim::thNum << std::endl;
	params << "dt:"<< sim::dt << std::endl;
	params << "tNum:"<< sim::tNum << std::endl;
	params << "plotSteps:"<< sim::plotSteps << std::endl;
	params << "rmin:"<< sim::rmin << std::endl;
	params << "thtd:"<< sim::thtd << std::endl;
	params << "rless:"<< sim::rless << std::endl;
#ifndef TOROIDAL
	params << "l:"<< sim::l << std::endl;
#endif
	params.close();

	//Append line with summary to results_summary
	ofstream summary ("results_summary",ios::app);
	summary << "results_" << timeStream.str() << ": ";
	summary << sim::rNum << "x" << sim::thNum << "x" << sim::dt << ", ";
	summary << "thtd:" << sim::thtd << ", ";
	summary << "rmin:" << sim::rmin << ", ";
	summary << "rless:" << sim::rless << ", ";
#ifndef TOROIDAL
	summary << "l:" << sim::l << ", ";
#endif
	//Specify compiler build options, T=TOROIDAL, O=PUREOHM
	summary << "COMP_OPT:";
#ifdef TOROIDAL
	summary << "T";
#endif
#ifdef PUREOHM
	summary << "O";
#endif
	summary << "." << endl;
	summary.close();

	return;
}		/* -----  end of function create_folder  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  print_header
 *  Description:  Prints block of text describing the simulation.
 * =====================================================================================
 */
	void
print_header ( )
{
	cout << endl;
	cout << "#################################################"<< endl;
	cout << "##########Axisymmetric Hall Evolution############"<< endl;
	cout << "#################################################"<< endl;
	cout << "PARAMETERS CHOSEN:"<< endl;
	cout << "-Number of Radial steps: "<< sim::rNum << endl;
	cout << "-Number of Angular steps: "<< sim::thNum << endl;
	cout << "-Size of timestep: "<< sim::dt << endl;
	cout << "-Number of Time steps: "<< sim::tNum << endl;
	cout << "-Save output every "<< sim::plotSteps << " steps" << endl;
	cout << "-Minimun radius: "<< sim::rmin << endl;
	cout << "-Ratio of hall to dissipative timescales: "<< sim::thtd << endl;
	cout << "-Radial points excluded in beta calculations: "<< sim::rless << endl;
#ifndef TOROIDAL
	cout << "-Multipoles used for external field: "<< sim::l << endl;
#endif
	cout << endl;
	cout << "Results stored in folder results_" << timeStream.str() << endl;
	cout << endl;
	cout << "Beggining simulation" << endl;
	return;
}		/* -----  end of function print_header  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  create_integrals_file
 *  Description:  Crease file integrals.dat where integrated quantities are logged
 * =====================================================================================
 */
	void
create_integrals_file ( )
{
	string filename="results_"+timeStream.str()+"/integrals.dat";
	integrals_file.open(filename.c_str());
#ifdef TOROIDAL
	integrals_file << "#t F_t E_T" << endl;
#else
	integrals_file << "#t F_t E_T E_Pi E_Pe" ;
	for(int n=1;n<=sim::l;n++){
		integrals_file << " E_" << n;
	}
	integrals_file << endl;
#endif
	return;
}		/* -----  end of function create_integrals_file  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  log_integrals_file
 *  Description:  Adds a line to the integrals file
 * =====================================================================================
 */
	void
log_integrals_file ( double t, double *integrals )
{
#ifdef TOROIDAL
	integrals_file << t << " " << integrals[0] << " " << integrals[1] << endl;
#else
	integrals_file << t << " " << integrals[0] << " " << integrals[1] << " " << integrals[2] << " " << integrals[3];
	//Log energy of multipoles outside the star
	for(int n=1;n<=sim::l;n++){
		integrals_file << " " << integrals[3+n];
	}
	integrals_file << endl;
#endif
	//Release memory
	delete[] integrals;
	integrals=NULL;
	return;
}		/* -----  end of function log_integrals_file  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  close_integrals_file
 *  Description:  
 * =====================================================================================
 */
	void
close_integrals_file ( )
{
	integrals_file.close();
	return;
}		/* -----  end of function close_integrals_file  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  log_field
 *  Description:  Creates log files for the functions alpha and beta at timestep k.
 * =====================================================================================
 */
	void
log_field ( int k )
{
#ifndef TOROIDAL
	//Construct filename for alpha log file, and open it
	stringstream AStream;
	AStream << "results_" << timeStream.str() << "/A_" << k;
	string filenameA=AStream.str();
	ofstream resultsA;
	resultsA.open(filenameA.c_str());
	for(int i=0;i<sim::rNum;i++){
		for(int j=0;j<sim::thNum;j++){
			resultsA << sim::A[i][j] << " ";
		}
		resultsA << endl;
	}
	resultsA.close();
#endif
	//Construct filename for alpha log file, and open it
	stringstream BStream;
	BStream << "results_" << timeStream.str() << "/B_" << k;
	string filenameB=BStream.str();
	ofstream resultsB;
	resultsB.open(filenameB.c_str());
	//Log values
	for(int i=0;i<sim::rNum;i++){
		for(int j=0;j<sim::thNum;j++){
			resultsB << sim::B[i][j] << " ";
		}
		resultsB << endl;
	}
	//Close files
	resultsB.close();

	return;
}		/* -----  end of function log_field  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  report_blowup
 *  Description:  Reports information related to a blowup. The arguments received are the timestep number "k", the radial step number "i" and the angular step number "j" where the abnormality ocurred.
 * =====================================================================================
 */
	void
report_blowup ( int k, int i, int j )
{
	cerr << "Blew up at step " << k << "in place " << i << "," << j << endl; 
	return;
}		/* -----  end of function report_blowup  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  report_completion
 *  Description:  Report information on the completion of the program. This could be just
 *  due to the end of the simulation, numerical blowups, or ctrl+c.
 * =====================================================================================
 */
	void
report_completion ( int info )
{
	if(info==SIGINT){
		cout<<"Interrupt signal captured!"<<endl;
	}else if(info==SIGUSR1){
		cout<<"Couldn't complete simulation!"<<endl;
	}else if(info==SIGUSR2){
		cout<<"Simulation completed normally"<<endl;
	}
	exit(info);
	return;
}		/* -----  end of function report_completion  ----- */

}		/* -----  end of namespace io  ----- */
