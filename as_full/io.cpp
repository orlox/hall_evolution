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
#include	"string.h"
#include	"io.h"
#include	"sim.h"

using namespace std;
namespace io{
//std::stringstream used to store the timestamp
stringstream timeStream;
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
	if(argc!=7){
		cerr<<"ERROR: Not enough arguments provided"<<endl;
		cerr<<"usage is:"<<endl;
		cerr<<"as_full_run <rNum(int)> <thNum(int)> <dt(float)> <tNum(int)> <tNumPlot(int)> <th/td(float)>"<<endl;
		return 1;
	}
	//check integrity of all parameters
	if(!io::is_integer(argv[1])||
			!io::is_integer(argv[2])||
			!io::is_float(argv[3])||
			!io::is_integer(argv[4])||
			!io::is_integer(argv[5])||
			!io::is_float(argv[6])){
		cerr<<"ERROR: Arguments not in the right format"<<endl;
		cerr<<"usage is:"<<endl;
		cerr<<"as_full_run <rNum(int)> <thNum(int)> <dt(float)> <tNum(int)> <tNumPlot(int)> <th/td(float)>"<<endl;
		return 1;
	}
	//store values for simulation
	sim::rNum=atoi(argv[1]);
	sim::thNum=atoi(argv[2]);
	sim::dt=atof(argv[3]);
	sim::tNum=atoi(argv[4]);
	sim::plotSteps=atoi(argv[5]);
	sim::thtd=atof(argv[6]);
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
	//check that each individual character is an integer and return 0 if one is not
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
	//check that each individual character is an integer and return 0 if one is not.
	//Also accepts a single point.
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
 *  Description:  Creates folder using as a distinct id the current timestamp. It also creates the file params.dat which stores the simulation parameters
 * =====================================================================================
 */
	void
create_folder ( )
{
	//store timestamp, will be needed also later to access folder
	timeStream << time(0);
	//create directory. THIS WONT WORK ON WINDOWS!!
	string mkdir="mkdir results_"+timeStream.str();
	system(mkdir.c_str());
	//create params.dat file, where simulation parameters are logged
	ofstream params;
	string filename="results_"+timeStream.str()+"/params.dat";
	params.open(filename.c_str());
	params<< "rNum:"<< sim::rNum << std::endl;
	params<< "thNum:"<< sim::thNum << std::endl;
	params<< "dt:"<< sim::dt << std::endl;
	params<< "tNum:"<< sim::tNum << std::endl;
	params<< "plotSteps:"<< sim::plotSteps << std::endl;
	params<< "rmin:"<< sim::rmin << std::endl;
	params<< "thtd:"<< sim::thtd << std::endl;
	params.close();

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
	cout<< endl;
	cout<< "#################################################"<< endl;
	cout<< "##########Axisymmetric Hall Evolution############"<< endl;
	cout<< "#################################################"<< endl;
	cout<< "PARAMETERS CHOSEN:"<< endl;
	cout<< "-Number of Radial steps: "<< sim::rNum << endl;
	cout<< "-Number of Angular steps: "<< sim::thNum << endl;
	cout<< "-Size of timestep: "<< sim::dt << endl;
	cout<< "-Number of Time steps: "<< sim::tNum << endl;
	cout<< "-Save output every "<< sim::plotSteps << " steps" << endl;
	cout<< "-Minimun radius: "<< sim::rmin << endl;
	cout<< "-Ratio of hall to dissipative timescales: "<< sim::thtd << endl;
	cout<< endl;
	cout<< "Results stored in folder results_" << timeStream.str() << endl;
	cout<< endl;
	return;
}		/* -----  end of function print_header  ----- */

}		/* -----  end of namespace io  ----- */
