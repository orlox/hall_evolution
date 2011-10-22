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

}		/* -----  end of namespace io  ----- */
