/*
 * =====================================================================================
 *
 *       Filename:  as_full.cpp
 *
 *    Description:  Main file of the program. This code follows the evolution of an axisymmetric field due to the hall effect.
 *
 *        Version:  1.0
 *        Created:  10/21/2011 08:09:09 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Pablo Marchant
 *        Company:  Pontifica Universidad Católica de Chile
 *
 * =====================================================================================
 */

#include	"io.h"
#include	"sim.h"
#include <signal.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  No description neccesary
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	//Read arguments from cli. If it fails, exit with error code 1!
	if(io::read_args(argc, argv))
		return 1;

	//set up the initial conditions
	sim::initial_conditions();

	//capture ctrl+c to avoid memory leaks
	signal(SIGINT, sim::release_memory);

	//create folder where results will be logged
	io::create_folder();

	//put in std::cout a block of text describing the simulation to be done
	io::print_header();

	//perform the simulation. If it fails, exit with error code 2!
	if(sim::simulate()){
		//release memory
		sim::release_memory(SIGUSR1);
	}
	//release memory
	sim::release_memory(SIGUSR2);
}				/* ----------  end of function main  ---------- */
