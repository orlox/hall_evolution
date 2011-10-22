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

#include	<stdlib.h>
#include	"io.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  No description neccesary
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	io::read_args(argc, argv);
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
