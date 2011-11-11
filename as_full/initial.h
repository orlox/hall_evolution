/*
 * =====================================================================================
 *
 *       Filename:  initial.h
 *
 *    Description:  Header file for initial.cpp, where the initial conditions for the simulation are set.
 *
 *        Version:  1.0
 *        Created:  10/22/2011 01:12:58 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Pablo Marchant
 *        Company:  Pontificia Universidad Católica de Chile
 *
 * =====================================================================================
 */
#ifndef INITIAL_H
#define INITIAL_H
namespace initial{
	//value for the inner radius of the shell where the field is contained
	extern double rmin;
	//functions that give the initial values of the magnetic field
	double A( double r, double th );
	double B( double r, double th );
	//functions that describe the structure of the star.
	//chi(r,th) should not be modified.
	double n( double r, double th );
	double n_rderivative( double r, double th );
	double n_thderivative( double r, double th );
	double chi( double r, double th );
	double chi_rderivative( double r, double th );
	double chi_thderivative( double r, double th );
	double eta( double r, double th );
}
#endif
