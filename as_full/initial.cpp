/*
 * =====================================================================================
 *
 *       Filename:  initial.cpp
 *
 *    Description:  This file contains the initial conditions for the simulation. In order to consider different initial conditions, this file must be modified and recompiled. All these functions are contained under the namespace "initial".
 *
 *    In order for the ratio t_h/t_D to make any sense, all quantities given here should have a characteristic value of 1.
 *
 *        Version:  1.0
 *        Created:  10/22/2011 01:17:38 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Pablo Marchant
 *        Company:  Pontificia Universidad Católica de Chile
 *
 * =====================================================================================
 */
#include "initial.h"
#include <math.h>
#include <gsl/gsl_sf.h>

namespace initial{
	//minimun radius of the shell containing the magnetic field
	double rmin=0.5;
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  A
 *  Description:  Gives the initial values of the function alpha that describes the poloidal field at a given radius and theta.
 *
 *  BOUNDARY CONDITIONS!!: A should be 1 in the axis of symmetry, and it should be such that the magnetic field is completely continuous along the surface of the star (i.e. no surface currents).
 * =====================================================================================
 */
#ifndef TOROIDAL
	double
A ( double r, double th )
{
	double A=.1501063239606539;
	double B=.1954706187148793;
	double k=4.057515676220868;
	return 10*r*(A*gsl_sf_bessel_jl(1,k*r)+B*gsl_sf_bessel_yl(1,k*r))*gsl_sf_legendre_Plm(1,1,cos(th))*sin(th);
	//return 0.1*(-480.0/47*(pow(r,4)/10-pow(r,2)/6)-17.0/94/r)*pow(sin(th),2);
}		/* -----  end of function Ai  ----- */
#endif

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  B
 *  Description:  Gives the initial values of the function beta that describes the toroidal field at a given radius and theta.
 *
 *  BOUNDARY CONDITIONS!!: B should be 0 in the axis of symmetry (because currents must remain bounded there) and at the surface of the star (to avoid surface currents).
 * =====================================================================================
 */
	double
B ( double r, double th )
{
	double k=6.572013199016351;
	double b=-2.12506938104384;
	return 3*(gsl_sf_bessel_jl(1,k*r)+b*gsl_sf_bessel_yl(1,k*r))*gsl_sf_legendre_Plm(1,1,cos(th))*r*sin(th);
}		/* -----  end of function B  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  n
 *  Description:  Gives the value of the electron density at a given radius and theta.
 * =====================================================================================
 */
	double
n ( double r, double th )
{
	return 1;
}		/* -----  end of function n  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  eta
 *  Description:  Gives the value of the resistivity at a given radius and theta.
 * =====================================================================================
 */
	double
eta ( double r, double th )
{
	return 1;
}		/* -----  end of function eta  ----- */

#ifndef PUREOHM
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  chi
 *  Description:  Gives the value of chi (which goes as 1/(r^2*sin^2(th)*n(r,th))) at a given radius and theta.
 *
 *  !!!THIS FUNCTION SHOULD NOT BE CHANGED!!!
 * =====================================================================================
 */
	double
chi ( double r, double th )
{
	return 1/(pow(r*sin(th),2)*n(r,th));
}		/* -----  end of function chi  ----- */
#endif
}		/* -----  end of namespace initial  ----- */
