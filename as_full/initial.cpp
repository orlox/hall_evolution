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
double k=6.572013199016351;
double b=-2.125069381043848;
#ifndef TOROIDAL
	double
A ( double r, double th )
{
	//test for rmin=0
	//return pow(sin(th),2)*(9.0/5.0*pow(r,2)-14.0/5.0*pow(r,4)+pow(r,6));
	//test for rmin=0.5
	return pow(sin(th),2)*(gsl_sf_bessel_jl(1,k*r)+b*gsl_sf_bessel_yl(1,k*r))*r;
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
	//return 0;
	return 20*(gsl_sf_bessel_jl(1,k*r)+b*gsl_sf_bessel_yl(1,k*r))*gsl_sf_legendre_Plm(1,1,cos(th))*r*sin(th);
}		/* -----  end of function B  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  rho
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

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  chi
 *  Description:  Gives the value of chi (which goes as 1/(r^2*sin^2(th)*n(r,th))) at a given radius and theta.
 *
 *  !!!THIS FUNCTION SHOULD NOT BE CHANGED!!!
 * =====================================================================================
 */
#ifndef PUREOHM
	double
chi ( double r, double th )
{
	return 1/(pow(r*sin(th),2)*n(r,th));
}		/* -----  end of function chi  ----- */
#endif
}		/* -----  end of namespace initial  ----- */
