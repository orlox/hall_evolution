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
	double rmin=0;
	double k=4.493409457909064;
	double b=-2.125069381043848;

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
	//test for rmin=0, with field outside the star
	return pow(sin(th),2)*(35.0/8.0*pow(r,2)-21.0/4.0*pow(r,4)+15.0/8.0*pow(r,6));
	//test for rmin=0.5
	//return pow(sin(th),2)*(gsl_sf_bessel_jl(1,k*r)+b*gsl_sf_bessel_yl(1,k*r))*r;
	return pow(sin(th),2)*(35.0/8.0*pow(r,2)-21.0/4.0*pow(r,4)+15.0/8.0*pow(r,6));
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
	return 30*(gsl_sf_bessel_jl(1,k*r))*r*sin(th);//+b*gsl_sf_bessel_yl(1,k*r))*gsl_sf_legendre_Plm(1,1,cos(th))*r*sin(th);
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
	return 1;//(1-pow(r,2));
}		/* -----  end of function n  ----- */

#ifndef PUREOHM
#ifndef TOROIDAL
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  n_rderivative
 *  Description:  Gives the radial derivative of n at a given radius and theta
 * =====================================================================================
 */
	double
n_rderivative ( double r, double th )
{
	return 0;//-2*r;
}		/* -----  end of function n_rderivative  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  n_thderivative
 *  Description:  Gives the theta derivative of n at a given radius and theta
 * =====================================================================================
 */
	double
n_thderivative ( double r, double th )
{
	return 0;
}		/* -----  end of function n_thderivative  ----- */
#endif
#endif

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

#ifndef TOROIDAL
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  chi_rderivative
 *  Description:  Gives the value of the r derivative of chi (which goes as -(2*n+r*dn/dr)/(r^3*sin^2(th)*n^2)) at a given radius and theta.
 *
 *  !!!THIS FUNCTION SHOULD NOT BE CHANGED!!!
 * =====================================================================================
 */
	double
chi_rderivative ( double r, double th )
{
	return -(2*n(r,th)+r*n_rderivative(r,th))/(pow(r,3)*pow(sin(th)*n(r,th),2));
}		/* -----  end of function chi_rderivative  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  chi_thderivative
 *  Description:  Gives the value of the theta derivative of chi (which goes as -(2*cos(th)*n+sin(th)*dn/dth)/(r^2*sin^3(th)*n^2)) at a given radius and theta.
 *
 *  !!!THIS FUNCTION SHOULD NOT BE CHANGED!!!
 * =====================================================================================
 */
	double
chi_thderivative ( double r, double th )
{
	return -(2*cos(th)*n(r,th)+sin(th)*n_thderivative(r,th))/(pow(sin(th),3)*pow(r*n(r,th),2));
}		/* -----  end of function chi_thderivative  ----- */
#endif
#endif
}		/* -----  end of namespace initial  ----- */
