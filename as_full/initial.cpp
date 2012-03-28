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
	double rmin=0.75;
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
	//#######################OHM EIGENMODES#######################
	//The following conditions correspond to Ohm eigenmodes
	//for constant resistivity, and are normalized so B_pol
	//is very close to 1 at its maximun. This values require
	//r_min=0.75.
	//MODE 1,1
	int l=1;
	double A=-0.55882;
	double B=-0.52004;
	double k=7.03266;
	//MODE 2,1
//	int l=1;
//	double A=-0.72288;
//	double B=-0.20659;
//	double k=19.12793;
	//MODE 1,2
//	int l=2;
//	double A=-0.52101;
//	double B=-0.04764;
//	double k=7.81795;
	//MODE 2,2
//	int l=2;
//	double A=-0.31188;
//	double B=0.39534;
//	double k=19.46616;

	return 1*r*(A*gsl_sf_bessel_jl(l,k*r)+B*gsl_sf_bessel_yl(l,k*r))*gsl_sf_legendre_Plm(l,1,cos(th))*sin(th);
	//######################END OHM EIGENMODES####################

	//Equilibrium field due to rigid body rotation of electrons for constant electron density
//	return 1/2.0*pow(sin(th),2)*((3*pow(rmin,5)-5*pow(rmin,3))/r+5*pow(r,2)-3*pow(r,4))/(2-5*pow(rmin,3)+3*pow(rmin,5));
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
	//#######################OHM EIGENMODES#######################
	//The following conditions correspond to Ohm eigenmodes
	//for constant resistivity, and are normalized so B_tor
	//is very close to 1 at its maximun. These values require
	//r_min=0.75.
	//MODE 1,1
	int l=1;
	double A=-2.00235;
	double B=10.81346;
	double k=12.67071;
	//MODE 2,1
//	int l=1;
//	double A=-1.88596;
//	double B=20.32767;
//	double k=25.18557;
	//MODE 1,2
//	int l=2;
//	double A=6.32868;
//	double B=3.82300;
//	double k=12.87682;
	//MODE 2,2
//	int l=2;
//	double A=13.11642;
//	double B=3.72583;
//	double k=25.29089;

	return r*(A*gsl_sf_bessel_jl(l,k*r)+B*gsl_sf_bessel_yl(l,k*r))*gsl_sf_legendre_Plm(l,1,cos(th))*sin(th);
	//######################END OHM EIGENMODES####################
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
