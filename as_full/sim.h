/*
 * =====================================================================================
 *
 *       Filename:  sim.h
 *
 *    Description:  Header file for sim.cpp where all the calculations related to the simulation are done
 *
 *        Version:  1.0
 *        Created:  10/21/2011 09:10:35 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Pablo Marchant
 *        Company:  Pontificia Universidad Católica de Chile
 *
 * =====================================================================================
 */
#ifndef SIM_H
#define SIM_H
namespace sim{
//values given as cli arguments.
extern int rNum;
extern int thNum;
extern double factor;
extern int tNum;
extern int plotSteps;
extern double thtd;
//physical values of functions defining the magnetic field.
extern double **B, **dBr, **dBth;
#ifndef TOROIDAL
	extern double **A, **Aaux, **gsA;
#endif
//minimun radius of the shell containing the magnetic field.
extern double rmin;
#ifndef TOROIDAL
	//Number of points used for the multipole fit outside the star
	extern int l;
	//Values of the coefficients that give the poloidal field outside the star.
	extern double *a;
#endif
//size of spatial (and temporal) steps.
extern double dr;
extern double dth;
extern double dt;
//function to set initial conditions
void initial_conditions();
//function that solves values that remain the same through all the simulation
void solve_repeated_values();
//function that performs the simulation
int simulate();
//function that solves new values of A
void solve_new_A();
//function that solves new values of B
void solve_new_B();
//function that solves integrated quantities
double* solve_integrals();
//function that solves boundary conditions on alpha
void solve_A_boundary();
//function that solves boundary conditions on beta
void solve_B_boundary();
//function releases memory when simulation is done
void release_memory(int signal);
}
#endif
