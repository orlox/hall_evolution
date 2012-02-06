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
extern double dt;
extern int tNum;
extern int plotSteps;
extern double thtd;
//physical values of functions defining the magnetic field.
extern double **B;
#ifndef TOROIDAL
extern double **A;
#endif
//minimun radius of the shell containing the magnetic field.
extern double rmin;
//Number of timesteps for which only Ohmic dissipation will be in effect. This enhances
//stability significantly.
extern int steps_less;
#ifndef TOROIDAL
//Number of points used for the multipole fit outside the star
extern int l;
#endif
//size of spatial steps.
extern double dr;
extern double dth;
//function to set initial conditions
void initial_conditions();
//function that solves values that remain the same through all the simulation
void solve_repeated_values();
//function that performs the simulation
int simulate();
//function that solves integrated quantities
double* solve_integrals();
//function that solves integrated quantities
void solve_A_boundary();
//function releases memory when simulation is done
void release_memory(int signal);
}
#endif
