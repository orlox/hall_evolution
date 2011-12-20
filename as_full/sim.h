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
//Lost energy calculated from an analytical expression. Should match the real loss.
extern double lost_energy;
//minimun radius of the shell containing the magnetic field.
extern double rmin;
//Number of points in the radial direction both at the surface and the inner boundary for which the values
//will be solved by interpolation, and not direct calculation through the time evolution equation.
extern int rless;
#ifndef TOROIDAL
#ifndef SIMPLE
//Number of points used for the multipole fit outside the star
extern int l;
#endif
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
