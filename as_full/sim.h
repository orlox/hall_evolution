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
//physical values that describe the structure of the star.
extern double **chi;
extern double **eta;
//arrays that contains the precalculated quantities neccesary to solve fluxes for B
extern double **hall_rflux;
extern double **hall_thflux;
extern double **res_rflux;
extern double **res_thflux;
//sines precalculated at each point in the grid
extern double *sines;
//minimun radius of the shell containing the magnetic field.
extern double rmin;
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
}
#endif
