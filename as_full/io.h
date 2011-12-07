/*
 * =====================================================================================
 *
 *       Filename:  io.h
 *
 *    Description:  Header file for io.cpp, where io operations are performed
 *
 *        Version:  1.0
 *        Created:  10/21/2011 08:19:23 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Pablo Marchant
 *        Company:  Pontificia Universidad Católica de Chile
 *
 * =====================================================================================
 */
#ifndef IO_H
#define IO_H
namespace io{
//Read arguments from command line.
int read_args(int argc, char *argv[]);
//Functions to check if string is valid int or float
int is_integer(char *value);
int is_float(char *value);
//Creates folder where results will be logged.
void create_folder();
//Prints block of text describing simulation to be done.
void print_header();
//Create file where integrated quantities are logged
void create_integrals_file();
//Prints a line in file with integrated quantities
void log_integrals_file(double t, double*B);
//Close file where integrated quantities are logged
void close_integrals_file();
//Creates two files for the complete profiles of alpha and beta
void log_field(int k);
//Prints info relating to a blowup
void report_blowup(int k,int i,int j);
//Prints info on the termination of the program
void report_completion(int info);
}
#endif
