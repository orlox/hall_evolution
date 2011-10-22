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
}
#endif
