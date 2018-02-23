/*
 * Util.h
 *
 *  Created on: 21 Feb 2018
 *      Author: mf
 */

#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_

#define IFMT "%8d"

/*
 * Compares the computed solution vector resX against the known vector expX.
 */
int compareSolution(double* resX, int m, double* expX, double eps);

/// Time in micro seconds.
long long current_timestamp();

/*
 * Read in vector from Matrix Market file.
 *
 * Assumption:
 *   - "dense" Matrix Market Format, i.e., also zero entries are provided.
 */
int readVec(int A_n, double *val, FILE* fp);

#endif /* SRC_UTIL_H_ */
