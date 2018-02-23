#include <stdio.h>

#define IFMT "%8d"
/*! brief
 *
 * <pre>
 * Output parameters
 * =================
 *   (nzval, rowind, colptr): (*rowind)[*] contains the row subscripts of
 *      nonzeros in columns of matrix A; (*nzval)[*] the numerical values;
 *  column i of A is given by (*nzval)[k], k = (*rowind)[i],...,
 *      (*rowind)[i+1]-1.
 * </pre>
 */
void dreadMM(FILE *fp, int *m, int *n, int *nonz, double **nzval, int **rowind, int **colptr);

///*
// * Read in vector from Matrix Market file.
// *
// * Assumption:
// *   - "dense" Matrix Market Format, i.e., also zero entries are provided.
// */
//int readVec(int A_n, double *val, FILE* fp);
