/*
Original work Copyright (c) 2003, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

Modified work Copyright 2018 Technische Universität Dresden, Germany

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
(3) Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * This file is part of SuperLUbench (https://github.com/mflehmig/SuperLUbench.git)
 */


#include "dreadMM.h"

#include <ctype.h>

#include "slu_ddefs.h"

#undef EXPAND_SYM

void dreadMM(FILE *fp, int *m, int *n, int *nonz, double **nzval, int **rowind, int **colptr)
{
  int j, k, jsize, nnz, nz, new_nonz;
  double *a, *val;
  int *asub, *xa, *row, *col;
  int zero_base = 0;
  char *p, line[512], banner[64], mtx[64], crd[64], arith[64], sym[64];
  int expand;

  /* File format:
   *    %%MatrixMarket matrix coordinate real general/symmetric/...
   *    % ...
   *    % (optional comments)
   *    % ...
   *    #rows    #non-zero
   *    Triplet in the rest of lines: row    col    value
   */

  /* 1/ read header */
  fgets(line, 512, fp);
  for (p = line; *p != '\0'; *p = tolower(*p), p++)
    ;

  if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, arith, sym) != 5)
  {
    printf("Invalid header (first line does not contain 5 tokens)\n");
    exit(-1);
  }

  if (strcmp(banner, "%%matrixmarket"))
  {
    printf("Invalid header (first token is not \"%%%%MatrixMarket\")\n");
    exit(-1);
  }

  if (strcmp(mtx, "matrix"))
  {
    printf("Not a matrix; this driver cannot handle that.\n");
    exit(-1);
  }

  if (strcmp(crd, "coordinate"))
  {
    printf("Not in coordinate format; this driver cannot handle that.\n");
    exit(-1);
  }

  if (strcmp(arith, "real"))
  {
    if (!strcmp(arith, "complex"))
    {
      printf("Complex matrix; use zreadMM instead!\n");
      exit(-1);
    }
    else if (!strcmp(arith, "pattern"))
    {
      printf("Pattern matrix; values are needed!\n");
      exit(-1);
    }
    else
    {
      printf("Unknown arithmetic\n");
      exit(-1);
    }
  }

  if (strcmp(sym, "general"))
  {
    printf("Symmetric matrix: will be expanded\n");
    expand = 1;
  }
  else
    expand = 0;

  /* 2/ Skip comments */
  while (banner[0] == '%')
  {
    fgets(line, 512, fp);
    sscanf(line, "%s", banner);
  }

  /* 3/ Read n and nnz */
#ifdef _LONGINT
  sscanf(line, "%ld%ld%ld",m, n, nonz);
#else
  sscanf(line, "%d%d%d", m, n, nonz);
#endif

  if (*m != *n)
  {
    printf("Rectangular matrix!. Abort\n");
    exit(-1);
  }

  if (expand)
    new_nonz = 2 * *nonz - *n;
  else
    new_nonz = *nonz;

  *m = *n;
  printf("m %lld, n %lld, nonz %lld\n", (long long) *m, (long long) *n, (long long) *nonz);
  //MS dallocateA_dist(*n, new_nonz, nzval, rowind, colptr); /* Allocate storage */
  dallocateA(*n, new_nonz, nzval, rowind, colptr); /* Allocate storage */
  a = *nzval;
  asub = *rowind;
  xa = *colptr;

  if (!(val = (double *) SUPERLU_MALLOC(new_nonz * sizeof(double))))
    ABORT("Malloc fails for val[]");
  if (!(row = (int *) SUPERLU_MALLOC(new_nonz * sizeof(int))))
    ABORT("Malloc fails for row[]");
  if (!(col = (int *) SUPERLU_MALLOC(new_nonz * sizeof(int))))
    ABORT("Malloc fails for col[]");

  for (j = 0; j < *n; ++j)
    xa[j] = 0;

  /* 4/ Read triplets of values */
  for (nnz = 0, nz = 0; nnz < *nonz; ++nnz)
  {
#ifdef _LONGINT
    fscanf(fp, "%lld%lld%lf\n", &row[nz], &col[nz], &val[nz]);
#else
    fscanf(fp, "%d%d%lf\n", &row[nz], &col[nz], &val[nz]);
#endif

    if (nnz == 0) /* first nonzero */
      if (row[0] == 0 || col[0] == 0)
      {
        zero_base = 1;
        printf("triplet file: row/col indices are zero-based.\n");
      }
      else
        printf("triplet file: row/col indices are one-based.\n");

    if (!zero_base)
    {
      /* Change to 0-based indexing. */
      --row[nz];
      --col[nz];
    }

    if (row[nz] < 0 || row[nz] >= *m || col[nz] < 0 || col[nz] >= *n /*|| val[nz] == 0.*/)
    {
      fprintf(stderr, "nz " IFMT ", (" IFMT ", " IFMT ") = %e out of bound, removed\n", nz, row[nz], col[nz], val[nz]);
      exit(-1);
    }
    else
    {
      ++xa[col[nz]];
      if (expand)
      {
        if (row[nz] != col[nz])
        { /* Excluding diagonal */
          ++nz;
          row[nz] = col[nz - 1];
          col[nz] = row[nz - 1];
          val[nz] = val[nz - 1];
          ++xa[col[nz]];
        }
      }
      ++nz;
    }
  }

  *nonz = nz;
  if (expand)
  {
    printf("new_nonz after symmetric expansion:\t" IFMT "\n", *nonz);
  }

  /* Initialize the array of column pointers */
  k = 0;
  jsize = xa[0];
  xa[0] = 0;
  for (j = 1; j < *n; ++j)
  {
    k += jsize;
    jsize = xa[j];
    xa[j] = k;
  }

  /* Copy the triplets into the column oriented storage */
  for (nz = 0; nz < *nonz; ++nz)
  {
    j = col[nz];
    k = xa[j];
    asub[k] = row[nz];
    a[k] = val[nz];
    ++xa[j];
  }

  /* Reset the column pointers to the beginning of each column */
  for (j = *n; j > 0; --j)
    xa[j] = xa[j - 1];
  xa[0] = 0;

  SUPERLU_FREE(val);
  SUPERLU_FREE(row);
  SUPERLU_FREE(col);

#ifdef CHK_INPUT
  int i;
  for (i = 0; i < *n; i++)
  {
    printf("Col %d, xa %d\n", i, xa[i]);
    for (k = xa[i]; k < xa[i+1]; k++)
    printf("%d\t%16.10f\n", asub[k], a[k]);
  }
#endif

}

static void dreadrhs(int m, double *b)
{
  FILE *fp, *fopen();
  int i;

  if (!(fp = fopen("b.dat", "r")))
  {
    fprintf(stderr, "dreadrhs: file does not exist\n");
    exit(-1);
  }
  for (i = 0; i < m; ++i)
    fscanf(fp, "%lf\n", &b[i]);
  /*fscanf(fp, "%d%lf\n", &j, &b[i]);*/
  /*        readpair_(j, &b[i]);*/
  fclose(fp);
}

//
//int readVec(int A_n, double *val, FILE* fp)
//{
//  int ret = 0;
//  int nnz, m, n;
//  char *p, line[512], banner[64], mtx[64], crd[64], arith[64], sym[64];
//  int expand;
//
//  /* File format:
//   *    %%MatrixMarket matrix coordinate real general/symmetric/...
//   *    % ...
//   *    % (optional comments)
//   *    % ...
//   *    #rows    #non-zero
//   *    Triplet in the rest of lines: row    col    value
//   */
//
//  /* 1/ read header */
//  fgets(line, 512, fp);
//  for (p = line; *p != '\0'; *p = tolower(*p), p++)
//    ;
//
//  if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, arith, sym) != 5)
//  {
//    printf("Invalid header (first line does not contain 5 tokens)\n");
//    exit(-1);
//  }
//
//  if (strcmp(banner, "%%matrixmarket"))
//  {
//    printf("Invalid header (first token is not \"%%%%MatrixMarket\")\n");
//    exit(-1);
//  }
//
//  if (strcmp(mtx, "matrix"))
//  {
//    printf("Not a matrix; this driver cannot handle that.\n");
//    exit(-1);
//  }
//
//  if (strcmp(crd, "coordinate"))
//  {
//    printf("Not in coordinate format; this driver cannot handle that.\n");
//    exit(-1);
//  }
//
//  if (strcmp(arith, "real"))
//  {
//    if (!strcmp(arith, "complex"))
//    {
//      printf("Complex matrix; use zreadMM instead!\n");
//      exit(-1);
//    }
//    else if (!strcmp(arith, "pattern"))
//    {
//      printf("Pattern matrix; values are needed!\n");
//      exit(-1);
//    }
//    else
//    {
//      printf("Unknown arithmetic\n");
//      exit(-1);
//    }
//  }
//
//  if (strcmp(sym, "general"))
//  {
//    printf("Symmetric matrix: will be expanded\n");
//    printf("Not implemented yet. Error.\n");
//    expand = 1;
//    exit(-3);
//  }
//
//  /* 2/ Skip comments */
//  while (banner[0] == '%')
//  {
//    fgets(line, 512, fp);
//    sscanf(line, "%s", banner);
//  }
//
//  /* 3/ Read n and nnz */
//#ifdef _LONGINT
//  sscanf(line, "%ld%ld%ld", &m, &n, &nnz);
//#else
//  sscanf(line, "%d%d%d", &m, &n, &nnz);
//#endif
//
//  // Compare dimension of vector and matrix. Do they match?
//  if (m != A_n)
//  {
//    printf("Dimension of vector and matrix do not match!. Abort\n");
//    exit(-1);
//  }
//  // Compare dimension of vector and matrix. Do they match?
//  if (m != nnz)
//  {
//    printf("NNZ != dimension!. Abort\n");
//    exit(-1);
//  }
//  // Is it really a vector, i.e., second dimension is equal to 1?
//  if (n != 1)
//  {
//    printf("Not a vector!. Abort\n");
//    exit(-1);
//  }
//
//  /* 4/ Read triplets of values */
//  int row, col, nz;
//  // Read first line and check for 1-based indices.
//#ifdef _LONGINT
//  fscanf(fp, "%lld%lld%lf\n", &row, &col, &val[0]);
//#else
//  fscanf(fp, "%d%d%lf\n", &row, &col, &val[0]);
//#endif
//  if (row == 0 || col == 0)
//  {
//    printf("triplet file: row/col indices are zero-based.\n");
//    exit(-3);
//  }
//  else
//  {
//    printf("Info: triplet file: row/col indices are one-based.\n");
//  }
//  // Read the rest of the vector.
//  for (nz = 1; nz < nnz; ++nz)
//  {
//#ifdef _LONGINT
//    fscanf(fp, "%lld%lld%lf\n", &row, &col, &val[nz]);
//#else
//    fscanf(fp, "%d%d%lf\n", &row, &col, &val[nz]);
//#endif
//    // +1 because of 1-based indices
//    if (row < 0 || row >= m + 1 || col != 1)
//    {
//      fprintf(stderr, "nz " IFMT ", (" IFMT ", " IFMT ") = %e out of bound, removed\n", nz, row, col, *val);
//      exit(-1);
//    }
//  }
//
//  return ret;
//}
//
