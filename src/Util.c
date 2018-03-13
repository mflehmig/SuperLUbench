/*
This file is part of SuperLUbench (https://github.com/mflehmig/SuperLUbench.git)

Copyright 2018 Technische Universität Dresden, Germany

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
(3) Neither the name of Lawrence Berkeley National Laboratory, U.S. Dept. of
Energy nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include "Util.h"


int compareSolution(double* resX, int m, double* expX, double eps)
{
  int ret = 0;
  double f;
  for (unsigned int i = 0; i < m; i++)
  {
//    if ((resX[i] < (expX[i] - eps)) || (resX[i] > (expX[i] + eps)))
//    {
//      printf("unexpected answer: expected(%ud)=%lg vs result=%lg\n", i, expX[i], resX[i]);
//      ret = 1;
//    }
    f = fabs(resX[i] - expX[i]);
    if (f > eps) {
      printf("\tunexpected answer: expected(%ud)=%.10f vs result=%.10f\n", i, expX[i], resX[i]);
      ret = 1;
    }
  }
  return ret;
}

/// Time in micro seconds.
long long current_timestamp()
{
  struct timeval te;
  gettimeofday(&te, NULL);
  long long usecs = te.tv_sec * 1000000LL + te.tv_usec;
  return usecs;
}

int readVec(int A_n, double *val, FILE* fp)
{
  int ret = 0;
  int nnz, m, n;
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
    printf("Not implemented yet. Error.\n");
    expand = 1;
    exit(-3);
  }

  /* 2/ Skip comments */
  while (banner[0] == '%')
  {
    fgets(line, 512, fp);
    sscanf(line, "%s", banner);
  }

  /* 3/ Read n and nnz */
#ifdef _LONGINT
  sscanf(line, "%ld%ld%ld", &m, &n, &nnz);
#else
  sscanf(line, "%d%d%d", &m, &n, &nnz);
#endif

  // Compare dimension of vector and matrix. Do they match?
  if (m != A_n)
  {
    printf("Dimension of vector and matrix do not match!. Abort\n");
    exit(-1);
  }
  // Compare dimension of vector and matrix. Do they match?
  if (m != nnz)
  {
    printf("NNZ != dimension!. Abort\n");
    exit(-1);
  }
  // Is it really a vector, i.e., second dimension is equal to 1?
  if (n != 1)
  {
    printf("Not a vector!. Abort\n");
    exit(-1);
  }

  /* 4/ Read triplets of values */
  int row, col, nz;
  // Read first line and check for 1-based indices.
#ifdef _LONGINT
  fscanf(fp, "%lld%lld%lf\n", &row, &col, &val[0]);
#else
  fscanf(fp, "%d%d%lf\n", &row, &col, &val[0]);
#endif
  if (row == 0 || col == 0)
  {
    printf("triplet file: row/col indices are zero-based.\n");
    exit(-3);
  }
  else
  {
    printf("Info: triplet file: row/col indices are one-based.\n");
  }
  // Read the rest of the vector.
  for (nz = 1; nz < nnz; ++nz)
  {
#ifdef _LONGINT
    fscanf(fp, "%lld%lld%lf\n", &row, &col, &val[nz]);
#else
    fscanf(fp, "%d%d%lf\n", &row, &col, &val[nz]);
#endif
    // +1 because of 1-based indices
    if (row < 0 || row >= m + 1 || col != 1)
    {
      fprintf(stderr, "nz " IFMT ", (" IFMT ", " IFMT ") = %e out of bound, removed\n", nz, row, col, *val);
      exit(-1);
    }
  }

  return ret;
}
