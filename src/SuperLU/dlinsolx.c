/*! \file
 Copyright (c) 2003, The Regents of the University of California, through
 Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from U.S. Dept. of Energy)

 All rights reserved.

 The source code is distributed under BSD license, see the file License.txt
 at the top-level directory.
 */

/*
 * -- SuperLU routine (version 3.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 1, 2008
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <slu_ddefs.h>
#include "../dreadMM.h"
#include "../Util.h"

void parse_command_line(int argc, char *argv[], int *lwork, double *u, yes_no_t *equil, trans_t *trans, char *fileA,
                        char *fileB, char *fileX, double* eps);
colperm_t getSuperLUOrdering();

int main(int argc, char *argv[])
{
  char equed[1];
  yes_no_t equil;
  trans_t trans;
  SuperMatrix A, L, U;
  SuperMatrix B, X;
  NCformat *Astore;
  NCformat *Ustore;
  SCformat *Lstore;
  GlobalLU_t Glu; /* facilitate multiple factorizations with
   SamePattern_SameRowPerm                  */
  double *a;
  int *asub, *xa;
  int *perm_r; /* row permutations from partial pivoting */
  int *perm_c; /* column permutation vector */
  int *etree;
  void *work;
  int info, lwork, nrhs, ldx;
  int i, m, n, nnz;
  double *rhsb, *rhsx, *xact;
  double *R, *C;
  double *ferr, *berr;
  double u, rpg, rcond;
  mem_usage_t mem_usage;
  superlu_options_t options;
  SuperLUStat_t stat;
  FILE *fp;
  char fileA[128];
  char fileB[128] = "";
  char fileX[128] = "";
  double eps = 5e-8;
  int permc_spec;

  /* Defaults */
  lwork = 0;
  nrhs = 1;
  equil = YES;
  u = 1.0;
  trans = NOTRANS;

  /* Set the default input options:
   options.Fact = DOFACT;
   options.Equil = YES;
   options.ColPerm = COLAMD;
   options.DiagPivotThresh = 1.0;
   options.Trans = NOTRANS;
   options.IterRefine = NOREFINE;
   options.SymmetricMode = NO;
   options.PivotGrowth = NO;
   options.ConditionNumber = NO;
   options.PrintStat = YES;
   */
  set_default_options(&options);
  options.ColPerm = getSuperLUOrdering();

  /* Can use command line input to modify the defaults. */
  parse_command_line(argc, argv, &lwork, &u, &equil, &trans, fileA, fileB, fileX, &eps);
  options.Equil = equil;
  options.DiagPivotThresh = u;
  options.Trans = trans;

  /* Add more functionalities that the defaults. */
  options.PivotGrowth = YES; /* Compute reciprocal pivot growth */
  options.ConditionNumber = YES;/* Compute reciprocal condition number */
  options.IterRefine = SLU_DOUBLE; /* Perform double-precision refinement */

  if (lwork > 0)
  {
    work = SUPERLU_MALLOC(lwork);
    if (!work)
    {
      ABORT("DLINSOLX: cannot allocate work[]");
    }
  }

  /* Read matrix A from a file in Harwell-Boeing format.*/
  // dreadhb(fp, &m, &n, &nnz, &a, &asub, &xa);
  fp = fopen(fileA, "r");
  if (fp == NULL)
  {
    perror("input error: could not open file\n");
    exit(-2);
  }
  dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
  fclose(fp);

  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
  Astore = A.Store;
  printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);

  if (!(rhsb = doubleMalloc(m * nrhs)))
    ABORT("Malloc fails for rhsb[].");
  if (!(rhsx = doubleMalloc(m * nrhs)))
    ABORT("Malloc fails for rhsx[].");
  dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
  dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
  xact = doubleMalloc(n * nrhs);
  ldx = n;

  /* Read matrix b and x from a file in Matrix Market format.*/
  if (strcmp(fileB, "") != 0 && strcmp(fileX, "") != 0)
  {
    printf("Read in b and x from provided files.\n");
    fp = fopen(fileB, "r");
    if (fp == NULL)
    {
      perror("input error: could not open file containing b\n");
      exit(-2);
    }
    readVec(n, rhsb, fp);
    fclose(fp);

    fp = fopen(fileX, "r");
    if (fp == NULL)
    {
      perror("input error: could not open file containing x\n");
      exit(-2);
    }
    readVec(n, xact, fp);
    fclose(fp);
  }
  else
  {
    // x = [1, 1, ...., 1]
    printf("Call dGenXTrue and dFillRHS.\n");
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(trans, nrhs, xact, ldx, &A, &B);
  }

  if (!(etree = intMalloc(n)))
    ABORT("Malloc fails for etree[].");
  if (!(perm_r = intMalloc(m)))
    ABORT("Malloc fails for perm_r[].");
  if (!(perm_c = intMalloc(n)))
    ABORT("Malloc fails for perm_c[].");
  if (!(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))))
    ABORT("SUPERLU_MALLOC fails for R[].");
  if (!(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))))
    ABORT("SUPERLU_MALLOC fails for C[].");
  if (!(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))))
    ABORT("SUPERLU_MALLOC fails for ferr[].");
  if (!(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))))
    ABORT("SUPERLU_MALLOC fails for berr[].");

  /* Initialize the statistics variables. */
  StatInit(&stat);

  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = 0: natural ordering
   *   permc_spec = 1: minimum degree ordering on structure of A'*A
   *   permc_spec = 2: minimum degree ordering on structure of A'+A
   *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
   */
  permc_spec = 1;
  get_perm_c(permc_spec, &A, perm_c);

  /* Solve the system and compute the condition number
   and error bounds using dgssvx.      */
  long long start = current_timestamp();
  dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr, &Glu,
         &mem_usage, &stat, &info);
  long long finish = current_timestamp();
  printf("dgssvx(): info %d\n", info);

  if (info == 0 || info == n + 1)
  {
    if (options.PivotGrowth == YES)
      printf("Recip. pivot growth = %e\n", rpg);
    if (options.ConditionNumber == YES)
      printf("Recip. condition number = %e\n", rcond);
    if (options.IterRefine != NOREFINE)
    {
      printf("Iterative Refinement:\n");
      printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
      for (i = 0; i < nrhs; ++i)
        printf("%8d%8d%16e%16e\n", i + 1, stat.RefineSteps, ferr[i], berr[i]);
    }
    Lstore = (SCformat *) L.Store;
    Ustore = (NCformat *) U.Store;
    printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
    printf("FILL ratio = %.1f\n", (float) (Lstore->nnz + Ustore->nnz - n) / nnz);

    printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu / 1e6, mem_usage.total_needed / 1e6);
    printf("Time for pdgssvx() [1e-6 s]: %lli\n", (finish - start));
    fflush(stdout);

    /* This is how you could access the solution matrix. */
    // double *sol = (double*) ((DNformat*) X.Store)->nzval;
    dinf_norm_error(nrhs, &X, xact);
  }
  else if (info > 0 && lwork == -1)
    printf("** Estimated memory: %d bytes\n", info - n);
  else
    printf("ERROR: dgssvx() returns info %d\n", info);

  if (options.PrintStat)
    StatPrint(&stat);
  StatFree(&stat);

  SUPERLU_FREE(rhsb);
  SUPERLU_FREE(rhsx);
  SUPERLU_FREE(xact);
  SUPERLU_FREE(etree);
  SUPERLU_FREE(perm_r);
  SUPERLU_FREE(perm_c);
  SUPERLU_FREE(R);
  SUPERLU_FREE(C);
  SUPERLU_FREE(ferr);
  SUPERLU_FREE(berr);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperMatrix_Store(&X);
  if (lwork == 0)
  {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
  }
  else if (lwork > 0)
    SUPERLU_FREE(work);

}

/*
 * Parse command line inputs.
 */
void parse_command_line(int argc, char *argv[], int *lwork, double *u, yes_no_t *equil, trans_t *trans, char *fileA,
                        char *fileB, char *fileX, double* eps)
{
  int c;
  extern char *optarg;

  while ((c = getopt(argc, argv, "hl:u:e:t:A:b:x:a:")) != EOF)
  {
    switch (c)
    {
      case 'h':
        printf("Options:\n");
        printf("\t-l <int> - length of work[*] array\n");
        printf("\t-u <int> - pivoting threshold\n");
        printf("\t-e <0 or 1> - equilibrate or not\n");
        printf("\t-t <0 or 1> - solve transposed system or not\n");
        exit(1);
        break;
      case 'l':
        *lwork = atoi(optarg);
        break;
      case 'u':
        *u = atof(optarg);
        break;
      case 'e':
        *equil = atoi(optarg);
        break;
      case 't':
        *trans = atoi(optarg);
        break;
      case 'A':  // file with matrix A
        memcpy(fileA, optarg, strlen(optarg) + 1);
        break;
      case 'b':  // file with matrix b
        memcpy(fileB, optarg, strlen(optarg) + 1);
        break;
      case 'x':  // file with matrix x
        memcpy(fileX, optarg, strlen(optarg) + 1);
        break;
      case 'a':
        *eps = atof(optarg);
        break;
      default:
        fprintf(stderr, "Invalid command line option %s.\n", optarg);
        break;
    }
  }
}


//------------------------------------------------------------------------------
/**
 * \brief Return column permutation specification for SuperLU (seq and mt).
 *
 * The user can specify the column permutation to use via environment variable
 * ORDERING, e.g., ORDERING=NATURAL. Thus, we can test the  different options
 * without recompiling ;-)
 *
 * Works for SuperLU and SuperLU_MT since the enums colperm_t of both library
 * versions have identical definitions.
 */
inline colperm_t getSuperLUOrdering()
{
  colperm_t ret = COLAMD; // Default value.
  printf("SuperLU Ordering: ");
  const char* env_p = getenv("ORDERING");
  // if(const char* env_p = getenv("ORDERING")) {
  if (env_p != NULL) {
    // std::string env(env_p);
    //if (env.compare("NATURAL") == 0) {
    if (strcmp(env_p, "NATURAL") == 0) {
      printf("NATURAL\n");
      return NATURAL;
    }
    if (strcmp(env_p, "MMD_ATA") == 0) {
      printf("MMD_ATA\n");
      return MMD_ATA;
    }
    if (strcmp(env_p,"MMD_AT_PLUS_A") == 0) {
      printf("MMD_AT_PLUS_A\n");
      return MMD_AT_PLUS_A;
    }
    if (strcmp(env_p,"COLAMD") == 0) {
      printf("COLAMD\n");
      return COLAMD;
    }
    if (strcmp(env_p,"METIS_AT_PLUS_A") == 0) {
      printf("METIS_AT_PLUS_A\n");
      return METIS_AT_PLUS_A;
    }
    if (strcmp(env_p,"PARMETIS") == 0) {
      printf("PARMETIS\n");
      return PARMETIS;
    }
    // ZOLTAN is only defined in SuperLU, not in SuperLU_MT. But since the user
    // guide does not hold any information about ZOLTAN, we do not use it.
//    if (strcmp(env_p,"ZOLTAN") == 0) {
//      //printf( "ZOLTAN\n";
//      m_warning(E_NULL, "Hqp_IpMatrix::getPermcSpec ZOLTAN is not enabled. "
//                        "Set to COLAMD.");
//    } else if (env.compare("MY_PERMC") == 0) {
//      //printf( "MY_PERMC\n";
//      ret = MY_PERMC;
//    }
    printf("Info: Not a valid Ordering. Use COLAMD.\n");
  }
  printf("COLAMD\n");
  return ret;
}

