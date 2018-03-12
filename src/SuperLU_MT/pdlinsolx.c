/*! \file
 Copyright (c) 2003, The Regents of the University of California, through
 Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from U.S. Dept. of Energy)

 All rights reserved.

 The source code is distributed under BSD license, see the file License.txt
 at the top-level directory.
 */

/*
 * -- SuperLU MT routine (version 3.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <slu_mt_ddefs.h>
#include "../Util.h"
#include "../dreadMM.h"

void parse_command_line(int argc, char *argv[], int_t *nprocs, int_t *lwork,
                        int_t *w, int_t *relax, double *u, fact_t *fact,
                        trans_t *trans, equed_t *equed, char *fileA,
                        char *fileB, char *fileX, int *reps);
colperm_t getSuperLUOrdering();
int getNumThreads();

int main(int argc, char *argv[])
{
  SuperMatrix A, L, U;
  SuperMatrix B, X;
  NCformat *Astore;
  SCPformat *Lstore;
  NCPformat *Ustore;
  int_t nprocs;
  fact_t fact;
  trans_t trans;
  yes_no_t usepr;
  equed_t equed;
  double *a_0, *a; /* Nonzero values of A. Since dgssvx may overwrite the
   values in A and thus a, we save the original values in a.*/
  int_t *asub, *xa;
  int_t *perm_c; /* column permutation vector */
  int_t *perm_r; /* row permutations from partial pivoting */
  void *work;
  superlumt_options_t superlumt_options;
  int_t info, lwork, nrhs, ldx, panel_size, relax;
  int_t m, n, nnz;
  double *rhsb_0, *rhsb, *rhsx, *xact;
  double *R, *C;
  double *ferr, *berr;
  double u, drop_tol, rpg, rcond;
  superlu_memusage_t slu_mem;
  char fileA[128];
  char fileB[128] = "";
  char fileX[128] = "";
  colperm_t permc_spec;
  int reps;
  char tmpfile[] = ".infnorm.txt";

  /* Default parameters to control factorization. */
  reps = 1;
  nprocs = getNumThreads();
  fact = EQUILIBRATE;
  trans = NOTRANS;
  equed = NOEQUIL;
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  u = 1.0;
  usepr = NO;
  drop_tol = 0.0;
  lwork = 0;
  nrhs = 1;

  permc_spec = getSuperLUOrdering();

  /* Command line options to modify default behaviour. */
  parse_command_line(argc, argv, &nprocs, &lwork, &panel_size, &relax, &u,
                     &fact, &trans, &equed, fileA, fileB, fileX, &reps);

  if (lwork > 0)
  {
    work = SUPERLU_MALLOC(lwork);
    printf("Use work space of size LWORK = " IFMT " bytes\n", lwork);
    if (!work)
      SUPERLU_ABORT("PDLINSOLX: cannot allocate work[]");
  }

  /* Read matrix A from a file in Matrix Market format.*/
  FILE *fp;
  fp = fopen(fileA, "r");
  if (fp == NULL)
  {
    perror("input error: could not open file\n");
    exit(-2);
  }
  dreadMM(fp, &m, &n, &nnz, &a_0, &asub, &xa);
  fclose(fp);

  // Create A
  if (!(a = doubleMalloc(nnz)))
    SUPERLU_ABORT("Malloc fails for a[].");
  memcpy(a, a_0, sizeof(double) * nnz);
  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
  Astore = A.Store;
  printf("Dimension " IFMT "x" IFMT "; # nonzeros " IFMT "\n", A.nrow, A.ncol,
         Astore->nnz);

  if (!(rhsb_0 = doubleMalloc(m * nrhs)))
    SUPERLU_ABORT("Malloc fails for rhsb_0[].");
  if (!(rhsb = doubleMalloc(m * nrhs)))
    SUPERLU_ABORT("Malloc fails for rhsb[].");
  if (!(rhsx = doubleMalloc(m * nrhs)))
    SUPERLU_ABORT("Malloc fails for rhsx[].");
  dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
  dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);

  xact = doubleMalloc(n * nrhs);
  ldx = n;

  /* Read vectors b and x from a file in Matrix Market format.*/
  if (strcmp(fileB, "") != 0 && strcmp(fileX, "") != 0)
  {
    printf("Read in b and x from provided files.\n");
    fp = fopen(fileB, "r");
    if (fp == NULL)
    {
      perror("input error: could not open file containing b\n");
      exit(-2);
    }
    readVec(n, rhsb_0, fp);
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

  if (!(perm_r = intMalloc(m)))
    SUPERLU_ABORT("Malloc fails for perm_r[].");
  if (!(perm_c = intMalloc(n)))
    SUPERLU_ABORT("Malloc fails for perm_c[].");
  if (!(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))))
    SUPERLU_ABORT("SUPERLU_MALLOC fails for R[].");
  if (!(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))))
    SUPERLU_ABORT("SUPERLU_MALLOC fails for C[].");
  if (!(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))))
    SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr[].");
  if (!(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))))
    SUPERLU_ABORT("SUPERLU_MALLOC fails for berr[].");

  superlumt_options.nprocs = nprocs;
  superlumt_options.fact = fact;
  superlumt_options.trans = trans;
  superlumt_options.refact = NO;  //refact;
  superlumt_options.panel_size = panel_size;
  superlumt_options.relax = relax;
  superlumt_options.usepr = usepr;
  superlumt_options.drop_tol = drop_tol;
  superlumt_options.diag_pivot_thresh = u;
  superlumt_options.SymmetricMode = NO;
  superlumt_options.PrintStat = YES;
  superlumt_options.perm_c = perm_c;
  superlumt_options.perm_r = perm_r;
  superlumt_options.work = work;
  superlumt_options.lwork = lwork;
  if (!(superlumt_options.etree = intMalloc(n)))
    SUPERLU_ABORT("Malloc fails for etree[].");
  if (!(superlumt_options.colcnt_h = intMalloc(n)))
    SUPERLU_ABORT("Malloc fails for colcnt_h[].");
  if (!(superlumt_options.part_super_h = intMalloc(n)))
    SUPERLU_ABORT("Malloc fails for colcnt_h[].");

  printf(
      "\n#Iter, info, #NNZ in L, #NNZ in U, L\\U [MB], Total [MB], ||X-Xtrue||/||X||, RPG, RCN, dgssvx()\n");
  double err;
  char dummy[128];
  int k;
  long long start, finish;
  // Solve the system 'reps' times.  Since A, b and x are overwritten by dgssvx(),
  // we need to recreate them in each iteration (using the same memory locations).
  for (k = 0; k < reps; ++k)
  {
    //! printf("\n\n-- Iteration # %i\n", k);

    // Fresh copy of a (and thus A), because values in a may be overwritten by dgssvx().
    memcpy(a, a_0, sizeof(double) * nnz);
    // Fresh copy of rhs b, because values in rhsb may be overwritten by dgssvx().
    memcpy(rhsb, rhsb_0, sizeof(double) * m);

    // Solve the system and compute the condition number
    // and error bounds using pdgssvx.
    start = current_timestamp();
    get_perm_c(permc_spec, &A, perm_c);

    pdgssvx(nprocs, &superlumt_options, &A, perm_c, perm_r, &equed, R, C, &L,
            &U, &B, &X, &rpg, &rcond, ferr, berr, &slu_mem, &info);
    finish = current_timestamp();
    //! printf("pdgssvx(): info " IFMT "\n", info);
    //! printf("Time for pdgssvx() [1e-6 s]: %lli\n", (finish - start));

    if (info == 0 || info == n + 1)
    {
//      printf("Recip. pivot growth = %e\n", rpg);
//      printf("Recip. condition number = %e\n", rcond);
//      printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
//      for (i = 0; i < nrhs; ++i)
//        printf(IFMT "%16e%16e\n", i + 1, ferr[i], berr[i]);

      Lstore = (SCPformat *) L.Store;
      Ustore = (NCPformat *) U.Store;
//      printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
//      printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
//      printf("No of nonzeros in L+U = " IFMT "\n",
//             Lstore->nnz + Ustore->nnz - n);
//      printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
//             superlu_memusage.for_lu / 1e6, superlu_memusage.total_needed / 1e6,
//             superlu_memusage.expansions);
//      printf("superlumt_options.nprocs = %d \n", superlumt_options.nprocs);
//      fflush(stdout);

      /* This is how you could access the solution matrix. */
      //    double *sol = (double*) ((DNformat*) X.Store)->nzval;
      //    int ret = compareSolution(sol, m, xact, eps);
      //    if (!ret)
      //      printf("Computed and provided solutions match within %g eps.\n", eps);
      //    else
      //      printf("Computed and provided solutions DO NOT match within %g eps.\n", eps);
//     dinf_norm_error(nrhs, &X, xact);
      // Hack: We want inf norm. But SuperLU method dinf_norm_error() does not
      //       store the norm value into a variable. Instead it calculates the
      //       norm value and prints it to stderr. So, we redirect the stderr
      //       to file and than read the value from file.
      int stdout_fd = dup(STDOUT_FILENO);
      freopen(tmpfile, "w", stdout);
      dinf_norm_error(nrhs, &X, xact);
      fclose(stdout);
      dup2(stdout_fd, STDOUT_FILENO);
      stdout = fdopen(STDOUT_FILENO, "w");
      close(stdout_fd);
      // tmpfile contains "||X - Xtrue||/||X|| = 6.454921e-10", so we read strings
      // until we get "=". Than, we read a double - the norm value!
      fp = fopen(tmpfile, "r");
      do
      {
        fscanf(fp, "%s", dummy);
      }
      while (strcmp(dummy, "="));
      fscanf(fp, "%le", &err);
      fclose(fp);
      // End Hack.

      printf("%d, %i, %i, %d, %d, %.3f, %.3f, %le, %e, %e, %lli \n",
             superlumt_options.nprocs, k, info, Lstore->nnz, Ustore->nnz,
             slu_mem.for_lu / 1e6, slu_mem.total_needed / 1e6, err, rpg, rcond,
             finish - start);
    }
    else if (info > 0 && lwork == -1)
      printf("** Estimated memory: " IFMT " bytes\n", info - n);

  }  // reps
  int ret = remove(tmpfile);
  // Todo: Handle error.

  SUPERLU_FREE(rhsb);
  SUPERLU_FREE(rhsx);
  SUPERLU_FREE(xact);
  SUPERLU_FREE(perm_r);
  SUPERLU_FREE(perm_c);
  SUPERLU_FREE(R);
  SUPERLU_FREE(C);
  SUPERLU_FREE(ferr);
  SUPERLU_FREE(berr);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperMatrix_Store(&X);
  SUPERLU_FREE(superlumt_options.etree);
  SUPERLU_FREE(superlumt_options.colcnt_h);
  SUPERLU_FREE(superlumt_options.part_super_h);
  if (lwork == 0)
  {
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);
  }
  else if (lwork > 0)
    SUPERLU_FREE(work);
}

/*
 * Parse command line options.
 */
void parse_command_line(int argc, char *argv[], int_t *nprocs, int_t *lwork,
                        int_t *w, int_t *relax, double *u, fact_t *fact,
                        trans_t *trans, equed_t *equed, char *fileA,
                        char *fileB, char *fileX, int *reps)
{
  int c;
  extern char *optarg;

  while ((c = getopt(argc, argv, "hp:l:w:s:u:f:t:e:A:b:x:R:")) != EOF)
  {
    switch (c)
    {
      case 'h':
        printf("Options:\n");
        printf("\t-p <int> - number of processes\n");
        printf("\t-l <int> - length of work[*] array\n");
        printf("\t-w <int> - panel size\n");
        printf("\t-s <int> - maximum size of relaxed supernodes\n");
        printf("\t-u <int> - pivoting threshold\n");
        printf("\t-f <FACTORED/DOFACT/EQUILIBRATE> - factor control\n");
        printf("\t-t <NOTRANS/TRANS/CONJ> - transpose or not\n");
        //printf("\t-r <NO/YES> - refactor or not\n");
        printf("\t-e <NOEQUIL/ROW/COL/BOTH> - equilibrate or not\n");
        printf("\t-A <FILE> - File holding matrix A in Matrix Market format\n");
        printf(
            "\t-b <FILE> - File holding right hand side vector b in Matrix Market format\n");
        printf(
            "\t-x <FILE> - File holding known solution vector x in Matrix Market format\n");
        printf("\t-R <NUM> - Number of iteratively solve Ax=b\n");
        exit(1);
        break;
      case 'p':
        *nprocs = atoi(optarg);
        break;
      case 'l':
        *lwork = atoi(optarg);
        break;
      case 'w':
        *w = atoi(optarg);
        break;
      case 's':
        *relax = atoi(optarg);
        break;
      case 'u':
        *u = atof(optarg);
        break;
      case 'f':
        *fact = (fact_t) atoi(optarg);
        break;
      case 't':
        *trans = (trans_t) atoi(optarg);
        break;
//      case 'r':
//        *refact = (yes_no_t) atoi(optarg);
//        break;
      case 'e':
        *equed = (equed_t) atoi(optarg);
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
      case 'R':  // repetitions
        *reps = atoi(optarg);
        break;
      default:
        fprintf(stderr, "Invalid command line option %s.\n", optarg);
        break;
    }
  }
}

/**
 * \brief Return value of env. variable OMP_NUM_THREADS.
 */
int getNumThreads()
{
  int nT;
  const char* env = getenv("OMP_NUM_THREADS");
  // if(const char* env_p = getenv("ORDERING")) {
  if (env != NULL)
  {
    nT = atoi(env);
    if (nT <= 0)
    {
      printf(
          "\nWARNING: Value of OMP_NUM_THREADS is less or equal to zero. Set "
          "number of threads to 1.\n");
      nT = 1;
    }
  }
  else
  {
    printf(
        "\nWARNING: Number of threads is not specified. Set it via environment"
        "variable OMP_NUM_THREADS=NUM. Use default value of 1 thread.\n");
    nT = 1;
  }
  return nT;
}

/**
 * \brief Return column permutation specification for SuperLU_MT.
 *
 * The user can specify the column permutation to use via environment variable
 * ORDERING, e.g., ORDERING=NATURAL. Thus, we can test the  different options
 * without recompiling ;-)
 */
inline colperm_t getSuperLUOrdering()
{
  colperm_t ret = COLAMD;  // Default value.
  printf("SuperLU Ordering: ");
  const char* env_p = getenv("ORDERING");
  // if(const char* env_p = getenv("ORDERING")) {
  if (env_p != NULL)
  {
    // std::string env(env_p);
    //if (env.compare("NATURAL") == 0) {
    if (strcmp(env_p, "NATURAL") == 0)
    {
      printf("NATURAL\n");
      return NATURAL;
    }
    if (strcmp(env_p, "MMD_ATA") == 0)
    {
      printf("MMD_ATA\n");
      return MMD_ATA;
    }
    if (strcmp(env_p, "MMD_AT_PLUS_A") == 0)
    {
      printf("MMD_AT_PLUS_A\n");
      return MMD_AT_PLUS_A;
    }
    if (strcmp(env_p, "COLAMD") == 0)
    {
      printf("COLAMD\n");
      return COLAMD;
    }
    if (strcmp(env_p, "METIS_AT_PLUS_A") == 0)
    {
      printf("METIS_AT_PLUS_A\n");
      return METIS_AT_PLUS_A;
    }
    if (strcmp(env_p, "PARMETIS") == 0)
    {
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

