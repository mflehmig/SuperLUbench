/*! \file
 Copyright (c) 2003, The Regents of the University of California, through
 Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from U.S. Dept. of Energy)

 All rights reserved.

 The source code is distributed under BSD license, see the file License.txt
 at the top-level directory.
 */

/*! @file 
 * \brief Driver program for pdgssvx_ABglobal example
 *
 * <pre>
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 * </pre>
 */

#include <unistd.h>
#include <superlu_ddefs.h>
#include "../Util.h"
#include "../dreadMM.h"

void parse_command_line(int argc, char *argv[], int* nprow, int* npcol,
                        char *fileA, char *fileB, char *fileX, int *reps);
colperm_t getSuperLUOrdering();

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * The driver program pddrive_ABglobal.
 *
 * This example illustrates how to use pdgssvx_ABglobal with the full
 * (default) options to solve a linear system.
 * 
 * Five basic steps are required:
 *   1. Initialize the MPI environment and the SuperLU process grid
 *   2. Set up the input matrix and the right-hand side
 *   3. Set the options argument
 *   4. Call pdgssvx_ABglobal
 *   5. Release the process grid and terminate the MPI environment
 *
 * On an IBM SP, the program may be run by typing
 *    poe pddrive_ABglobal -r <proc rows> -c <proc columns> <input_file> -procs <p>
 * </pre>
 */
int main(int argc, char *argv[])
{
  superlu_dist_options_t options;
  SuperLUStat_t stat;
  SuperMatrix A;
  ScalePermstruct_t ScalePermstruct;
  LUstruct_t LUstruct;
  gridinfo_t grid;
  double *berr;
  double *a_0, *a, *b_0, *b, *xact;
  int_t *asub_0, *asub, *xa;
  int_t m, n, nnz;
  int_t nprow, npcol;
  int iam, info, ldb, ldx, nrhs;
  char trans[1];
  FILE *fp;
  char fileA[128];
  char fileB[128];
  char fileX[128];
  int reps = 1;

  nprow = 1; /* Default process rows.      */
  npcol = 1; /* Default process columns.   */
  nrhs = 1; /* Number of right-hand side. */

  /* ------------------------------------------------------------
   INITIALIZE MPI ENVIRONMENT.
   ------------------------------------------------------------*/
  MPI_Init(&argc, &argv);

  /* Parse command line argv[]. */
  parse_command_line(argc, argv, &nprow, &npcol, fileA, fileB, fileX, &reps);

  /* ------------------------------------------------------------
   INITIALIZE THE SUPERLU PROCESS GRID.
   ------------------------------------------------------------*/
  superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);

  /* Bail out if I do not belong in the grid. */
  iam = grid.iam;
  if (iam >= nprow * npcol)
    goto out;

  /* ------------------------------------------------------------
   PROCESS 0 READS THE MATRIX A, AND THEN BROADCASTS IT TO ALL
   THE OTHER PROCESSES.
   ------------------------------------------------------------*/
  if (!iam)
  {
    /* Read the matrix stored on disk in Matrix-Market format. */
    fp = fopen(fileA, "r");
    dreadMM_dist(fp, &m, &n, &nnz, &a_0, &asub_0, &xa);
    fclose(fp);

    printf("Input matrix file: %s\n", fileA);
    printf("\tDimension\t" IFMT "x" IFMT "\t # nonzeros " IFMT "\n", m, n, nnz);
    printf("\tProcess grid\t%d X %d\n", (int) grid.nprow, (int) grid.npcol);

    /* Broadcast matrix A to the other PEs. */
    MPI_Bcast(&m, 1, mpi_int_t, 0, grid.comm);
    MPI_Bcast(&n, 1, mpi_int_t, 0, grid.comm);
    MPI_Bcast(&nnz, 1, mpi_int_t, 0, grid.comm);
    MPI_Bcast(a_0, nnz, MPI_DOUBLE, 0, grid.comm);
    MPI_Bcast(asub_0, nnz, mpi_int_t, 0, grid.comm);
    MPI_Bcast(xa, n + 1, mpi_int_t, 0, grid.comm);
  }
  else
  {
    // Receive matrix A from PE 0.
    MPI_Bcast(&m, 1, mpi_int_t, 0, grid.comm);
    MPI_Bcast(&n, 1, mpi_int_t, 0, grid.comm);
    MPI_Bcast(&nnz, 1, mpi_int_t, 0, grid.comm);

    // Allocate storage for compressed column representation.
    dallocateA_dist(n, nnz, &a_0, &asub_0, &xa);

    MPI_Bcast(a_0, nnz, MPI_DOUBLE, 0, grid.comm);
    MPI_Bcast(asub_0, nnz, mpi_int_t, 0, grid.comm);
    MPI_Bcast(xa, n + 1, mpi_int_t, 0, grid.comm);
  }

  // Allocate memory for a. Memory for a_O is already allocated.
  if (!(a = doubleMalloc_dist(nnz)))
    ABORT("Malloc fails for a[]");
  memcpy(a, a_0, sizeof(double) * nnz);
  if (!(asub = intMalloc_dist(nnz)))
    ABORT("Malloc fails for asub[]");
  memcpy(asub, asub_0, sizeof(int) * nnz);

  // Create compressed column matrix for A.
  dCreate_CompCol_Matrix_dist(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D,
                              SLU_GE);

  // Allocate memory for b and x.
  if (!(b_0 = doubleMalloc_dist(m)))
    ABORT("Malloc fails for b_0[]");
  if (!(b = doubleMalloc_dist(m)))
    ABORT("Malloc fails for b[]");
  if (!(xact = doubleMalloc_dist(n)))
    ABORT("Malloc fails for xact[]");

  // Read matrix b and x from a file in Matrix Market format.
  if (strcmp(fileB, "") != 0 && strcmp(fileX, "") != 0)
  {
    // Only masters reads the files. Than broadcasts to the other PEs.
    if (!iam)
    {
      printf("Read in b and x from provided files.\n");
      fp = fopen(fileB, "r");
      if (fp == NULL)
      {
        perror("input error: could not open file containing b\n");
        exit(-2);
      }
      readVec(n, b_0, fp);
      fclose(fp);

      fp = fopen(fileX, "r");
      if (fp == NULL)
      {
        perror("input error: could not open file containing x\n");
        exit(-2);
      }
      readVec(n, xact, fp);
      fclose(fp);

      MPI_Bcast(b_0, m, MPI_DOUBLE, 0, grid.comm);
      MPI_Bcast(xact, m, MPI_DOUBLE, 0, grid.comm);
    }
    else
    {
      MPI_Bcast(b_0, m, MPI_DOUBLE, 0, grid.comm);
      MPI_Bcast(xact, m, MPI_DOUBLE, 0, grid.comm);
    }
  } // Read b and x

//  *trans = 'N';
//  ldx = n;
//  ldb = m;
//  dGenXtrue_dist(n, nrhs, xact, n);
//  dFillRHS_dist(trans, nrhs, xact, n, &A, b, m);

  if (!(berr = doubleMalloc_dist(nrhs)))
    ABORT("Malloc fails for berr[].");

  /* ------------------------------------------------------------
   NOW WE SOLVE THE LINEAR SYSTEM.
   ------------------------------------------------------------*/

  /* Set the default input options:
   options.Fact = DOFACT;
   options.Equil = YES;
   options.ColPerm = METIS_AT_PLUS_A;
   options.RowPerm = LargeDiag;
   options.ReplaceTinyPivot = YES;
   options.Trans = NOTRANS;
   options.IterRefine = DOUBLE;
   options.SolveInitialized = NO;
   options.RefineInitialized = NO;
   options.PrintStat = YES;
   */
  set_default_options_dist(&options);
  options.ColPerm = getSuperLUOrdering();

  if (!iam)
  {
    print_sp_ienv_dist(&options);
    print_options_dist(&options);
  }

  /* Initialize ScalePermstruct and LUstruct. */
  ScalePermstructInit(m, n, &ScalePermstruct);
  LUstructInit(n, &LUstruct);
  // Solve the system 'reps' times.  Since A, b and x are overwritten by pdgssvx_ABglobal(),
  // we need to recreate them in each iteration (using the same memory locations).
  int k;
  long long start, finish;
  for (k = 0; k < reps; ++k)
  {
    printf("\n\n-- Iteration # %i\n", k);

    // Fresh copy of a (and thus A), because values in a may be overwritten by pdgssvx_ABglobal().
    memcpy(a, a_0, sizeof(double) * nnz);
    // Fresh copy of asub (and thus A), because values in asub may be overwritten by pdgssvx_ABglobal().
    memcpy(asub, asub_0, sizeof(int) * nnz);
    // Fresh copy of rhs b, because values in rhsb may be overwritten by pdgssvx_ABglobal().
    memcpy(b, b_0, sizeof(double) * m);

    /* Initialize the statistics variables. */
    PStatInit(&stat);

    start = current_timestamp();
    /* Call the linear equation solver. */
    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b, m, 1, &grid,
                     &LUstruct, berr, &stat, &info);
    finish = current_timestamp();

    /* Compare value of info with A->ncol (=m). */
    if (!iam)
    {
      if (info == 0)
      {
        printf("SUCCESS: pdgssvx_ABglobal() returns info %d which means FINE\n", info);
        printf("Time for pdgssvx_ABglobal() [1e-6 s]: %lli\n", (finish - start));
      }
      if (info > 0 && info <= m)
        printf(
            "FAILED: pdgssvx() returns info %d which means U(%d,%d) is exactly zero\n",
            info, info, info);
      if (info > 0 && info > m)
        printf(
            "FAILED: pdgssvx() returns info %d which means memory allocation failed\n",
            info);
    }

    /* Check the accuracy of the solution (only if pdgssvx_ABglobal succeeded). */
    if (!iam && info == 0)
      dinf_norm_error_dist(n, nrhs, b, m, xact, n, &grid);

    PStatPrint(&options, &stat, &grid); /* Print the statistics. */

  } // reps

  /* ------------------------------------------------------------
   DEALLOCATE STORAGE.
   ------------------------------------------------------------*/
  PStatFree(&stat);
  Destroy_CompCol_Matrix_dist(&A);
  Destroy_LU(n, &grid, &LUstruct);
  ScalePermstructFree(&ScalePermstruct);
  LUstructFree(&LUstruct);
  SUPERLU_FREE(b);
  SUPERLU_FREE(xact);
  SUPERLU_FREE(berr);

  /* ------------------------------------------------------------
   RELEASE THE SUPERLU PROCESS GRID.
   ------------------------------------------------------------*/
  out: superlu_gridexit(&grid);

  /* ------------------------------------------------------------
   TERMINATES THE MPI EXECUTION ENVIRONMENT.
   ------------------------------------------------------------*/
  MPI_Finalize();
}

/*
 * Parse command line options.
 */
void parse_command_line(int argc, char *argv[], int* nprow, int* npcol,
                        char *fileA, char *fileB, char *fileX, int *reps)
{
  int c;
  extern char *optarg;

  while ((c = getopt(argc, argv, "hr:c:A:b:x:R:")) != EOF)
  {
    switch (c)
    {
      case 'h':
        printf("Options:\n");
        printf("\t-r <int>: process rows    (default %4d)\n", *nprow);
        printf("\t-c <int>: process columns (default %4d)\n", *npcol);
        exit(0);
        break;
      case 'r':
        *nprow = atoi(optarg);
        break;
      case 'c':
        *npcol = atoi(optarg);
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
 * \brief Return column permutation specification for SuperLU_DIST.
 *
 * The user can specify the column permutation to use via environment variable
 * ORDERING, e.g., ORDERING=NATURAL. Thus, we can test the  different options
 * without recompiling ;-)
 *
 */
inline colperm_t getSuperLUOrdering()
{
  colperm_t ret = NATURAL;  // Default value.
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
  printf("NATURAL\n");
  return ret;
}

