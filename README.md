# SuperLUbench
Benchmark and compare SuperLU, SuperLU_MT and SuperLU_Dist


## Read in CSR Format
Create `supermatrix` in compressed column format:
```
dCreate_CompRow_Matrix(&A, _dim, _dim, _nnz, _v->ve, _jv->ive, _iv->ive, SLU_NR, SLU_D, SLU_GE);
```

The linear systems from HQP are in CSR format. *Thus, provide an reading function.*



readMM --> nzval, int_t **rowind, int_t **colptr



## Read Matrix Market Format
SuperLU only provides methods to read in matrices in HB and Triplet format.

The source code directory of SuperLU_MT provides the file dreadmt.c. But it is not included in the packaing library of (at least my) Linux distro. It seems to be unsupported.
This method reads the matrix from stdin. Argh.

SuperLU_DIST provides a method to read in matrices in Matrix Market Format. This seems to be supported and maintained. Can it be used for the other two flavours?



## Iterative Solving Ax=b
For benchmarking reasons, it is necessary to solve Ax=b multiple times.

-R NUM
TODO More Information

### Additional Information for Developers
The used functions `dgssvx()`, `pdgssvx()` and `pdgssvx_ABglobal()` within SuperLUbench may overwrite the structures holding A and b. Thus, the values are stored in additional arrays and restored before iteratively calling the solve functions.

In SuperLU and SuperLU_MT it seems to be sufficient to recreate the arrays `double *a` and `double *rhsb`.

Additional to this two arrays, the array `asub` (holding the row indices) needs to be recreated in SuperLU_DIST since `pdgssvx_ABglobal()` may overwrite it. I'm not sure, if this also holds for the methods `dgssvx()` and `pdgssvx()`, respectively.