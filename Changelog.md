## SeDuMi 1.3.7 (2023-03-18)

> https://github.com/sqlp/sedumi/compare/v1.3.6...v1.3.7

- Bugs fixed:
  - [#79: sedumi.m: use Matlab line continuation](https://github.com/sqlp/sedumi/issues/79)

## SeDuMi 1.3.6 (2023-03-17)

> https://github.com/sqlp/sedumi/compare/v1.3.5...v1.3.6

- Bugs fixed:
  - [#73: fix indentation to get rid of compilation warnings](https://github.com/sqlp/sedumi/issues/73)
  - [#74: fix complex interface](https://github.com/sqlp/sedumi/issues/74)
  - [#75: remove eigs warning for function handle matrices](https://github.com/sqlp/sedumi/issues/75)
  - [#76: install_sedumi.m: remove try-catch block](https://github.com/sqlp/sedumi/issues/76)
  - [#77: sedumi_version.m: make SeDuMi version detectable](https://github.com/sqlp/sedumi/issues/77)
  - [#78: add unit test for complex numbers](https://github.com/sqlp/sedumi/issues/78)

Many thanks to @araujoms for his major contributions to this release.

## SeDuMi 1.3.5 (2021-06-15)

> https://github.com/sqlp/sedumi/compare/v1.3.4...v1.3.5

- Merge of "MEX-free implementation of pretransfo and posttransfo" branch by
  Michael C. Grant.
- Improved `install_sedumi.m` function.
- Use GitHub Actions for CI tests.
- Bugs fixed:
  - [#64: function handle problem in minpsdeig.m](https://github.com/sqlp/sedumi/issues/64)
  - [#18: Typo in eigK.m](https://github.com/sqlp/sedumi/issues/18)
  - [#15: results are not apparently deterministic?](https://github.com/sqlp/sedumi/issues/15)

## SeDuMi 1.3.4 (2020-01-09)

> https://github.com/sqlp/sedumi/compare/v1.3.2...v1.3.4
>
> Note: version 1.3.3 was not officially released.

- Many many bug fixes in seven years.


## SeDuMi 1.3.2 (2013-07-24, Michael Grant)

> https://github.com/sqlp/sedumi/compare/v1.3.1...v1.3.2

- Merged changes from CVX and SeDuMi/Octave into the distribution.
- Further changes to suppress Matlab Code Analyzer warnings.
- Removed `matrix.h` dependency from `givens.h` for Octave compatibility;
  fixed assertion problem with `givensrot.c`.


## SeDuMi 1.3.1 (2013-04-04)

> https://github.com/sqlp/sedumi/compare/v1.3.0...v1.3.1

- Corrected bugs in `sedumi.m` and `getada.m`.
- Updated .m files as per mlint suggestions.
- Commented `length(b)` vs `length(c)` check in pretransfo.
- Recompiled windows mex files, updated `install_sedumi.m` for MATLAB v8


## New features in SeDuMi 1.1 since 1.05R5

- 2010-02-15 Fix: global variables are used for some large data structures
  (ADA, At).

- 2010-02-15 Fix: `getada1` and `getada2` have been backported to Matlab,
  which actually increases speed, due to parallel BLAS.

- 2009-11-25 Fix: computation in `getsymbada.m` has been reorganized for speed.

- 2009-11-15 Fix: the `psdeig`, `psdinvscale`, `psdjmulv`, `psdfactor`,
  `psdscale` and `triumtriu` mex files were replaced with a pure Matlab
  equivalent.

- 2009-09-05 Fix: SeDuMi now uses wall clock time instead of cputime, as the
  CPU time is misleading in parallel processing environments.
  It is reported in `info.walltime`. `info.cpusec` contains the total cpu time,
  which may be more than the wall time.

- 2009-08-05 Fix: Based on feedback from Didier Henrion and others, the rank
  and infeasibility detection code has been modified.

- 2009-05-10 Fix: Michael C. Grant modified some files to compile under 64-bit
  versions.

- 2009-05-10 Fix: A new installer script is provided by Michael C. Grant.

- 2009-05-10 Fix: Some BLAS calls are fixed in `blkaux.c`.

- 2008-04-09 Fix: BLAS is used to improve performance.

- 2008-04-09 Fix: SeDuMi now works on 64-bit operating systems.

- 2006-10-14 Fix: A bugfix was contributed by Michael C. Grant regarding the
  handling of 2-dimensional Lorentz cones.

- 2006-10-10 Fix: A new, completely rewritten SDPA format reader is provided in
  the conversion folder.  It is fast and it follows the SDPA standard.

- 2006-10-10 Fix: Michael Grant corrected a bug in the postprocessing phase.

- 2006-10-10 Fix: The semidefinite preprocessing assumed that there are no
  free variables.

- 2006-10-10 Fix: Matlab R2006b does not support `fprintf` with `fid=0`,
  (i.e., no output).  A workaround is provided.  Thanks to Johan Löfberg.

- 2006-10-10 Fix: The solutions were incorrect if the error measures were turned on.

- 2006-07-12 Fix: A bug was discovered by Paul Tseng when using free variables
  together with rotated Lorentz cones.

- 2005-06-22 Fix: A bug about free variable handling has been reported by
  Johan Löfberg, the fixed version is now posted.

- 2005-06-21 Fix: A bug about complex variables has been reported by
  Michael C. Grant, the correct version is now available.

- 2005-06-10 Fix: `&&` and `||` were replaced by `&` and `|` to make the
  package compatible with Matlab 6.0 (R12).

- 2005-05-10 Feature: The binaries are now built from Matlab, see `Install.txt`
  for details.

- 2005-03-02 Fix: The `'See also'` references at the end of the help portions
  are now clickable.

- 2005-03-02 Fix: If ADA is actually dense then the symbolic Cholesky
  factorization and the minimum degree ordering are not performed.
  Also in this case `SYMBADA` is created directly as a fully dense matrix of
  ones stored as sparse.

- 2005-03-02 Feature: Included a simple heuristic to turn step-differentiation on/off.

- 2005-03-02 Feature: New options to control preprocessing: `pars.free`, `pars.sdp`.

- 2005-03-02 Feature: Preprocessing routines are now included.
  Diagonal SDP blocks are converted into nonnegative variables.
  Free variables are handled in a quadratic cone,
  split free variables are detected in the linear part.
  The sparsity structure of A is corrected at the end of preprocessing.

- 2005-03-02 Fix: New default values.
  `pars.eps = 1e-8`, `pars.stepdif = 2`, `pars.cg.qprec = 1`.

- 2005-02-02 Feature: If `pars.errors=1` then SeDuMi outputs the six errors
  measures defined in the 7th DIMACS challenge paper that is available at
  <http://plato.asu.edu/dimacs/node3.html>.
  These are in `info.err1`, ..., `info.err6`.

- 2005-01-20 Fix: Minor speed improvement in the `sdmaux*` files (loop unrolling).

- 2005-01-10 Fix: Minor speed improvement in `popK.m`.

- 2005-01-10 Fix: `blkchol` is now invoked directly from `SeDuMi.m` without
  `sparchol`.

- 2005-01-10 Fix: A small change in `Amul.m` resulted in drastic speed
  improvement for middle-sized problems.
  Additionally, `sum(c.*x)` is much faster than `c'*x` if `c` and `x` are sparse.
  The moral is that taking the transpose of a sparse matrix is slow.

- 2004-01-03 Fix: Fixed a small bug in `wregion.m` concerning division by zero.
  (Johan Löfberg)

- 2004-11-29 Fix: Many if-then statements checking for data validity in the C
  code were replaced by `mxAssert` statements.
  The only exceptions are memory failure checks.

- 2004-11-29 Feature: Timing is changed.
  So far SeDuMi reported the CPU time spent in the main IPM loop.
  Now `info.timing` is a vector containing the time spent in the preprocessing,
  IPM and postprocessing, respectively.
  The old `info.cpusec` contains the total time of the algorithm.
  This does not make much difference now, since 99% of the time is spent in
  the IPM loop, but if more stuff is added to the pre- and postprocessing,
  this can be an issue.

- 2004-11-11 Fix: Cleaned up the code with M-Lint.
  All unused variables and never used assignments are removed.
  All `|` and `&` are replaced with `||` and `&&` whenever it was possible.
  This speeds up conditional statements in the preprocessing phase.
  Unused .m and .c files are removed from the distribution.

- 2004-11-03 Feature: pars.stopat can now be a vector.
  The algorithm enters debug mode at each iteration specified in `pars.stopat`.
  The order does not matter, invalid values are ignored.
