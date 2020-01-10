In this folder there are some example problems that you can try with SeDuMi.
For more problems see the Homepage of Hans D. Mittelmann
http://plato.asu.edu/sub/testcases.html.

- **arch0.mat**: A middle-sized SDP problem from the SDPLIB set maintained by
  Biran Borchers.  This particular problem is a truss topology design problem,
  contributed by Katsuki Fujisawa.  For details see

  > T. Nakamura and M. Ohsaki. A Natural Generator of Optimum Topology of
  > Plane Trusses for Specified Fundamental Frequency.  Computer Methods in
  > Applied Mechanics and Engineering 94(1992):113-129.
  > DOI: [10.1016/0045-7825(92)90159-H](https://doi.org/10.1016/0045-7825(92)90159-H).

  Optimal value is `5.66517e-01`.


- **control07.mat**: Again from SDPLIB, a problem from control and system theory
  contributed by Katsuki Fujisawa.  For details see

  > K. Fujisawa, M. Fukuda, M. Kojima, and K. Nakata.  Numerical Evaluation
  > of SDPA (Semidefinite Programming Algorithm)  Technical Report B-330,
  > Department of Mathematical and Computing Sciences, Tokyo Institute of
  > Technology, September, 1997.

  Optimal value is `2.06251e+01`.


- **nb.mat**:  This problem is by Robert Vanderbei, it was used at the
  [7th DIMACS computational Challenge](http://archive.dimacs.rutgers.edu/Challenges/Seventh/).
  It is a middle-sized mixed second-order/linear problem.

  Optimal value is `-0.05070309`.


- **OH_2Pi_STO-6GN9r12g1T2.mat**: This is a middle-sized SDP problem from
  electronic structure calculation.  For details see

  > Z. Zhao, B. J. Braams, M. Fukuda, M. L. Overton, and J. K. Percus,
  > "The reduced density matrix method for electronic structure calculations
  > and the role of three-index representability", October, 2003.
  > DOI: [10.1063/1.1636721](https://doi.org/10.1063/1.1636721).

  at http://mf.c.titech.ac.jp/mituhiro/software.html.


- **trto3.mat**: A problem by Kocvara, from single-load truss topology design.
  Normally formulated as LP, here reformulated as SDP for testing purposes.
  For details see

  > A. Ben-Tal and A. Nemirovski. Lectures on Modern Convex Optimization.
  > MPS-SIAM Series on Optimization. SIAM Philadelphia, 2001.
  > DOI: [10.1137/1.9780898718829](https://doi.org/10.1137/1.9780898718829).

  > M. Kocvara and J. Zowe. How mathematics can help in design of mechanical
  > structures. In D.F. Griffiths and G.A. Watson, eds., Numerical Analysis
  > 1995, Longman, Harlow, 1996, pp. 76--93.
