# SeDuMi: Optimization over symmetric cones

#### [Click here](https://github.com/sqlp/sedumi/releases/latest) to download the latest SeDuMi bundle. These bundles now include pre-compiled MATLAB and Octave MEX files files for Windows, Linux, and macOS (Intel and Apple Silicon).

**SeDuMi (Self-Dual-Minimization)** is a Matlab/GNU Octave package for solving
convex optimization problems involving linear equations and inequalities,
second-order cone constraints, and semidefinite constraints (linear matrix
inequalities).

Please note that this is an *unofficial* repository for SeDuMi.
The [official SeDuMi site](https://sedumi.ie.lehigh.edu/) is hosted by the
[CORAL Lab](https://coral.ise.lehigh.edu/) at the
[Department of Industrial and Systems Engineering](https://engineering.lehigh.edu/ise)
at [Lehigh University](https://www.lehigh.edu/).
This repository was originally not intended to remain a permanent fork
from the last official SeDuMi release.
However, this repository contains a still maintained versions of SeDuMi.

This version of SeDuMi is distributed under the
[GNU General Public License 2.0](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).


## Installation instructions

Download and add the extracted SeDuMi directory to the Matlab load path
and run

    install_sedumi

This version of SeDuMi comes with pre-compiled binaries for Matlab (MS Windows,
Linux, and macOS).  In case the binaries don't work for your system,
try to run

    install_sedumi -rebuild

Note that you may have to install a compiler for Matlab first.
Using GNU Octave, the binaries are compiled automatically.


## Development History

The original version was developed by Jos. F. Sturm, who sadly passed away
in 2003.  Development continued under the direction of
[Prof. Tamás Terlaky](https://coral.ise.lehigh.edu/terlaky/) by former
Ph.D. students Imre Pólik and Oleksandr Romanko.

Of course, like many open-source projects, it has benefited considerably from
contributions by others, including the developers of
[YALMIP](https://yalmip.github.io/) and [CVX](http://cvxr.com),
two modeling frameworks for optimization that use SeDuMi as a solver.
The authors of these packages co-administer this repo,
along with [Jonathan Currie](http://www.i2c2.aut.ac.nz/) from AUT University.

For a list of changes for each SeDuMi version, check the file
[Changelog.md](https://github.com/sqlp/sedumi/blob/master/Changelog.md).


## Citation

If you find this software useful, please cite it in your publication as follows:

```
@article{doi:10.1080/10556789908805766,
  author = {Sturm, Jos F.},
  title = {Using SeDuMi 1.02, A {MATLAB} toolbox for optimization over symmetric cones},
  journal = {Optimization Methods and Software},
  volume = {11},
  number = {1-4},
  pages = {625-653},
  year  = {1999},
  doi = {10.1080/10556789908805766},
  URL = {https://doi.org/10.1080/10556789908805766}
}
```

Account of original sources:

The files `ordmmd.c` and `symfct.c` are C-versions of Fortran modules by
Esmond G. Ng and Barry W. Peyton, Oak Ridge National Laboratory, 1994,
from which large parts are based on SPARSPAK-A RELEASE III by Joseph W.H. Liu,
University of Waterloo, 1984.

All other files in this distribution are by Jos F. Sturm, Tilburg University,
2001, updated and further developed by Imre Polik, McMaster University,
Hamilton, Canada.


## Report problems

You are welcome to submit bug reports or request for help on the
[GitHub issue page](https://github.com/sqlp/sedumi/issues).
We cannot guarantee that they will be addressed in a timely fashion,
we will do our best.

### Development notes

As of April 2024, this repository uses
[GitHub Actions](https://github.com/features/actions) to compile
MEX files for Linux, Windows, and macOS (both Intel and Apple
Silicon variants). Whenever a new Git tag is pushed to the
repository, these actions automatically create `.zip` and `.tgz`
bundles of that version of the code, including those compiled
MEX files, and publishes those bundles to the
[Releases](https://github.com/sqlp/sedumi/releases) page.

If you wish to contribute fixes or improvements to this repository, please feel free to submit a
[pull request](https://github.com/sqlp/sedumi/pulls).
