/* ************************************************************
   MODULE sdmaux*.c  -- Several low-level subroutines for the
   mex-files in the Self-Dual-Minimization package.

% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
%
% Copyright (C) 2001 Jos F. Sturm (up to 1.05R5)
%   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% Affiliation SeDuMi 1.03 and 1.04Beta (2000):
%   Dept. Quantitative Economics, Maastricht University, the Netherlands.
%
% Affiliations up to SeDuMi 1.02 (AUG1998):
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA

 ************************************************************ */

#include "blksdp.h"
/* ************************************************************
   TIME-CRITICAL PROCEDURE -- scalarmul
   Computes  r = alpha * x  using BLAS.
   ************************************************************ */
void scalarmul(double *r, double alpha,const double *x,mwIndex n)
{
    blasint one=1,nn=n;
    FORT(dcopy)(&nn,(double*)x,&one,r,&one);
    FORT(dscal)(&nn,&alpha,r,&one);
}

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- addscalarmul
   Computes  r += alpha * x  using BLAS.
   ************************************************************ */
void addscalarmul(double *r, double alpha,const double *x,mwIndex n)
{
    blasint one=1,nn=n;
    FORT(daxpy)(&nn,&alpha,(double*)x,&one,r,&one);
}

/* ************************************************************
   TIME-CRITICAL PROCEDURE -- subscalarmul(x,alpha,y,n)
   Computes x -= alpha * y using BLAS.
   ************************************************************ */
void subscalarmul(double *x, double alpha, const double *y, mwIndex n)
{
    blasint one=1,nn=n;
    double minusalpha=-alpha;
    FORT(daxpy)(&nn,&minusalpha,(double*)y,&one,x,&one);
}
