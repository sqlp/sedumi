function lab = eig(x,K)

% [lab,q,f] = eigK(x,K)
%
% EIGK  Computes the spectral values ("eigenvalues") or even the complete
%       spectral decomposition  of a vector x with respect to a self-dual
%       homogeneous cone K.
%
% > LAB = EIGK(x,K) This yield the spectral values of x with respect to
%       the self-dual homogeneous cone that you describe in the structure
%       K. Up to 4 fields can be used, called K.l, K.q, K.r, and K.s, for
%       Linear, Quadratic, Rotated, and Semi-definite. Type `help sedumi' 
%       for more details on this structure.
%
%       The length of the vector LAB is the order of the cone. Note that
%       x in K if and only if LAB>=0, and x in int(K) if and only if LAB>0.
%
% See also sedumi, mat, vec, eyeK.

% Complete rewrite for SeDuMi 1.3 by Michael C. Grant
% Copyright (C) 2013 Michael C. Grant.
%
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
%

xi = 0;
lab_c = [];
if isfield(K,'f'),
    xi = xi + K.f;
end
if isfield(K,'l') && K.l > 0,
    lab_c{end+1} = x(xi+1:xi+K.l);
    xi = xi+K.l;
end
if isfield(K,'q') && ~isempty(K.q),
    lab = zeros(2*length(K.q),1);
    li = 0;
    for k = 1:length(K.q)
        kk = K.q(k);
        x0 = x(xi+1);
        lab(li+1:li+2) = x0+[-1;+1]*norm(x(xi+2:xi+kk));
        xi = xi + kk;
        li = li + 2;
    end
    lab_c{end+1} = lab * sqrt(0.5);
end
if isfield(K,'r') && ~isempty(K.r),
    % This is a simpler formula than the one found in eigK.c. In theory
    % there could be cancellation error in the smaller eigenvalue. But 
    % the rotated Lorentz vector is not used internally where this
    % cancellation error might matter.
    lab = zeros(2*length(K.r),1);
    li = 0;
    for k = 1:length(K.r)
        kk = K.r(k);
        x1 = xx(xi+1);
        x2 = xx(xi+2);
        lab(li+1:li+2) = x1+x2+[-1;+1]*norm([x1-x2;2*x(xi+3:xi+kk)]);
        xi = xi + kk;
        li = li + 2;
    end
    lab_c{end+1} = lab * 0.5;
end
if isfield(K,'s') && ~isempty(K.s),
    lab = zeros(sum(K.s),1);
    li = 0;
    Ks = K.s;
    Kq = K.s .* K.s;
    nc = length(Ks);
    % When used internally, Hermitian terms are broken apart into real and
    % imaginary halves, so we need to catch this.
    if isfield(K,'rsdpN'),
        nr = K.rsdpN;
    else
        nr = nc;
    end
    for i = 1 : nc,
        ki = Ks(i);
        qi = Kq(i);
        XX = x(xi+1:xi+qi); xi=xi+qi;
        if i > nr,
            XX = XX + 1i*x(xi+1:xi+qi); xi=xi+qi;
        end
        XX = reshape(XX,ki,ki);
        lab(li+1:li+ki) = 0.5 * eig( XX + XX' );
        li = li + ki;
    end
    lab_c{end+1} = lab;
end
lab = vertcat(lab_c{:});

