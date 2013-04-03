%                                                 y = triumtriu(r,u,K)
% TRIUMTRIU  Computes y = r * u
%   Both r and u should be upper triangular.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function y = triumtriu(r,u,K)
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

Ks=K.s;
startindices=K.sblkstart-K.mainblks(end)+1;
y=zeros(K.blkstart(end)-K.mainblks(end),1);
for k = 1:length(Ks)
%This works, but I don't like the lot of 0 matrices.
    temp=triu(reshape(r(startindices(k):startindices(k+1)-1),Ks(k),Ks(k)),0)*triu(reshape(u(startindices(k):startindices(k+1)-1),Ks(k),Ks(k)),0);
    y(startindices(k):startindices(k+1)-1)=temp+triu(temp,1)';
end