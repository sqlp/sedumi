function [At,b,c,K,prep,origcoeff] = pretransfo(At,b,c,K,pars)

% [At,b,c,K,prep] = pretransfo(At,b,c,K)
%
% PRETRANSFO  Checks data and then transforms into internal SeDuMi format.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

% Significant rewrite for SeDuMi 1.3 by Michael Grant.
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

% -------------------------------------------------------------------------
% Make sure that all fields exist in K, and verify that they are valid
% -------------------------------------------------------------------------

if ~isfield(K,'f') || isempty(K.f)
    K.f = 0;
elseif numel(K.f) ~= 1 || K.f ~= floor(K.f) || K.f < 0 || ~isreal(K.f)
    error('K.f should be nonnegative integer')
end
if ~isfield(K,'l') || isempty(K.l)
    K.l = 0;
elseif K.l ~= floor(K.l) || K.l < 0 || ~isreal(K.l)
    error('K.l should be nonnegative integer')
end
if ~isfield(K,'q') || ~nnz(K.q)
    K.q = zeros(1,0);
    L_q = 0;
    N_q = 0;
else
    K.q = K.q(:)';
    L_q = length(K.q);
    N_q = sum(K.q);
    if any(K.q ~= floor(K.q)) || any(K.q<2) || ~isreal(K.q)
        error('K.q should contain only integers bigger than 1')
    end
end
if ~isfield(K,'r') || ~nnz(K.r)
    K.r = zeros(1,0);
    L_r = 0;
    N_r = sum(K.r);
else
    K.r = K.r(:)';
    L_r = length(K.r);
    N_r = sum(K.r);
    if any(K.r ~= floor(K.r)) || any(K.r<3) || ~isreal(K.r)
        error('K.r should contain only integers bigger than 2')
    end
end
if ~isfield(K,'s') || ~nnz(K.s)
    K.s = zeros(1,0);
    L_s = 0;
    N_s = 0;
else
    K.s = K.s(:)';
    L_s = length(K.s);
    N_s = sum(K.s.^2);
    if any(K.s ~= floor(K.s)) || any(K.s<1) || ~isreal(K.s)
        error('K.s should contain only positive integers')
    end
end
N_fl = K.f + K.l;
N_qr = N_q + N_r;
L_qr = L_q + L_r;
if ~isfield(K,'ycomplex') || isempty(K.ycomplex)
    K.ycomplex = zeros(1,0);
else
    K.ycomplex = sort(K.ycomplex(:))';
    if any(K.ycomplex ~= floor(K.ycomplex)) || any(K.ycomplex<1)
        error('K.ycomplex should be a list containing only positive integers');
    elseif any(K.ycomplex > numel(b))
        error('Elements of K.ycomplex are out of range');
    elseif ~all(diff(K.ycomplex))
        K.ycomplex(find(~diff(K.ycomplex))+1) = [];
    end
end
if ~isfield(K,'xcomplex') || isempty(K.xcomplex)
    K.xcomplex = zeros(1,0);
else
    K.xcomplex = sort(K.xcomplex(:))';
    if any(K.xcomplex ~= floor(K.xcomplex)) || any(K.xcomplex<1)
        error('K.xcomplex should be a list containing only positive integers');
    elseif any(K.xcomplex > N_fl + N_qr)
        error('Elements of K.xcomplex are out of range');
    elseif ~all(diff(K.xcomplex))
        K.xcomplex(find(~diff(K.xcomplex))+1) = [];
    end
end
if ~isfield(K,'scomplex') || isempty(K.scomplex)
    K.scomplex = zeros(1,0);
else
    K.scomplex = sort(K.scomplex(:))';
    if any(K.scomplex~=floor(K.scomplex)) || any(K.scomplex<1)
        error('K.scomplex should be a list containing only positive integers');
    elseif any(K.scomplex>L_s)
        error('Elements of K.xcomplex are out of range');
    elseif ~all(diff(K.scomplex))
        K.scomplex(find(~diff(K.scomplex))+1) = [];
    end
end
N = N_fl + N_qr + N_s;

% -------------------------------------------------------------------------
% Verify the size and validity of At, b, and c
% N = # variables, m = # constraints
% -------------------------------------------------------------------------

% SeDuMi assumes that if At is not consistent with K but At' is, that the
% user supplied the transpose of the coefficient matrix. This introduces a
% rare ambiguity in the case where At happens to be square. Past versions 
% of SeDuMi would not have allowed this, instead rejecting the case where
% m >= N; however, with complex variables, the situation is not quite as
% clear, and there may technically be cases where m >= N is acceptable.
if ndims(At) > 2 %#ok
    error('A must be a matrix');
elseif nnz(isnan(At)) || nnz(isinf(At)),
    error('A contains NaN or Inf');
elseif size(At,1) == N,
    % nothing
elseif size(At,2) == N,
    At = At';
else
    error('(At,K) size mismatch');
end
if all(size(b)>1)
    error('Parameter b must be a vector');
elseif any(isnan(b)) || any(isinf(b)),
    error('b contains NaN or Inf');
elseif length(b) ~= size(At,2),
    error('(At,b) size mismatch');
else
    b = b(:);
end
if all(size(c)>1)
    error('Parameter c must be a vector');
elseif any(isnan(c)) || any(isinf(c))
    error('c contains NaN or Inf');
elseif length(c) ~= N
    error('(c,K) size mismatch');
else
    c = c(:);
end

% -------------------------------------------------------------------------
% Save the standardized data for further use if needed
% -------------------------------------------------------------------------

if isfield(pars,'errors') && pars.errors==1
    origcoeff.At=At;
    origcoeff.c=c;
    origcoeff.b=b;
    origcoeff.K=K;
else
    origcoeff=[];
end

% -------------------------------------------------------------------------
% Flag diagonal SDP blocks for removal
% -------------------------------------------------------------------------
if ~isfield(pars,'sdp')
    pars.sdp = 1;    % Enable/disable SDP preprocessing
end
if pars.sdp && L_s,
    % Some serious MATLAB trickery here if I do say so myself. Here's what
    % we are doing. "spattern" contains the indices of the nonzero elements
    % of the dual SDP variables. These are global indices over all SDPs, so
    % we use "sblk" to tell us which SDP the index belongs to. The "rem"
    % statement determines whether each element is off-diagonal. If even
    % one off-digonal element is detected for SDP #k, then sdiag(k)=false.
    ssiz        = (K.s).^2;
    strt        = cumsum([1,ssiz(1:end-1)]);
    sblk        = zeros(1,N_s);
    sblk(strt)  = 1;
    sblk        = cumsum(sblk);
    spattern    = find(c(N_fl+N_qr+1:end)|any(At(N_fl+N_qr+1:end,:),2))';
    sblk        = sblk(spattern);
    sblk        = sblk(rem(spattern-strt(sblk),K.s(sblk)+1)~=0);
    sdiag       = true(1,L_s);
    sdiag(sblk) = false;
else
    % Even if we disable SDP processing, we're going to move 1x1 SDPs into
    % the nonnegative variable block. It just doesn't make sense to deploy
    % all of that SDP machinery for nonnegative variables.
    sdiag = K.s == 1;
end
if any(sdiag),
    K.scomplex(sdiag(K.scomplex)) = [];
    prep.sdiag = K.s(sdiag);
end

% -------------------------------------------------------------------------
% Handle K.ycomplex by splitting apart the complex constraints into pairs
% of real constraints
% -------------------------------------------------------------------------

if ~isempty(K.ycomplex)
    b  = [ real(b) ; imag(b(K.ycomplex)) ];
    At = [ At, 1j * At(:,K.ycomplex) ];
else
    b  = real(b);
end

% -------------------------------------------------------------------------
% Find the locations of the the complex data, so we can convert into
% SeDuMi's internal format, which uses only MATLAB's real representation.
% -------------------------------------------------------------------------
% Strictly speaking, nonnegative variables, the first variable in a Lorentz 
% cone, and the first two variables in a rotated Lorentz cone are real. But
% But SeDuMi has allowed them all to be specified as complex; the imaginary
% portions are interpreted as free variables. We have kept that behavior.

% This code replaces whichcpx.c in its entirety. It actually fixes a bug in
% rotated Lorentz cone handling that was probably never exercised.
if isempty(K.xcomplex),
    K.fcplx = zeros(1,0);
    K.qcplx = zeros(1,0);
    K.rcplx = zeros(1,0);
    K.qc = zeros(1,L_q);
    K.rc = zeros(1,L_r);
else
    xc = K.xcomplex;
    tt = xc <= N_fl;
    K.fcplx = xc(tt);
    xc = xc(~tt) -  N_fl;
    tt = xc <= N_q;
    K.qcplx = xc(tt);
    xc = xc(~tt) - N_q;
    tt = xc <= N_r;
    K.rcplx = xc(tt);
    if isempty(K.qcplx),
        K.qc = zeros(1,L_q);
    else
        ndxs = cumsum([1,K.q(1:end-1)]);
        t2 = any(bsxfun(@eq,K.qcplx,ndxs'),1);
        K.fcplx = [ K.fcplx, K.qcplx(t2) + N_fl ];
        K.qcplx(t2) = [];
        t2 = sum(bsxfun(@gt,K.qcplx,ndxs'),1);
        K.qc = full(sparse(1,t2,1,1,L_q));
    end
    if ~isempty(K.rcplx),
        K.rc = zeros(1,L_r);
    else
        ndxs = cumsum([1,K.r(1:end-1)]);
        t2 = any(bsxfun(@eq,K.rcplx,[ndxs,ndxs+1]'),1);
        K.fcplx = [ K.fcplx, K.rcplx(t2) + N_fl + N_q ];
        K.rcplx(t2) = [];
        t2 = sum(bsxfun(@gt,K.rcplx,ndxs'),1);
        K.rc = full(sparse(1,t2,1,1,L_r));
    end
end
N_isdp = sum(K.s(K.scomplex).^2);
K.cdim = length(K.xcomplex) + N_isdp;
K.rsdpN = length(K.s) - length(K.scomplex);
K.f = K.f;
N_fc = length(K.fcplx);
N_f  = K.f + N_fc;
N_qr = N_q;
N_qc = length(K.qcplx);
N_q  = N_qr + N_qc;
N_rr = N_r;
N_rc = length(K.rcplx);
N_r  = N_rr + N_rc;

% -------------------------------------------------------------------------
% We have significantly rewritten this section of the code. This section
% constructs a sparse matrix that represents the following transformations:
%   --- Free variables split into differences of nonnegative variables; OR
%   --- Free variables placed in a Lorentz-cone
%   --- Rotated lorentz cones translated to standard Lorentz cones
%   --- Lorentz cones rearranged to trace block + norm-bound blocks
%   --- Conversion of diagonal SDPs to nonnegative variables
%   --- SDP coefficients moved to the lower triangle for increased sparsity
% This code replaces the rotlorenz, qreshape, and vectril MEX files.
% -------------------------------------------------------------------------

% This reserves space for diagonal SDPs
newL = sum(K.s(sdiag)); 
newQ = [];
ii = {}; jj = {}; vv = {};

% Translate free variables into split variables or a Lorentz cone
if N_f,
    if ~isfield( pars, 'free' ) || pars.free == 2 && ( L_q + L_r + L_s ) || pars.free == 1,
        ii{end+1} = K.l + L_qr + 2 : K.l + L_qr + 1 + N_f;
        jj{end+1} = [ 1 : K.f, K.fcplx ];
        vv{end+1} = [ ones(1,K.f), -1j * ones(1,N_fc) ];
        newQ = N_f + 1;
        prep.freeQ = N_f;
    else
        ii{end+1} = 1 : 2 * N_f;
        jt = [ [ 1 : N_f ; 1 : N_f ], [ K.fcplx, K.fcplx ] ];
        jj{end+1} = jt(:)';
        vt = [ [ ones(1,K.f) ; -ones(1,K.f) ], [ -1j * ones(1,N_fc), 1j * ones(1,N_fc) ] ];
        vv{end+1} = vt(:)';
        newL = newL + 2 * N_f;
        prep.freeL = N_f;
    end
end

% Copy nonnegative variables without change
if K.l,
    ii{end+1} = newL + 1 : newL + K.l;
    jj{end+1} = K.f + 1 : K.f + K.l;
    vv{end+1} = ones(1,K.l);
    newL = newL + K.l;
end

% Rearrange Lorentz cones to trace block + norm-bound blocks
tr_off = newL + length(newQ);
nb_off = newL + L_qr + sum(newQ);
if N_q,
    ndxs = cumsum([1,K.q(1:end-1)+K.qc(1:end-1)]);
    it        = zeros(1,N_q);
    it(ndxs)  = tr_off + 1 : tr_off + L_q;
    it(it==0) = nb_off + 1 : nb_off + ( N_q - L_q );
    ii{end+1} = it(:)';
    if any(K.qc),
        ndxc = K.qcplx + ( 1 : length(K.qcplx) );
        jt = [ K.f + K.l, ones(1,N_q-1) ];
        jt(ndxc) = 0;
        jt = cumsum(jt);
    else
        jt = K.f + K.l + 1 : K.f + K.l + N_q;
        ndxc = [];
    end
    jj{end+1} = jt;
    vt = ones(1,N_q);
    vt(ndxc)  = -1j;
    vv{end+1} = vt;
    tr_off    = tr_off + L_q;
    nb_off    = nb_off + N_q - L_q;
    newQ      = [ newQ, K.q + K.qc ];
end

% Transform rotated Lorentz cones to standard Lorentz cones, and rearrange
% to trace block + norm-bound blocks.
if N_r,
    ndxr = cumsum([1,K.r(1:end-1)+K.rc(1:end-1)]);
    ndxp = ndxr + 2*(0:L_r-1);
    it         = zeros(1,N_r+2*L_r);
    it(ndxp)   = tr_off + 1 : tr_off + L_r;
    it(ndxp+1) = -1;
    it(ndxp+2) = it(ndxp);
    it(it==0)  = nb_off + 1 : nb_off + ( N_r - L_r );
    it(ndxp+1) = it(ndxp+3);
    ii{end+1}  = it;
    jt = [ K.f + K.l + N_qr, ones(1,N_r+2*L_r-1) ];
    jt([ndxp+1,ndxp+3]) = 0;
    if any(K.rc),
        ndxc = K.rcplx + ( 1 : length(K.rcplx) );
        ndxc = ndxc + 2 * sum(bsxfun(@gt,ndxc,ndxr(:)));
        jt(ndxc) = 0;
    else
        ndxc = [];
    end
    jj{end+1} = cumsum(jt);
    vt = ones(1,N_r+2*L_r);
    vt([ndxp,ndxp+1,ndxp+2]) = sqrt(0.5);
    vt(ndxp+3) = -sqrt(0.5);
    vt(ndxc) = -1j;
    vv{end+1} = vt;
    nb_off = nb_off + N_r - L_r;
    newQ = [ newQ, K.r + K.rc ];
end

% For diagonal SDPs, use the space reserved above to convert them to
% nonnegative variablees. Otherwise, replace the coefficients with
% tril(X) + tril(X',-1)
if K.rsdpN,
    ioff = nb_off;
    joff = N_fl + N_qr;
    loff = 0;
    for k = 1 : L_s + L_z,
        kk = K.s(k);
        kq = kk * kk;
        if sdiag(k),
            ii{end+1} = loff + 1 : loff + kk; %#ok
            jj{end+1} = joff + 1 : kk + 1 : joff + kq; %#ok
            vv{end+1} = ones(1,kk); %#ok
            loff = loff + kk;
        elseif ~K.scomplex(k),
            cc = 0 : kk - 1; rr = cc';
            is = bsxfun(@max,rr,cc) + bsxfun(@min,rr,cc) * kk + 1;
            is = is(:)';
            ii{end+1} = ioff + is; %#ok
            jj{end+1} = joff + 1 : joff + kq; %#ok
            vv{end+1} = ones(1,kq); %#ok
        end
        ioff = ioff + kq;
        joff = joff + kq;
        if k > K.rsdpN,
            ioff = ioff + kq;
            joff = joff + kq;
        end
    end
end

% Replace the Hermitian semidefinite coefficients with tril(X)+tril(X',-1)
if K.rsdpN < length(K.s),
    joff = N_fl + N_qr;
    for k = 1 : L_s + L_z,
        kk = K.s(k);
        kq = kk * kk;
        if K.scomplex(k),
            cc = 0 : kk - 1; rr = cc';
            is = bsxfun(@max,rr,cc) + bsxfun(@min,rr,cc) * kk + 1;
            js = joff + 1 : joff + kq;
            is = is(:)';
            ii{end+1} = ioff + is; %#ok
            jj{end+1} = joff + 1 : joff + kq; %#ok
            vv{end+1} = ones(1,kq); %#ok
            vs = -1j * sign(bsxfun(@minus,rr,cc));
            tt = vs ~= 0;
            ii{end+1} = kq + is(tt); %#ok
            jj{end+1} = kq + js(tt); %#ok
            vv{end+1} = vs(tt); %#ok
        end
        ioff = ioff + kq;
        joff = joff + kq;
    end
end

% Update free, nonnegative, and Lorentz variable counts
K.f = 0;
K.l = newL;
K.q = newQ;
K.s = [ K.s(~K.scomplex&~sdiag), K.s(K.scomplex&~sdiag) ];
K.rsdpN = nnz(~K.scomplex&~sdiag);
K.N = K.l + sum(K.q) + sum(K.s(1:K.rsdpN).^2)+2*sum(K.s(K.rsdpN+1:end).^2);

% Create the artificial (x0,z0) variable for the self-dual model by
% appending a zero row to At and c. This is accomplished in QR by adding
% 1 to all of the row indices created above.
K.N  = K.N + 1;
K.l  = K.l + 1;
K.m  = length(b);

% Transform At, c
% The transformation matrix QR does not satisfy QR'*QR=I. But it does, in 
% fact, serve as a reverse transformation:
%    --- For Lorentz cones, QR applies a permutation; QR' reverses it.
%    --- For rotated Lorentz cones, QR applies a unitary, self-adjoint
%        rotation on the first two variables, so QR' reverses it.
%    --- For split free variables, QR creates the positive and negative
%        parts; QR' combines them back together.
%    --- For free variables placed in a Lorentz cone, QR moves them to the
%        Lorentz cone block, adding an extra epigraph variable. QR' moves
%        them back and drops the extra variable.
%    --- For semidefinte cones, QR adds the strict upper triangle to the
%        strict lower triangle; QR' copies the strict lower triangle to the
%        strict upper triangle, ensuring symmetry.
%    --- QR adds a row for the self-dual variable; QR' removes it.
[dummy,ndxs] = sort(cellfun(@(x) x(1),jj)); %#ok
QR = sparse( horzcat(ii{ndxs})+1, horzcat(jj{ndxs}), horzcat(vv{ndxs}), K.N, length(c) );
At = real( sparse( QR * At ) );
c  = real( sparse( QR * c  ) );
b  = sparse( b );
prep.QR = QR;
clear ii jj vv

% -------------------------------------------------------------------------
% Now K has field K.{l,q,s}
% Generate a more detailed description of cone K:
% Let K.blkstart(i):K.blkstart(i+1)-1 give the index range of block i.
% Compute maxN=max(order), len=sum(order) for LORENTZ, real PSD, herm PSD
% yields: K.{l,q,s,rsdpN,blkstart,rLen,hLen,qMaxn,rMaxn,hMaxn}
% -------------------------------------------------------------------------
Ksr = K.s(1:K.rsdpN);
Ksc = K.s(K.rsdpN+1:end);
K.blkstart = cumsum([K.l+1,length(K.q)+length(K.r),K.q-1,Ksr.^2,2*Ksc.^2]);
K.rLen = sum(Ksr);
K.hLen = sum(Ksc);
K.qMaxn = max([0,K.q]);
K.rMaxn = max([0,Ksr]);
K.hMaxn = max([0,Ksc]);
K.mainblks = K.blkstart(cumsum([1 1 length(K.q)]));
K.qblkstart = K.blkstart(2:2+length(K.q));  % Also include blkend
K.sblkstart = K.blkstart(2+length(K.q):end);
K.lq = K.mainblks(end)-1;
