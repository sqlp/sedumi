function lab = mineigK(x,K)

% lab = mineigK(x,K)
%
% MINEIGK  Computes the minimum eigenvalue of a vector x with respect to a 
%          self-dual homogenous cone K.
%
% See also sedumi, mat, vec, eyeK.

% New function by Michael C. Grant
% Copyright (C) 2013 Michael C. Grant.

xi = 0;
if isfield(K,'f'),
    xi = xi + K.f;
end
if isfield(K,'l') && K.l > 0,
    lab = min(x(xi+1:xi+K.l));
    xi = xi+K.l;
else
    lab = Inf;
end
if isfield(K,'q') && ~isempty(K.q),
    scl = sqrt(0.5);
    for k = 1:length(K.q)
        kk = K.q(k);
        lab = min(lab,scl*(x(xi+1)-norm(x(xi+2:xi+kk))));
        xi = xi + kk;
    end
end
if isfield(K,'r') && ~isempty(K.r),
    % This is a simpler formula than the one found in eigK.c. In theory
    % there could be cancellation error in the smaller eigenvalue. But 
    % the rotated Lorentz vector is not used internally where this
    % cancellation error might matter.
    for k = 1:length(K.r)
        kk = K.r(k);
        x1 = xx(xi+1);
        x2 = xx(xi+2);
        lab = min(lab,0.5*(x1+x2-norm([x1-x2;2*x(xi+3:xi+kk)])));
    end
end
if isfield(K,'s') && ~isempty(K.s),
    Ks = K.s;
    Kq = K.s .* K.s;
    nc = length(Ks);
    OPTS.disp=0;
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
        XX = XX + XX';
        if ki > 500
            lab=min(lab,0.5*eigs(XX,1,'SA',OPTS));
        else
            lab=min(lab,0.5*min(eig(XX)));
        end
    end
end

