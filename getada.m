function absd=getada(Alq,K,d,DAt,pars)
%Computes the Newton system
%incorporates getada1,2,3
% INPUT:
% ADA: only the sparsity pattern is used, not the values
% Alq: Structure containing the dense and sparse rows (actually cols) of the 
%      coefficient matrix, (actually A'), with the dense columns (actually rows) removed
% K: The structure of the cone
% d: The scaling elements
% DAt: The scaled A matrix used in getada2
% pars: parameter structure (for density thresholds)

global ADA_sedumi_
%This is why ADA_sedumi_ is not in the function calling sequence, it's global
m = size(ADA_sedumi_,1);
if spars(ADA_sedumi_) > pars.spars_thold.ADA
    %We convert it to a dense 0 matrix
    ADA_sedumi_ = zeros(m,m);
else
    %Create a sparse 0 matrix to accomodate the number of nonzeros needed
    ADA_sedumi_ = sparse([],[],[],m,m,nnz(ADA_sedumi_));
end
if spars(DAt.q) > pars.spars_thold.DAt_q
    %TODO: This conversion will have to move to getdatm.m
    DAt.q = full(DAt.q);
end
scalingvector = [d.l; -d.det; zeros(K.mainblks(3)-K.mainblks(2),1)];
for i = 1:length(d.det)
    scalingvector(K.qblkstart(i):K.qblkstart(i+1)-1) = d.det(i);
end
ipos=find(scalingvector>0); % find positive and negative entries, will compute separately below
ineg=find(scalingvector<0);

%getada1
% compute dense*dense, dense*sparse, and sparse*sparse separately for speed and efficiency
% for dense matrices, matlab optimizes A'*A to take advantage of symmetry
% To take advantage, must separate out positive and negative scalingvector entries
if ~isempty(Alq.dense) % need to check because zero-dim multiples can cost cycles
   %Add=Alq_d'*bsxfun(@times,scalingvector,Alq_d); % will not take advantage of symmetry
   Adp=bsxfun(@times,sqrt(scalingvector(ipos)),Alq.dense(ipos,:));
   Adn=bsxfun(@times,sqrt(-scalingvector(ineg)),Alq.dense(ineg,:));
   Add=Adp'*Adp-Adn'*Adn;
   clear Adp Adn;
   ADA_sedumi_(Alq.drows,Alq.drows)=Add;
   clear Add;
end;
if (~isempty(Alq.dense) && ~isempty(Alq.sparse))
   Ads=Alq.dense'*(diag(sparse(scalingvector))*Alq.sparse);
   % Matlab does not take advantage of symmetry for sparse matrix multiplies, so do in one call
   ADA_sedumi_(Alq.drows,~Alq.drows)=Ads;
   ADA_sedumi_(~Alq.drows,Alq.drows)=Ads';
   clear Ads;
end;
if ~isempty(Alq.sparse)
   Ass=Alq.sparse'*diag(sparse(scalingvector))*Alq.sparse;
   ADA_sedumi_(~Alq.drows,~Alq.drows)=Ass;
   clear Ass;
end;
   
%getada2
ADA_sedumi_ = ADA_sedumi_+DAt.q'*DAt.q;
clear DAt

ADA_sedumi_ = sparse(ADA_sedumi_);
absd = full(diag(ADA_sedumi_));
