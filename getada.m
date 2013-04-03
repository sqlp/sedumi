function absd=getada(A,K,d,DAt)
%Computes the Newton system
%incorporates getada1,2,3
% INPUT:
% ADA: only the sparsity pattern is used, not the values
% A: The coefficient matrix, actually A', with the dense columns (actually rows) removed
% K: The structure of the cone
% d: The scaling elements
% DAt: The scaled A matrix used in getada2

global ADA
%This is why ADA is not in the function calling sequence, it's global
m = size(ADA,1);
if spars(ADA) > 0.3
    %We convert it to a dense 0 matrix
    ADA = zeros(m,m);
else
    %Create a sparse 0 matrix to accomodate the number of nonzeros needed
    ADA = sparse([],[],[],m,m,nnz(ADA));
end
if spars(DAt.q > 0.2)
    %TODO: This conversion will have to move to getdatm.m
    DAt.q = full(DAt.q);
end
Alq = A(1:K.mainblks(3)-1,:);
if spars(Alq) > 0.2
    Alq = full(Alq);
end
scalingvector = [d.l; -d.det; zeros(K.mainblks(3)-K.mainblks(2),1)];
for i = 1:length(d.det)
    scalingvector(K.qblkstart(i):K.qblkstart(i+1)-1) = d.det(i);
end
%getada2
ADA = ADA+DAt.q'*DAt.q;
clear DAt
%getada1
ADA = ADA+Alq'*diag(sparse(scalingvector))*Alq;
clear Alq
ADA = sparse(ADA);
absd = full(diag(ADA));