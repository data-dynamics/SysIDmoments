
function ZZ = momentPowers(n,d)
% n: # of variables
% d: degree of polynomial
% this file is part of SOSTOOLS, modified by NO
% generates the moments of at most degree d in R^n in lexicographical order
% and outputs an M(d) by n list containing powers
% here M(d) = nchoosek(n+d,n)
ZZ = sparse(1,n);
for i = 1:n
    ss = size(ZZ,1);
    ZZ = sprepmat(ZZ,d+1,1);
    for j = 0:d
        ZZ(ss*j+1:ss*j+ss,i) = j;
    end;
    idx = find(sum(ZZ,2) <= d);   % Throw away invalid monomials
    ZZ = ZZ(idx,:);
end;
ZZ = full(ZZ);
%idx = find(sum(ZZ,2) == d);
%Z = ZZ(idx,:);