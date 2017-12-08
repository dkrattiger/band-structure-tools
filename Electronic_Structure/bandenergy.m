function E_b = bandenergy(E,f,w)

% E: Band structure energies
% f: Band occupation values f(i) =  1 if ith band is occupied, = 0  if ith band is empty
% w: k-point weights


[n_bands,n_kap] = size(E);
if nargin <3
    w = ones(n_kap,1)/n_kap;
end

if nargin <2    
    f = ones(n_bands,1);
end

% column vector-ize each array;
f = f(:);
w = w(:);

% use vector-matrix-vector multiplication to sum over both indices
E_b = f'*E*w;
    