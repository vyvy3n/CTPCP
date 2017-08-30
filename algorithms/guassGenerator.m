% Guassian Sampling Generator
function G = guassGenerator(M, q)

% Tensor Sampling based on Gaussian measurements
%
% Input:
%       M       -    d1*d2*d3 tensor, the tensor to be sampled from
%       q       -    int, the number of sampled obeservations
%
% Output:
%       G{i}    - the sampling tensor 1,...,q

if nargin < 3
    q = floor(0.1*numel(X));
end

for i=1:q
    G{i}=normrand(0,1,size(M));
end

return G


