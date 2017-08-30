% Guassian Sampling
function GM = guassSamp(M, q)

% Tensor Sampling based on Gaussian measurements
%
% Input:
%       M       -    d1*d2*d3 tensor, the tensor to be sampled from
%       q       -    int, the number of sampled obeservations
%
% Output:
%       GM = [g_1,...,g_q] - the sampling vector

if nargin < 3
    q = floor(0.1*numel(X));
end
for i=1:q
    g{i}
[n1, n2, n3] = size(M);


