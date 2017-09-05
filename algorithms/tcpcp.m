function [L,S,obj,err,iter] = tcpcp(dim,g,GM,lambda,opts)
%
% Tensor Compressive Principal Component Pursuit Algomrith
%
% Description: to solve Tensor Compressive Principal Component Analysis 
% based on Tensor Nuclear Norm problem by ADMM
%
% min_{X} ||L||_*+lambda*||S||_1, s.t. P_G(L+S)=P_G(L+S),
%
% where M is the original matrix,
%       P_G is sampling based on Guassian Measurement.
%
% ---------------------------------------------
% Input:
%       X       -    d1*d2*d3 tensor
%       lambda  -    >0, parameter
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%
% Output:
%       L       -    d1*d2*d3 tensor
%       S       -    d1*d2*d3 tensor
%       obj     -    objective function value
%       err     -    residual 
%       iter    -    number of iterations


tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end

L = zeros(dim);
S = zeros(dim);
P = L;
Q = S;
Z1 = L;
Z2 = S;
m = prod(dim);

iter = 0;
for iter = 1 : max_iter
    Lk = L;
    Sk = S;
    % update P
    [P,tnnP] = prox_tnn(L+Z1,1/mu);
    % update Q
    S = prox_l1(S+Z2,lambda/mu);
    % update L
    GM2 = GM'*GM;
    L = reshape(pinv(GM2+eye(m))*(GM'*g+P(:)-GM2*S(:)),dim);
    % update S
    S = reshape(pinv(GM2+eye(m))*(GM'*g+Q(:)-GM2*L(:)),dim);
    % dual update difference
    dZ1 = L-P;
    dZ2 = S-Q;
    chgL = max(abs(Lk(:)-L(:)));
    chgS = max(abs(Sk(:)-S(:)));
    chg = max([ chgL chgS max(abs(dZ1(:))) max(abs(dZ2(:))) ]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnP+lambda*norm(S(:),1);
            err = norm(dZ1(:))+norm(dZ2(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    if chg < tol
        break;
    end 
    % dual update Z1
    Z1 = Z1 + dZ1;
    % dual update Z2
    Z2 = Z2 + dZ2;
    mu = min(rho*mu,max_mu);    
end
obj = tnnP+lambda*norm(S(:),1);
err = norm(dZ1(:))+norm(dZ2(:));
