function [P,Q,obj,err,iter] = tcpcp(dim,g,GM,GM2,lambda,opts)

% Tensor Compressive Principal Component Pursuit Algomrith
%
% Description: to solve Tensor Compressive Principal Component Analysis 
% based on Tensor Nuclear Norm problem by ADMM
%
% min_{L,S} ||L||_*+lambda*||S||_1, s.t. P_G(L+S)=P_G(L+S),
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

tol         = 1e-10; 
max_iter    = 500;
rho         = 1.05;
mu          = 0.001;
max_mu      = 1e10;
penalty     = 0.005;
DEBUG       = 0;
muOnly      = 1;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'penalty');     penalty = opts.penalty;      end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end
if isfield(opts, 'muOnly');      muOnly = opts.muOnly;          end

X = reshape(GM\g,dim);%
L = X;
S = zeros(dim);
P = L;
Q = S;
Z1 = zeros(dim);
Z2 = zeros(dim);
Z3 = zeros(size(g));
m = prod(dim);

iter = 0;
for iter = 1 : max_iter
    Pk = P; Qk = Q;
    % update P
    [P,tnnP] = prox_tnn(X-S+Z1/mu,1/mu);
    % update Q
    Q = prox_l1(X-L+Z2/mu,lambda/mu);
    % update L
    if muOnly==1
        [ll,~] = cgs((GM2+eye(m)),(GM'*g+P(:)-Z1(:)/mu-GM2*S(:)-GM'*Z3/mu),1e-6,200);
    else
        [ll,~] = cgs((penalty*GM2+mu*eye(m)),(penalty*GM'*g+mu*P(:)-Z1(:)-penalty*GM2*S(:)),1e-6,200);
    end
    L = reshape(ll,dim);
    % update S
    if muOnly==1
        [ss,~] = cgs((GM2+eye(m)),(GM'*g+Q(:)-Z2(:)/mu-GM2*L(:)-GM'*Z3/mu),1e-6,200); 
    else
        [ss,~] = cgs((penalty*GM2+mu*eye(m)),(penalty*GM'*g+mu*Q(:)-Z2(:)-penalty*GM2*L(:)),1e-6,200);
    end
    S = reshape(ss,dim);
    % dual update difference
    dZ1 = L-P;
    dZ2 = S-Q;
    if muOnly==1
        dZ3 = GM*(L(:)+S(:))-g;
    else
        dZ3 = Z3;
    end
    chgP = max(abs(Pk(:)-P(:)));
    chgQ = max(abs(Qk(:)-Q(:)));
    chg = max([ chgP chgQ max(abs(dZ1(:))) max(abs(dZ2(:))) max(abs(dZ3(:)))]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnP+lambda*norm(Q(:),1);
            err = norm(dZ1(:))+norm(dZ2(:))+norm(dZ3(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)...
                    ', norm(Z1)=' num2str(norm(abs(Z1(:)))) ', norm(Z2)=' num2str(norm(abs(Z2(:))))...
                    ', norm(P)=' num2str(norm(abs(P(:)))) ', norm(Q)=' num2str(norm(abs(Q(:))))...
                    ', norm(L)=' num2str(norm(abs(L(:)))) ', norm(S)=' num2str(norm(abs(S(:))))]); 
        end
    end
    
    if chg < tol
        break;
    end 
    % dual update Z1
    Z1 = Z1 + mu*dZ1;
    % dual update Z2
    Z2 = Z2 + mu*dZ2;
    if muOnly==1
        Z3 = Z3 + mu*dZ3;
    end
    X = L+S;
    mu = min(rho*mu,max_mu);    
end
obj = tnnP+lambda*norm(S(:),1);
err = norm(dZ1(:))+norm(dZ2(:))+norm(dZ3(:));
