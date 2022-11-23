function [s, history] = GL(L, b,lambda1,rho1,alpha)
% group_lasso  Solve group lasso problem via ADMM
%
% [x, history] = group_lasso(A, b, p, lambda, rho, alpha);
%
% solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda sum(norm(x_i))
%
% The input p is a K-element vector giving the block sizes n_i, so that x_i
% is in R^{n_i}.
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%

% %-------------------- Global constants and defaults---------------------
QUIET    = 0;
MAX_ITER = 5000;
ABSTOL   = 1e-4;%绝对容忍�?
RELTOL   = 1e-2;%相对容忍�?
% % --------------------Data preprocessing------------------------------------
[dc,ds] = size(L);
n_dic=size(b,2);

% % ------------------------ADMM solver------------------------------------
s = zeros(ds,n_dic);
w = zeros(ds,n_dic);
x = zeros(ds,n_dic);




if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm1', 'eps pri1', 'r norm2', 'eps pri2' ,...
      's norm1', 'eps dual1', 's norm2', 'eps dual2','objective');
end

R = inv(1/dc*(L'*L)+rho1*eye(ds));

for k = 1:MAX_ITER
    
    % s-update
   
    s = R*((1/dc)*L'*b-x+rho1*w);
    
    % w-update
    wold =w;
    w_hat = alpha*s+(1-alpha)*wold;
    w = update_2_2(w_hat,rho1,lambda1,x);
    
    % x-update
    x=x+rho1*(w_hat-w);
    
    % diagnostics, reporting, termination checks
     history.r1_norm(k)  =norm(s-w);%主残差r1
     history.s1_norm(k)  = norm(-rho1*(w - wold));%对偶残差s1
     
     history.eps_pri1(k) = sqrt(ds)*ABSTOL + RELTOL*max(norm(s), norm(-w));
     history.eps_dual1(k)= sqrt(ds)*ABSTOL + RELTOL*norm(x); %不加 rho1 
    
    

    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n',...
        k,  history.r1_norm(k), history.eps_pri1(k), ...
            history.s1_norm(k), history.eps_dual1(k));
        fprintf('\n')
    end

    if ( history.r1_norm(k) < history.eps_pri1(k) && ...
         history.s1_norm(k) < history.eps_dual1(k))
     %fprintf('%d\n',k);
         break;
    end

end


end

function Z=update_2_2(vs_hat,rho,lambda,y)
W=vs_hat+(1/rho)*y;
lamb_rho=lambda/rho;
z1=max(0,(1-lamb_rho./sqrt(sum(W.*W,2))));
Z=W.*z1;
end