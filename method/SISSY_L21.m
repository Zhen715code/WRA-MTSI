function [s, history] = l2_l2_l2(L, b, v,  lambda1,lambda2,rho1,rho2,alpha)
% 0.0001,0.001,0.01,1,0.01
%
% Solves the following problem via ADMM:
%
% minimize  1/(2 dc)|| Ls - b ||_2_2 + lambda1|| s ||_2+lambda2|| vs ||_2
%
% The solution is returned in the vector s.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% lambda1,lambda2 is the augmented Lagrangian parameter.
%
% 


% %-------------------- Global constants and defaults---------------------
QUIET    = 0;
MAX_ITER = 5000;
ABSTOL   = 1e-4;%绝对容忍度
RELTOL   = 1e-2;%相对容忍度
% % --------------------Data preprocessing------------------------------------
[dc,ds] = size(L);
[dp,ds]=size(v);
n_dic=size(b,2);

% % ------------------------ADMM solver------------------------------------
s = zeros(ds,n_dic);
u = zeros(dp,n_dic);
w = zeros(ds,n_dic);
x = zeros(ds,n_dic);
y = zeros(dp,n_dic);


if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm1', 'eps pri1', 'r norm2', 'eps pri2' ,...
      's norm1', 'eps dual1', 's norm2', 'eps dual2','objective');
end

R = inv(1/dc*(L'*L)+rho2*(v'*v)+rho1*eye(ds));

for k = 1:MAX_ITER
    
    % s-update
   
    s = R*((1/dc)*L'*b-x+rho2*v'*u+rho1*w-v'*y);
    
    % u-update
    uold = u;
    vs_hat = alpha*v*s+(1-alpha)*uold;
    u = update_2_2(vs_hat,rho2,lambda2,y);
    
    % w-update
    wold =w;
    w_hat = alpha*s+(1-alpha)*wold;
    w = update_2_2(w_hat,rho1,lambda1,x);
    
    % x-update
    x=x+rho1*(w_hat-w);
    
     % y-update
     y = y+rho2*(vs_hat-u);

    % diagnostics, reporting, termination checks
     history.r1_norm(k)  =norm(s-w);%主残差r1
     history.s1_norm(k)  = norm(-rho1*(w - wold));%对偶残差s1
     
     history.eps_pri1(k) = sqrt(ds)*ABSTOL + RELTOL*max(norm(s), norm(-w));
     history.eps_dual1(k)= sqrt(ds)*ABSTOL + RELTOL*norm(x); %不加 rho1 
    
     history.r2_norm(k)  =  norm(v*s-u);%主残差r2
     history.s2_norm(k)  = norm(-rho2*v'*(u - uold));   
 
     history.eps_pri2(k) = sqrt(dp)*ABSTOL + RELTOL*max([norm(v*s), norm(-u)]);
     history.eps_dual2(k)= sqrt(ds)*ABSTOL + RELTOL*norm(v'*y); %不加 rho2 
    

    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n',...
        k,  history.r1_norm(k), history.eps_pri1(k), ...
            history.r2_norm(k), history.eps_pri2(k), ...
            history.s1_norm(k), history.eps_dual1(k), ...
            history.s2_norm(k), history.eps_dual2(k));
        fprintf('\n')
    end

    if ( history.r1_norm(k) < history.eps_pri1(k) && ...
         history.s1_norm(k) < history.eps_dual1(k) && ...
         history.r2_norm(k) < history.eps_pri2(k) && ...
         history.s2_norm(k) < history.eps_dual2(k)  )
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
