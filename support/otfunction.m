function [marginals,q,b,wes_dis] = otfunction(P, M, b,epsilon, gamma, tol,maxiter)
[n_features, n_tasks ]= size(P);
frac = gamma / (gamma + epsilon);
support = any(P,2);
P = P(support,:);
M=M(support,support);
M = exp(M);
if isempty(b)
    b = ones(size(P,1), n_tasks);
end
Kb = M*b;
weights =ones(n_tasks,1) / n_tasks;
q = ones(size(P,1),1);
qold =ones(size(P,1),1);
cstr = 1.;
qmax_old = 1.;
for i=1 :maxiter
        a = (P ./ Kb).^frac;
        Ka =M'*a;
        q = (Ka .^ (1 - frac))*weights;
        q = q .^ (1 / (1 - frac));
        Q = q;
        
        qmax = max(max(q));
        if i > 2
            cstr = max(abs(q - qold)) / max([qmax_old, qmax, 1.]);
        end
        qold = q;
        qmax_old = qmax;
        b_old = b;
        b = (Q ./ Ka).^ frac;
        Kb = M*b;
        if abs(cstr) < tol && i > 2
            break
        end      
end
wes_dis=zeros(1,n_tasks);
for t=1:n_tasks
    wes_dis(1,t)=norm(a(:,t)*b(:,t)'.*M.*M,1);
end
marginals = (a .* Kb)';
m = zeros(n_tasks, n_features);
m(:, support) = marginals;
marginals = m;


