function [s,m,q,b,d]= WRA_MTSI(L, Y_all, V, M,b, gamma,epsilon,lambda1, lambda2, rho,mu,eta1,eta2,s,m,alpha,n_subject)
%add rho update
[dc,ds] = size(L);
n_dic=1;
p=size(V,1);
s_m=zeros(ds,n_subject);
max_iter_ot = 20;
tol_ot = 1e-4;
ABSTOL = 1e-4;
RELTOL = 1e-2;
MAX_ITER =300;
CoreNum=6; %设定机器CPU核心数量
if isempty(gcp('nocreate'))
    parpool(CoreNum);
end

parfor j=1:n_subject
%for j=1:n_subject
    rou=rho;
    rou_old = rou;
    Z=zeros(p,n_dic);
    U=zeros(p,n_dic);
    LL=sum(L.*L,1);
    VV=sum(V.*V,1);
    mu_gamma=mu*gamma/n_subject;
    mu_gamma_lamb=lambda1+mu_gamma;
    mu_gamma_m=mu_gamma*m(j,:);
    dc1=1/dc;
    R1=dc1*LL+rou*VV;
    Y=Y_all((j-1)*dc+1:j*dc,:);
    for k=1:MAX_ITER
        % s-update
        s1=s(:,j);
        %       tic;
        [s(:,j)]=update_s(L, Y,V,Z,U,eta1,eta2,rou,s1,dc1,ds,...
            R1,mu_gamma_lamb,mu_gamma_m);
        s1=s(:,j);
        % z-update
        Z_old=Z;
        V_s1=V*s1;
        V_s=alpha*V_s1+(1-alpha)*Z;
        Z=update_Z(V_s,Z,U,rou,lambda2);
        % u-update
        
        U=U + rou * (V_s-Z);
        
        % Stopping condition
        r1_norm  = norm(V_s1-Z);%主残差r1
        s1_norm  =norm(-rou * V'*(Z -Z_old));%对偶残差s1
        
        eps_pri1 = sqrt(p)*ABSTOL + RELTOL*max(norm(V_s1),norm(-Z));
        eps_dual1= sqrt(ds)*ABSTOL + RELTOL*norm(V'*U);
        
        
        fprintf('%5d\t%5d\t%f\t%f\t%f\t%f\n',j,k,r1_norm, eps_pri1,s1_norm, eps_dual1);
        if ( r1_norm < eps_pri1 && ...
                s1_norm < eps_dual1)
            %             fprintf('%5d\t%5d\n',j,k)
            break;
        end
        
        %---------------- rou update ---------------------------------%
         if mod(k,10) == 0
             if r1_norm > 10 *  s1_norm
                 rou = 2*rou;
            elseif s1_norm > 10 * r1_norm
                 rou = rou/2;
             end
             if rou ~= rou_old
                 R1=dc1*LL+rou*VV;
                 fprintf('%s%f\t\n','rou:',rou);
             end
             rou_old = rou;
         end
        
    end
    s_m(:,j)=sqrt(sum(s(:,j).*s(:,j),2));
end

if all(any(s_m,1),1)
    fprintf('m--------\n')
    [m,q,b,d] = otfunction(s_m, M, b,epsilon,gamma,tol_ot,max_iter_ot);
end
end
function [s]=update_s(L, Y,V,Z,U,eta1,eta2,rho,s,dc1,ds,R1,mu_gamma_lamb,mu_gamma_m)
s_2_2=sum(s.*s,2)';
Ls=L*s;
Vs=V*s;
P1=R1+mu_gamma_lamb*(1./(sqrt(s_2_2)+eta1*ones(1,ds)))-mu_gamma_m.*(1./(s_2_2+eta2*ones(1,ds)));
for i=1:ds
    L_S=Ls-L(:,i)*s(i,:);
    V_S=Vs-V(:,i)*s(i,:);
    P2=dc1*L(:,i)'*(Y-L_S)-V(:,i)'*U+rho*V(:,i)'*(Z-V_S);
    s(i,:)=P2/P1(:,i);
    Ls=L_S+L(:,i)*s(i,:);
    Vs=V_S+V(:,i)*s(i,:);
end
end

function Z=update_Z(V_s,~,U,rho,lambda2)
W=V_s+(1/rho)*U;
lamb_subject=lambda2/rho;
z1=max(0,(1-lamb_subject./sqrt(sum(W.*W,2))));
Z=W.*z1;
end
