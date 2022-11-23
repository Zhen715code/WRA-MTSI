function [thet1,m1,b]= MTW(X, Y, M,thet1,thet2,b, epsilon, gamma,lamb, mu,m1,n_subject)
[n_chanel,n_sources] = size(X);
max_iter_ot = 20;
tol_ot = 1e-7;
MAX_ITER = 1000;
ftol=1e-7;
thet1=squeeze(thet1);
thet2=squeeze(thet2);

% CoreNum=6; %设定机器CPU核心数量
% if isempty(gcp('nocreate'))
%     parpool(CoreNum);
% end

parfor j=1:n_subject
% for j=1:n_subject
    thet_old =thet1(:, j) ;
    Y1=Y((j-1)*n_chanel+1:j*n_chanel,:);
    for k=1:MAX_ITER
        thet1(:, j)=update_thet(thet1(:, j),thet2(:, j),Y1,X, mu,gamma, lamb,m1(j,:),n_chanel,n_subject,n_sources);
        dx = max(abs(thet1(:, j) - thet_old))/ max([1, max(thet_old),max(thet1(:, j))]);
%         fprintf('%4s%4d\t%f\n','dx:',dx)
%         fprintf('\n');
        thet_old=thet1(:, j);
        if  dx<ftol
            %fprintf('%4d%4d\n',j,k);
            break;
        end
    end
end

if all(any(thet1>1e-10,1),1)
%     fprintf('\n')
    %fprintf('m--------\n')
    [m1,~,b] = otfunction(thet1, M, b,epsilon,gamma,tol_ot,max_iter_ot);
end
end

function thet_1=update_thet(thet_1,thet_2,Y,X,mu ,gamma,lamb,m,n_chanel,n_subject,n_source)
mu_gamma=mu * gamma /n_subject;
w = (mu_gamma +lamb) *ones(1,n_source);
X_1=X;
X_thet_Y=X_1*(thet_2)+Y;
c=-mu_gamma*m;
Xthet=X_1*thet_1;
for i=1:n_source
    X_thet=Xthet-X_1(:, i)*thet_1(i,:);
    b=w(:,i)+X_1(:, i)'*(X_thet-X_thet_Y)/n_chanel;
    a=norm(X_1(:, i),2)^2/n_chanel;
    thet_1(i,:)=(-b+sqrt(b^2-4*a*c(:,i)))/(2*a);
    Xthet=X_thet+X_1(:, i)*thet_1(i,:);
end
end
