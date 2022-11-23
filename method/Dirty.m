function [S] = Dirty(Y,X, alpha,beta)
[dc,ds]=size(X);
n_sub=size(Y,1);
S_c=zeros(ds,n_sub);
S_s=zeros(ds,n_sub);
%Y=Y1;
max_iter=1000;
tol=1e-4;
S=S_c+S_s;
L_S=zeros(n_sub,dc);
R=Y-L_S;
X_ALL=zeros(n_sub,dc,ds);
for xx=1:n_sub
    X_ALL(xx,:,:)=X;
end
X=X_ALL;
clear X_ALL;
LS=max(squeeze(sum(X.*X,2)),[],1)';
idx=find(LS==0);
LS(idx)=(min(LS(find(LS>0))));
alpha=alpha*dc;
beta=beta*dc;
for i=1:max_iter
    w_max=0;
    d_w_max=0;
    for j=1:ds
        if LS(j,:)==0
            continue;
        end
        grad=zeros(1,n_sub);
        tmp1=zeros(1,n_sub);
        tmp2=zeros(1,n_sub);
        
        normtmp=0;
        for t=1:n_sub
            for n=1:dc
                grad(:,t)=grad(:,t)+X(t,n,j)*R(t,n);
            end
            grad(:,t)=grad(:,t)./LS(j,:);
            tmp1(:,t)=grad(:,t)+S_c(j,t);
            tmp2(:,t)=grad(:,t)+S_s(j,t);
            normtmp=normtmp+tmp1(:,t).*tmp1(:,t);
        end
        normtmp=sqrt(normtmp);
        thresholdl2=0;
        if normtmp
            thresholdl2=max(1-alpha/(LS(j,:)*normtmp),0);
        end
        tmp1=tmp1*thresholdl2;
        thresholdl1=beta/LS(j,:);
        a=abs(tmp2)-thresholdl1;
        a(find(a<0))=0;
        tmp2=sign(tmp2).*a;
        new_S=tmp1+tmp2;
        if any(S(j,:))
            for t=1:n_sub
                R(t,:)=R(t,:)+squeeze(X(t,:,j)).*S(j,t);
            end
        end
        d_w_j=max(abs(S(j,:)-new_S));
        d_w_max=max(d_w_max,d_w_j);
        w_max=max(w_max,max(abs(tmp1+tmp2)));
        S_c(j,:)=tmp1;
        S_s(j,:)=tmp2;
        S(j,:)=new_S;
        if any(S(j,:))
            for t=1:n_sub
                R(t,:)=R(t,:)-squeeze(X(t,:,j)).*S(j,t);
            end
        end
    end
    if (w_max == 0 || d_w_max / w_max < tol)
            break
    end
end
end
