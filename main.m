clear all;
VertConn1 = VariationEdge(Cortex.VertConn);
n_sub=6;
dt=500;
[dc,ds] = size(Gain);
Path = '';
File = dir(fullfile(Path, '*.mat'));
FileNames = {File.name};
file_length = length(FileNames);
%%
%======================================%
%         Leadfield Matrix normalization
%=====================================%
LfvW = (mean(Gain.^2,1)).^0.5;
Gain = Gain.*kron(ones(dc,1),1./LfvW);
%% Whiten measurements and lead field matrix
noice_Path = '';
[Cov_noise,W]=computed_noise(noice_Path);
dc=size(W,1);
EEGdata=zeros(n_sub*dc,dt);
L_all=zeros(n_sub,dc,ds);
ratio=1e-5;
%%
for index =1:n_sub
    file = char(FileNames(index));
    file_path = strcat(Path, file);
    disp(file_path);
    load(file_path);
    Y=F([1:32,34:42,44:64],2251:2750);
    Y = Y./ratio;
    Y = W*Y;
    EEGdata((index-1)*dc+1:index*dc,:)=Y;
    L_all(index,:,:)=W*Gain;
end
%%
wMNE_s_group=zeros(n_sub,ds,dt);
LORETA_s_group=zeros(n_sub,ds,dt);
L2L2L2_s_group=zeros(n_sub,ds,dt);
SBL_s_group=zeros(n_sub,ds,dt);
GL_s_group=zeros(n_sub,ds,dt);
Dirty_s_group=zeros(n_sub,ds,dt);
MTW_s_group=zeros(n_sub,ds,dt);
MTWVSSI_s_group=zeros(n_sub,ds,dt);
 
wMNE_ev=zeros(1,n_sub);
LORETA_ev=zeros(1,n_sub);
L2L2L2_ev=zeros(1,n_sub);
SBL_ev=zeros(1,n_sub);
GL_ev=zeros(1,n_sub);
Dirty_ev=zeros(1,n_sub);
MTW_ev=zeros(1,n_sub);
MTWVSSI_ev=zeros(1,n_sub);

wMNE_thr=zeros(1,n_sub);
LORETA_thr=zeros(1,n_sub);
L2L2L2_thr=zeros(1,n_sub);
SBL_thr=zeros(1,n_sub);
GL_thr=zeros(1,n_sub);
Dirty_thr=zeros(1,n_sub);
MTW_thr=zeros(1,n_sub);
MTWVSSI_thr=zeros(1,n_sub);
MTWVSSI_m=zeros(10,n_sub,ds);
sigma=zeros(1,n_sub);
for sub=1:n_sub
    Y=EEGdata((sub-1)*dc+1:sub*dc,:);
    Gain1=squeeze(L_all(sub,:,:));
%     %  wMNE
%     [Kernel,par] = MNE(Y,Gain3,Gain1,Cortex,'wMNE','reg',1);
%     S_wMNE = Kernel*Y*ratio;
%     wMNE_s_group(sub,:,:)=S_wMNE;
%     wMNE_thr(1,sub)=ThresholdSelect(S_wMNE, 'ifplot',0);
% %     S_wMNE=squeeze(wMNE_s_group(sub,:,:));
%     wMNE_ev(1,sub)= 1 - norm(Y - Gain1*S_wMNE./ratio,'fro')^2/norm(Y,'fro')^2;
%     
%     %  LORETA
%     [Kernel,par] = MNE(Y,Gain3,Gain1,Cortex,'LORETA','reg',1);
%     S_LORETA = Kernel*Y*ratio;
%     LORETA_s_group(sub,:,:)=S_LORETA;
%     LORETA_thr(1,sub)=ThresholdSelect(S_LORETA, 'ifplot',0);
%     %     S_LORETA=squeeze(LORETA_s_group(sub,:,:));
%     LORETA_ev(1,sub)=1 - norm(Y - Gain1*S_LORETA./ratio,'fro')^2/norm(Y,'fro')^2;
%     
%     %  SBL
%     [Kernel,par] = SBL(Y,Gain1,'epsilon',1e-4,'flags',1,'prune',[1,1e-6]);
%     S_SBL = Kernel*Y*ratio;
%     SBL_s_group(sub,:,:)=S_SBL;
%     SBL_thr(1,sub)=ThresholdSelect(S_SBL, 'ifplot',0);
%     %     S_SBL=squeeze(SBL_s_group(sub,:,:));
%     SBL_ev(1,sub)=1 - norm(Y - Gain1*S_SBL./ratio,'fro')^2/norm(Y,'fro')^2;
    
    %%   SISSY_L21
    [Dic1] = TBFSelection(Y,0,'threshold','Permutation');
    sigma=0.01*norm(Y*Dic1',1)/sqrt(dc);
    lamb=1;
    mu=1e-2; 
    [S_L2L2L2,R] =  SISSY_L21(Gain1, Y*Dic1',VertConn1,lamb*mu,lamb,1e1,1e3,1.5);%80,0.1  (L, b, v,  lambda1,lambda2,rho1,rho2,alpha)
    S_L2L2L2=S_L2L2L2*Dic1*ratio;
    L2L2L2_s_group(sub,:,:)=S_L2L2L2;
    L2L2L2_thr(1,sub)=ThresholdSelect(S_L2L2L2, 'ifplot',0);
    %     S_L2L2L2=squeeze(L2L2L2_s_group(sub,:,:))*ratio;
     L2L2L2_ev(1,sub)=1 - norm(Y - Gain1*S_L2L2L2./ratio,'fro')^2/norm(Y,'fro')^2;
    
end
%%  SVD for TBFs
Y_all2=EEGdata;
[Dic] = TBFSelection(Y_all2,0,'threshold','Permutation');
Y_DIC=Y_all2*Dic';
n_dic=size(Dic,1);
[dc,ds] = size(Gain);
%%  GL 
fprintf('%s\n','-----------GL-------------');
S_GL=zeros(n_dic,ds,n_sub);
lambda=1e-1;
rho=1e1;
for t=1:n_dic
    Y2=zeros(dc,n_sub);
    for sub=1:n_sub
        Y1=Y_DIC(:,t);
        Y2(:,sub)=Y1((sub-1)*dc+1:sub*dc,:);
        L=squeeze(L_all(sub,:,:));
    end
    S_GL(t,:,:) =GL(L, Y2,lambda,rho,1.8);
end
for sub=1:n_sub
    s_GL=S_GL(:,:,sub)'*Dic*ratio;
    GL_s_group(sub,:,:)=s_GL;
    GL_thr(1,sub)=ThresholdSelect(s_GL, 'ifplot',0);
    
    Y=EEGdata((sub-1)*dc+1:sub*dc,:)*ratio;
    GL_ev(1,sub)=1 - norm(Y - Gain1*s_GL,'fro')^2/norm(Y,'fro')^2;
end

%%  Dirty
fprintf('%s\n','-----------Dirty--------------');
S_Dirty=zeros(n_dic,ds,n_sub);
alpha=10;
beta=8;
for t=1:n_dic
    S_Dirty(t,:,:) = Dirty(Y_DIC(:,t),L_all, alpha,beta);
end
for sub=1:n_sub
    s_dirty=S_Dirty(:,:,sub)'*Dic*ratio;
    Dirty_s_group(sub,:,:)=s_dirty;
    Dirty_thr(1,sub)=ThresholdSelect(s_dirty, 'ifplot',0);
    Y=EEGdata((sub-1)*dc+1:sub*dc,:)*ratio;
   Dirty_ev(1,sub)=1 - norm(Y - Gain1*s_dirty,'fro')^2/norm(Y,'fro')^2;
end

%%  MTW
fprintf('%s\n','-----------MTW--------------');
M = D/median(reshape(D,1,[]));
epsilon=100. /ds;
gamma=max(max(M));
M=(- M / epsilon);
MTW_lamb=1e1;
MTW_mu=1e1;
tol=5;
MAX_ITER=100;
theta=zeros(n_dic,ds,n_sub);
thetaold=zeros(n_dic,ds,n_sub);
thet1 =ones(n_dic,ds, n_sub)/ds ;
thet2 =ones(n_dic,ds, n_sub)/ds ;
m1 = ones(n_sub, ds)/ds;
m2 = ones(n_sub, ds)/ds;
b1=[];
b2=[];
for i=1:MAX_ITER
    fprintf('%4s%4d\t%f\n','step:',i);
    fprintf('\n');
    dx2=0;
    for t=1:n_dic
        [thet1(t,:,:),m1,b1]=MTW( L_all, Y_DIC(:,t),M,thet1(t,:,:),thet2(t,:,:),b1,epsilon, gamma,MTW_lamb, MTW_mu,m1);
        [thet2(t,:,:),m2,b2]=MTW(L_all,-Y_DIC(:,t), M, thet2(t,:,:),thet1(t,:,:),b2,epsilon, gamma,MTW_lamb, MTW_mu,m2);
        theta(t,:,:) = thet1(t,:,:) - thet2(t,:,:);
        dx2 =max(dx2,(norm(squeeze(theta(t,:,:)) - squeeze(thetaold(t,:,:)), 1)^ 2)/ (norm(squeeze(theta(t,:,:)), 1) ^2));
        thetaold(t,:,:)=theta(t,:,:);
    end
    fprintf('%s%4d\n','theta_dx2:',dx2);
    MTW_theta=theta;
    for sub=1:n_sub
        s_mtw=MTW_theta(:,:,sub)'*Dic*ratio;
        MTW_s_group(sub,:,:)=s_mtw;
        MTW_thr(1,sub)=ThresholdSelect(s_mtw, 'ifplot',0);
        Y=EEGdata((sub-1)*dc+1:sub*dc,:)*ratio;
        MTW_ev(1,sub)=1 - norm(Y - Gain1*s_mtw,'fro')^2/norm(Y,'fro')^2;
    end
     if dx2<1e-4
        break;
    end
end

%%  WRA_MTSI
MTWVSSI_barycenter=zeros(ds,15);
sigma0=zeros(1,n_sub);
b=[];
fprintf('%s\n','-----------MTWVSSI--------------');
%³¬²Î
MTWVSSI_lambda1=1e-3;
MTWVSSI_lambda2=5;
MTWVSSI_rho=1e5;
MTWVSSI_mu=1e-3;

m= zeros(n_sub, ds);
s=ones(ds,n_dic,n_sub)/1e3;
wes_dis=zeros(15,n_sub);
s_old=s;
max_iter_ot = 20;
tol_ot = 1e-4;
s_m=zeros(15,n_sub,ds);
for j=1:n_sub
    s1=squeeze(Dirty_s_group(j,:,:))*Dic'/ratio;
    s(:,:,j)=s1;
    sigma0(1,j)=0.01*norm(Y_DIC((j-1)*dc+1:j*dc,:),1)/sqrt(dc);
end
s_old=s;
sigma=sigma0;
% if all(any(s_m,1),1)
%     fprintf('m--------\n')
%     [m,q,b] = otfunction(s_m, M,b, epsilon,gamma,tol_ot,max_iter_ot);
% end

Gain1=squeeze(L_all(sub,:,:));
for i=1:30
    fprintf('%4s%4d\t%f\n','step:',i);
    fprintf('\n');
    [s,m,q,b,d]=WRA_MTSI(Gain1, Y_DIC,VertConn1, M,b, gamma,epsilon,MTWVSSI_lambda1, MTWVSSI_lambda2, MTWVSSI_rho,MTWVSSI_mu,eps,eps,s,m,1.6,n_sub); 
    MTWVSSI_barycenter(:,i)=q;
    MTWVSSI_m(i,:,:)=m;
    dx2=(norm(s(:,:,1) - s_old(:,1), 1)^ 2)/ (norm(s(:,:,1), 1) ^2);
    for kk=2:n_sub
        dx2=max(dx2,(norm(s(:,:,kk) - s_old(:,:,kk), 1)^ 2)/ (norm(s(:,:,kk), 1) ^2));
    end
    s_old = s;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    fprintf('%s%4d\n','s_dx2:',dx2)
    wes_dis(i,:)=d;
    for j=1:n_sub
        s_mtwvsi1=reshape(s(:,:,j),ds,[])*Dic*ratio;
        MTWVSSI_s_group(j,:,:)=s_mtwvsi1;
        MTWVSSI_thr(1,j)=ThresholdSelect(s_mtwvsi1, 'ifplot',0);
        Y=EEGdata((j-1)*dc+1:j*dc,:)*ratio;
        MTWVSSI_ev(1,j)=1 - norm(Y - Gain1*s_mtwvsi1,'fro')^2/norm(Y,'fro')^2;
        s_m(i,j,:)=sqrt(sum(s_mtwvsi1.*s_mtwvsi1,2))/ratio;
    end
    if dx2<1e-5
        break;
    end
end
