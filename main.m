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
MTWVSSI_s_group=zeros(n_sub,ds,dt);
MTWVSSI_ev=zeros(1,n_sub);
MTWVSSI_thr=zeros(1,n_sub);
MTWVSSI_m=zeros(10,n_sub,ds);

%%  SVD for TBFs
Y_all2=EEGdata;
[Dic] = TBFSelection(Y_all2,0,'threshold','Permutation');
Y_DIC=Y_all2*Dic';
n_dic=size(Dic,1);
[dc,ds] = size(Gain);

M = D/median(reshape(D,1,[]));
epsilon=100. /ds;
gamma=max(max(M));
M=(- M / epsilon);
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
