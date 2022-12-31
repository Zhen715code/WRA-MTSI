clear;
clc;
load_path='\data\';
save_path='\result\result_16t\';
%% % % ----------------------------62*6002--------------------------------%%
load([load_path,'Cortex.mat']);% Cortex
load([load_path,'Gain.mat']);% Lead-filed matrix
load([load_path,'Gain3.mat']);
load([load_path,'GridLoc.mat']);
load([load_path,'Distance_6003.mat'],'D');
% load('\seed_simulate_all_2t.mat');
% load('\seed_simulate_all_8t.mat');
load('\seed_simulate_all_16t.mat');
% load('\seed_simulate_all_overlap1.mat');
% load('\seed_simulate_all_overlap4.mat');
%load('\seed_simulate_all_overlap7_4.mat');
% load('\seed_simulate_all_overlap10.mat');
VertConn1 = VariationEdge(Cortex.VertConn);

ds=size(Cortex.Vertices,1);
%%
dc=62;
dt=500;
n_subject=4;
SNIR=0;
SNR=0;
overlap=0.7;
n_seed = 3; % number of patches
area=10;
n=50;

% Real_all=zeros(n,n_subject,ds,dt);
% save([save_path,'Real_all'],'Real_all');
% clear Real_all;
% WRAdata.WRA_s=zeros(n,n_subject,ds,dt);
% save([save_path,'WRAdata'],'WRAdata');
% clear WRAdata;

MAX_ITER=15;
WRA_auc=zeros(1,n_subject);
WRA_sd=zeros(1,n_subject);
WRA_dle=zeros(1,n_subject);
WRA_rmse=zeros(1,n_subject);
%%
for k=1:n
    fprintf('%4s%4d\t%f\n','step--------------------------:',k);
    
    Real_s_group=zeros(n_subject,ds,dt);
   
    WRA_s_group=zeros(n_subject,ds,dt);
    ActiveVoxSeed=squeeze(ActiveVoxSeed_all(k,:,:));
    seedvox=squeeze(seedvox_all(k,:,:));
    
    Y_all=zeros(n_subject,dc,dt);
    Y_all2=zeros(n_subject*dc,dt);
    L_all=zeros(n_subject,dc,ds);
    
    for i=1:n_subject
        [L,B,s_real]=Simulate_data(Cortex,Gain,SNR,ActiveVoxSeed(i,:),seedvox(i,:),SNIR);
        Y_all(i,:,:)=B;
        Y_all2((i-1)*dc+1:i*dc,:)=B;
        Real_s_group(i,:,:)=s_real;
        L_all(i,:,:)=L;
    end
    %     load([save_path,'Real_all.mat'],'Real_all');
    %     Real_all(k,:,:,:)= Real_s_group;
    %     save([save_path,'Real_all'],'Real_all');
    %     clear Real_all;
    %%  SVD for TBFs
    [Dic] = TBFSelection(Y_all2,0,'threshold','Permutation');
    Y_DIC=Y_all2*Dic';
    n_dic=size(Dic,1);
    StimTime=251;

    %%  WRA
    fprintf('%s\n','-----------WRA--------------');
    M = D/median(reshape(D,1,[]));
    epsilon=100. /ds;
    gamma=max(max(M));
    M=(- M / epsilon);
    WRA_lambda1=1e-2;
    WRA_lambda2=1e-2;
    WRA_rho=1;
    WRA_mu=1e-3;
    m= zeros(n_subject, ds) / ds;
    s=zeros(ds,n_dic,n_subject);
    s_m=zeros(ds,n_subject);
    s_old=s;
    max_iter_ot = 20;
    tol_ot = 1e-4;
    b=[];
    for j=1:n_subject
        s_m(:,j)=sqrt(sum(squeeze(s(:,:,j)).*squeeze(s(:,:,j)),2));
    end
    if all(any(s_m,1),1)
        fprintf('m--------\n')
        [m,~,b] = otfunction(s_m, M,b, epsilon,gamma,tol_ot,max_iter_ot);
    end
    
    for i=1:20
        fprintf('%4s%4d\t%f\n','step--------------------------:',k);
        fprintf('%4s%4d\t%f\n','step:',i);
        fprintf('\n');
        [s,m,~,b,d]=WRA_MTSI(L, Y_DIC, VertConn1, M,b, gamma,epsilon,WRA_lambda1, WRA_lambda2, WRA_rho,WRA_mu,eps,eps,s,m,1.8,n_subject);
        dx2=(norm(s(:,:,1) - s_old(:,1), 1)^ 2)/ (norm(s(:,:,1), 1) ^2);
        for kk=2:n_subject
            dx2=max(dx2,(norm(s(:,:,kk) - s_old(:,:,kk), 1)^ 2)/ (norm(s(:,:,kk), 1) ^2));
        end
        s_old = s;
        fprintf('%s%4d\n','s_dx2:',dx2)
        for j=1:n_subject
            s_wra1=reshape(s(:,:,j),ds,[])*Dic;
            s_real=reshape(Real_s_group(j,:,:),ds,[]);
            Roc =ROCextent(s_real(:,StimTime:end),s_wra1(:,StimTime:end),Cortex.VertConn,seedvox(j,:));
            WRA_auc(1,j) = median(Roc.mean);
            [WRA_sd(1,j),WRA_dle(1,j),WRA_rmse(1,j)]= PerformanceMetric(GridLoc,s_wra1(:,StimTime:end),s_real(:,StimTime:end), ActiveVoxSeed(j,:),'interval');
            WRA_s_group(j,:,:)=s_wra1;
        end
        fprintf('%f\t%f\t%f\t%f\t%f\n',mean(WRA_auc),mean(WRA_sd),mean(WRA_dle),mean(WRA_rmse));
        if dx2<1e-4
            break;
        end
    end
    clear WRA_lambda1  WRA_lambda2  WRA_rho  WRA_mu eps eps s m dx2 s_old s_real Roc s_wra1;
       
    %load([save_path,'WRAdata'],'WRAdata');
    WRAdata.WRA_auc(k,:)=WRA_auc;
    WRAdata.WRA_sd(k,:)=WRA_sd;
    WRAdata.WRA_dle(k,:)=WRA_dle;
    WRAdata.WRA_rmse(k,:)=WRA_rmse;
    %WRAdata.WRA_s(k,:,:,:)=WRA_s_group;
    save([save_path,'WRAdata'],'WRAdata');
    clear WRA_s_group;
end
