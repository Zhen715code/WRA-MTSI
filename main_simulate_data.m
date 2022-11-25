clear;
clc;
load_path='';
save_path='';
%% % % ----------------------------62*6002--------------------------------%%
load([load_path,'Cortex.mat']);% Cortex
load([load_path,'Gain.mat']);% Lead-filed matrix
load([load_path,'Gain3.mat']);
load([load_path,'GridLoc.mat']);
load([load_path,'Distance_6003.mat'],'D');
load([load_path,'seed_all.mat']);%Destrieux分区
load([load_path,'seed_simulate_all.mat']);
VertConn1 = VariationEdge(Cortex.VertConn);
n_subject=size(sub_list,2);
ds=size(Cortex.Vertices,1);
%%
dc=62;
dt=500;
SNIR=0;
SNR=0;
overlap=0.7;
n_seed = 3; % number of patches
area1=10;
%% 指标
n=50;%蒙特卡洛次数
% MTWVSSIdata.MTWVSSI_s=zeros(n,n_subject,ds,dt);
% save([save_path,'MTWVSSIdata'],'MTWVSSIdata');
% clear MTWVSSIdata;

MAX_ITER=15;%mtwvssi外层迭代次数
MTWVSSI_auc=zeros(1,n_subject);
MTWVSSI_sd=zeros(1,n_subject);
MTWVSSI_dle=zeros(1,n_subject);
MTWVSSI_rmse=zeros(1,n_subject);
MTWVSSI_se=zeros(1,n_subject);
%%
for k=1:n
    fprintf('%4s%4d\t%f\n','step--------------------------:',k);
    Real_s_group=zeros(n_subject,ds,dt);
    wMNE_s_group=zeros(n_subject,ds,dt);
    LORETA_s_group=zeros(n_subject,ds,dt);
    L2L2L2_s_group=zeros(n_subject,ds,dt);
    GL_s_group=zeros(n_subject,ds,dt);
    Dirty_s_group=zeros(n_subject,ds,dt);
    MTW_s_group=zeros(n_subject,ds,dt);
    MTWVSSI_s_group=zeros(n_subject,ds,dt);
    Y_all=zeros(n_subject,dc,dt);
    Y_all2=zeros(n_subject*dc,dt);
    L_all=zeros(n_subject,dc,ds);
    seedvox_all=squeeze(seed_simulate_all(k,:,:));
    ActiveVoxSeed_all=cell(n_subject,n_seed);
    [L_all(1,:,:),B1,Real_s_group(1,:,:),ActiveVoxSeed_all(1,:),StimTime]=Simulate_data2(Cortex,Gain,SNR,area1,seedvox_all(1,:),SNIR);
    Y_all(1,:,:)=B1;
    Y_all2(1:dc,:)=B1;
    for i=2:n_subject
        area2=10;
        [L,B,s_real,ActiveVoxSeed,~]=Simulate_data3(Cortex,Gain,SNR,area1,area2,n_seed,ActiveVoxSeed_all(1,:),seedvox_all(i,:),SNIR,overlap);
        Y_all(i,:,:)=B;
        Y_all2((i-1)*dc+1:i*dc,:)=B;
        Real_s_group(i,:,:)=s_real;
        L_all(i,:,:)=L;
        ActiveVoxSeed_all(i,:)=ActiveVoxSeed;
    end
    %load([save_path,'Real_all.mat'],'Real_all');
    %Real_all(k,:,:,:)= Real_s_group;
    %save([save_path,'Real_all'],'Real_all');
    %clear Real_all;
    %%  SVD for TBFs
    [Dic] = TBFSelection(Y_all2,0,'threshold','Permutation');
    Y_DIC=Y_all2*Dic';
    n_dic=size(Dic,1);
    clear Y_all2 B1 B L seedvox ActiveVoxSeed overlap_list s_real seedvox newvertsSeed1;
   
    M = D/median(reshape(D,1,[]));
    epsilon=100. /ds;
    gamma=max(max(M));
    M=(- M / epsilon);
    %%  WRA_MTSI
    fprintf('%s\n','-----------MTWVSSI--------------');
    M = D/median(reshape(D,1,[]));
    epsilon=100. /ds;
    gamma=max(max(M));
    M=(- M / epsilon);
    %超参
    MTWVSSI_lambda1=1e-2;
    MTWVSSI_lambda2=1e-2;
    MTWVSSI_rho=1;
    MTWVSSI_mu=1e-4;
    
    m= zeros(n_subject, ds) / ds;
    s=permute(theta,[2 1 3]);
    s_m=zeros(ds,n_subject);
    s_old=s;
    max_iter_ot = 20;
    tol_ot = 1e-4;
    for j=1:n_subject
        s_m(:,j)=sqrt(sum(s(:,j).*s(:,j),2));
    end
    if all(any(s_m,1),1)
        fprintf('m--------\n')
        m = otfunction(s_m, M, epsilon,gamma,tol_ot,max_iter_ot);
    end
    
    for i=1:30
        fprintf('%4s%4d\t%f\n','step--------------------------:',k);
        fprintf('%4s%4d\t%f\n','step:',i);
        fprintf('\n');
        [s,m]=WRA_MTSI(L_all, Y_DIC, VertConn1, M, gamma,epsilon,MTWVSSI_lambda1, MTWVSSI_lambda2, MTWVSSI_rho,MTWVSSI_mu,eps,eps,s,m,1.8);
        dx2=(norm(s(:,:,1) - s_old(:,1), 1)^ 2)/ (norm(s(:,:,1), 1) ^2);
        for kk=2:n_subject
            dx2=max(dx2,(norm(s(:,:,kk) - s_old(:,:,kk), 1)^ 2)/ (norm(s(:,:,kk), 1) ^2));
        end
        s_old = s;
        fprintf('%s%4d\n','s_dx2:',dx2)
        for j=1:n_subject
            s_mtwvsi1=reshape(s(:,:,j),ds,[])*Dic;
            s_real=reshape(Real_s_group(j,:,:),ds,[]);
            Roc =ROCextent(s_real(:,StimTime:end),s_mtwvsi1(:,StimTime:end),Cortex.VertConn,seedvox_all(j,:));
            MTWVSSI_auc(1,j) = median(Roc.mean);
            [MTWVSSI_sd(1,j),MTWVSSI_dle(1,j),MTWVSSI_rmse(1,j),MTWVSSI_se(1,j)]= PerformanceMetric(GridLoc,s_mtwvsi1(:,StimTime:end),s_real(:,StimTime:end), ActiveVoxSeed_all(j,:),'interval');
            MTWVSSI_s_group(j,:,:)=s_mtwvsi1;
            %MTWVSSI_ev(1,j) =  1 - norm(squeeze(Y_all(j,:,:)) - squeeze(L_all(j,:,:))*s_mtwvsi1,'fro')^2/norm(squeeze(Y_all(j,:,:)),'fro')^2;
        end
        fprintf('%f\t%f\t%f\t%f\t%f\n',mean(MTWVSSI_auc),mean(MTWVSSI_sd),mean(MTWVSSI_dle),mean(MTWVSSI_rmse),mean(MTWVSSI_se));
        
    end

    load([save_path,'MTWVSSIdata'],'MTWVSSIdata');
    MTWVSSIdata.MTWVSSI_auc(k,:)=MTWVSSI_auc;
    MTWVSSIdata.MTWVSSI_sd(k,:)=MTWVSSI_sd;
    MTWVSSIdata.MTWVSSI_dle(k,:)=MTWVSSI_dle;
    MTWVSSIdata.MTWVSSI_rmse(k,:)=MTWVSSI_rmse;
    %MTWVSSIdata.MTWVSSI_s(k,:,:,:)=MTWVSSI_s_group;
    MTWVSSIdata.MTWVSSI_se(k,:)=MTWVSSI_se;
    save([save_path,'MTWVSSIdata'],'MTWVSSIdata');
    %clear MTWVSSI_s_group MTWVSSIdata;
    if dx2<1e-5
        break;
    end
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
