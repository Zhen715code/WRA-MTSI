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
wMNE_auc=zeros(1,n_subject);
wMNE_sd=zeros(1,n_subject);
wMNE_dle=zeros(1,n_subject);
wMNE_rmse=zeros(1,n_subject);
wMNE_ev=zeros(1,n_subject);

% Real_all=zeros(n,n_subject,ds,dt);
% save([save_path,'Real_all'],'Real_all');
% clear Real_all;
% wMNEdata.wMNE_s=zeros(n,n_subject,ds,dt);
% save([save_path,'wMNEdata'],'wMNEdata');
% clear wMNEdata;
% LORETAdata.LORETA_s=zeros(n,n_subject,ds,dt);
% save([save_path,'LORETAdata'],'LORETAdata');
% clear LORETAdata;
% L2L2L2data.L2L2L2_s=zeros(n,n_subject,ds,dt);
% save([save_path,'L2L2L2data'],'L2L2L2data');
% clear L2L2L2data;
% GLdata.GL_s=zeros(n,n_subject,ds,dt);
% save([save_path,'GLdata'],'GLdata');
% clear GLdata;
% Dirtydata.Dirty_s=zeros(n,n_subject,ds,dt);
% save([save_path,'Dirtydata'],'Dirtydata');
% clear Dirty_s;
% MTWdata.MTW_s=zeros(n,n_subject,ds,dt);
% save([save_path,'MTWdata'],'MTWdata');
% clear MTWdata;
% MTWVSSIdata.MTWVSSI_s=zeros(n,n_subject,ds,dt);
% save([save_path,'MTWVSSIdata'],'MTWVSSIdata');
% clear MTWVSSIdata;

LORETA_auc=zeros(1,n_subject);
LORETA_sd=zeros(1,n_subject);
LORETA_dle=zeros(1,n_subject);
LORETA_rmse=zeros(1,n_subject);
LORETA_se=zeros(1,n_subject);

L2L2L2_auc=zeros(1,n_subject);
L2L2L2_sd=zeros(1,n_subject);
L2L2L2_dle=zeros(1,n_subject);
L2L2L2_rmse=zeros(1,n_subject);
L2L2L2_se=zeros(1,n_subject);

GL_auc=zeros(1,n_subject);
GL_sd=zeros(1,n_subject);
GL_dle=zeros(1,n_subject);
GL_rmse=zeros(1,n_subject);
GL_se=zeros(1,n_subject);

Dirty_auc=zeros(1,n_subject);
Dirty_sd=zeros(1,n_subject);
Dirty_dle=zeros(1,n_subject);
Dirty_rmse=zeros(1,n_subject);
Dirty_se=zeros(1,n_subject);

MTW_auc=zeros(1,n_subject);
MTW_sd=zeros(1,n_subject);
MTW_dle=zeros(1,n_subject);
MTW_rmse=zeros(1,n_subject);
MTW_se=zeros(1,n_subject);

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
    for sub=1:n_subject
        Y=squeeze(Y_all(sub,:,:));
        L=squeeze(L_all(sub,:,:));
        s_real=squeeze(Real_s_group(sub,:,:));
        ActiveVoxSeed=ActiveVoxSeed_all(sub,:);
        seedvox=seedvox_all(sub,:);
        %%  wMNE
        [Kernel,par] = MNE(Y,Gain3,L,Cortex,'MNE','reg',1);
        S_wMNE = Kernel*Y;
        [wMNE_sd(1,sub),wMNE_dle(1,sub),wMNE_rmse(1,sub)]= PerformanceMetric(GridLoc,S_wMNE(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval');
        Roc = ROCextent(s_real(:,StimTime:end),S_wMNE(:,StimTime:end),Cortex.VertConn,seedvox);
        wMNE_auc(1,sub) = median(Roc.mean);
        wMNE_s_group(sub,:,:)=S_wMNE;
        wMNE_ev(1,sub) =  1 - norm(Y - L*S_wMNE,'fro')^2/norm(Y,'fro')^2;
        %%  LORETA
        [Kernel,par] = MNE(Y,Gain3,L,Cortex,'LORETA','reg',1);
        S_LORETA = Kernel*Y;
        [LORETA_sd(1,sub),LORETA_dle(1,sub),LORETA_rmse(1,sub)]= PerformanceMetric(GridLoc,S_LORETA(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval');
        Roc = ROCextent(s_real(:,StimTime:end),S_LORETA(:,StimTime:end),Cortex.VertConn,seedvox);
        LORETA_auc(1,sub) = median(Roc.mean);
        LORETA_s_group(sub,:,:)=S_LORETA;
        LORETA_ev(1,sub) =  1 - norm(Y - L*S_LORETA,'fro')^2/norm(Y,'fro')^2;
        %%  SISSY_L21
        [Dic1] = TBFSelection(Y,0,'threshold','Permutation');
        [S_L2L2L2,R] = SISSY_L21(L, Y*Dic1',VertConn1,1e-1,1e-2,1e-1,1,0.05);%80,0.1  (L, b, v,  lambda1,lambda2,rho1,rho2,alpha)
        S_L2L2L2=S_L2L2L2*Dic1;
        [L2L2L2_sd(1,sub),L2L2L2_dle(1,sub),L2L2L2_rmse(1,sub)] ...
            = PerformanceMetric(GridLoc,S_L2L2L2(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed,'interval');
        Roc = ROCextent(s_real(:,StimTime:end),S_L2L2L2(:,StimTime:end),Cortex.VertConn,seedvox);
        L2L2L2_auc(1,sub) = median(Roc.mean);
        L2L2L2_s_group(sub,:,:)=S_L2L2L2;
        L2L2L2_ev(1,sub) =  1 - norm(Y - L*S_L2L2L2,'fro')^2/norm(Y,'fro')^2;
    end
    clear Y L s_real ActiveVoxSeed seedvox Kernel par S_LORETA Roc Dic1 S_L2L2L2 S_wMNE R;
    fprintf('%s\n','-----------GL-------------');
    S_GL=zeros(n_dic,ds,n_subject);
    lambda=1e-2;
    rho=1e2;
    for t=1:n_dic
        Y2=zeros(dc,n_subject);
        for sub=1:n_subject
            Y1=Y_DIC(:,t);
            Y2(:,sub)=Y1((sub-1)*dc+1:sub*dc,:);
            L=squeeze(L_all(sub,:,:));
        end
        S_GL(t,:,:) =GL(L, Y2,lambda,rho,1.8);
        %S_GL(t,:,:) = dirty(Y_DIC(:,t),L_all, alpha,beta);
    end
    for sub=1:n_subject
        s_GL=S_GL(:,:,sub)'*Dic;
        s_real=squeeze(Real_s_group(sub,:,:));
        Roc =ROCextent(s_real(:,StimTime:end),s_GL(:,StimTime:end),Cortex.VertConn,seedvox_all(sub,1),'flag',2);
        GL_auc(1,sub) = median(Roc.mean);
        [GL_sd(1,sub),GL_dle(1,sub),GL_rmse(1,sub),GL_se(1,sub)]= PerformanceMetric(GridLoc,s_GL(:,StimTime:end),s_real(:,StimTime:end), ActiveVoxSeed_all(sub,:));
        GL_s_group(sub,:,:)=s_GL;
        %GL_ev(1,sub) =  1 - norm(squeeze(Y_all(sub,:,:)) - squeeze(L_all(sub,:,:))*s_GL,'fro')^2/norm(squeeze(Y_all(sub,:,:)),'fro')^2;
    end
    clear alpha beta s_GL S_GL s_real Roc;
    
    %%  Dirty
    fprintf('%s\n','-----------Dirty--------------');
    S_Dirty=zeros(n_dic,ds,n_subject);
    alpha=0.3;
    beta=0.2;
    for t=1:n_dic
        S_Dirty(t,:,:) = Dirty(Y_DIC(:,t),L_all, alpha,beta);
    end
    for sub=1:n_subject
        s_dirty=S_Dirty(:,:,sub)'*Dic;
        s_real=squeeze(Real_s_group(sub,:,:));
        Roc =ROCextent(s_real(:,StimTime:end),s_dirty(:,StimTime:end),Cortex.VertConn,seedvox_all(sub,1),'flag',2);
        Dirty_auc(1,sub) = median(Roc.mean);
        [Dirty_sd(1,sub),Dirty_dle(1,sub),Dirty_rmse(1,sub),Dirty_se(1,sub)]= PerformanceMetric(GridLoc,s_dirty(:,StimTime:end),s_real(:,StimTime:end), ActiveVoxSeed_all(sub,:));
         Dirty_s_group(sub,:,:)=s_dirty;
         %Dirty_ev(1,sub) =  1 - norm(squeeze(Y_all(sub,:,:)) - squeeze(L_all(sub,:,:))*s_dirty,'fro')^2/norm(squeeze(Y_all(sub,:,:)),'fro')^2;
    end
    clear alpha beta s_dirty s_real S_Dirty Roc;
    
%     fprintf('%s\n','-----------wMNE-------------');
%     fprintf('%f\t%f\t%f\t%f\t%f\n',mean(wMNE_auc),mean(wMNE_sd),mean(wMNE_dle),mean(wMNE_rmse),mean(wMNE_ev));
%     fprintf('%s\n','-----------LORETA-------------');
%     fprintf('%f\t%f\t%f\t%f\t%f\n',mean(LORETA_auc),mean(LORETA_sd),mean(LORETA_dle),mean(LORETA_rmse),mean(LORETA_ev));
%     fprintf('%s\n','-----------L2L2L2-------------');
%     fprintf('%f\t%f\t%f\t%f\t%f\n',mean(L2L2L2_auc),mean(L2L2L2_sd),mean(L2L2L2_dle),mean(L2L2L2_rmse),mean(L2L2L2_ev));
    fprintf('%s\n','-----------GL-------------');
    fprintf('%f\t%f\t%f\t%f\t%f\n',mean(GL_auc),mean(GL_sd),mean(GL_dle),mean(GL_rmse),mean(GL_se));
    fprintf('%s\n','-----------Dirty------------');
    fprintf('%f\t%f\t%f\t%f\t%f\n',mean(Dirty_auc),mean(Dirty_sd),mean(Dirty_dle),mean(Dirty_rmse),mean(Dirty_se));
    
    %%  MTW
    fprintf('%s\n','-----------MTW--------------');
    M = D/median(reshape(D,1,[]));
    epsilon=100. /ds;
    gamma=max(max(M));
    M=(- M / epsilon);
    MTW_lamb=1e-1;
    MTW_mu=5;
    tol=5;
    MAX_ITER=200;
    theta=zeros(n_dic,ds,n_subject);
    thetaold=zeros(n_dic,ds,n_subject);
    thet1 =ones(n_dic,ds, n_subject)/ds ;
    thet2 =ones(n_dic,ds, n_subject)/ds ;
    m1 = ones(n_subject, ds)/ds;
    m2 = ones(n_subject, ds)/ds;
    for i=1:MAX_ITER
        fprintf('%4s%4d\t%f\n','step:',i);
        fprintf('\n');
        dx2=0;
        for t=1:n_dic
            [thet1(t,:,:),m1]=MTW( L_all, Y_DIC(:,t),M,thet1(t,:,:),thet2(t,:,:),epsilon, gamma,MTW_lamb, MTW_mu,m1);
            [thet2(t,:,:),m2]=MTW(L_all,-Y_DIC(:,t), M, thet2(t,:,:),thet1(t,:,:),epsilon, gamma,MTW_lamb, MTW_mu,m2);
            theta(t,:,:) = thet1(t,:,:) - thet2(t,:,:);
            dx2 =max(dx2,(norm(squeeze(theta(t,:,:)) - squeeze(thetaold(t,:,:)), 1)^ 2)/ (norm(squeeze(theta(t,:,:)), 1) ^2));
            thetaold(t,:,:)=theta(t,:,:);
        end
        fprintf('%s%4d\n','theta_dx2:',dx2);
        MTW_theta=theta;
        for sub=1:n_subject
            s_real=squeeze(Real_s_group(sub,:,:));
            s_mtw=MTW_theta(:,:,sub)'*Dic;
            Roc =ROCextent(s_real(:,StimTime:end),s_mtw(:,StimTime:end),Cortex.VertConn,seedvox_all(sub,1),'flag',2);
            MTW_auc(1,sub) = median(Roc.mean);
            [MTW_sd(1,sub),MTW_dle(1,sub),MTW_rmse(1,sub),MTW_se(1,sub)]= PerformanceMetric(GridLoc,s_mtw(:,StimTime:end),s_real(:,StimTime:end), ActiveVoxSeed_all(sub,:));
            MTW_s_group(sub,:,:)=s_mtw;
            %MTW_ev(1,sub) =  1 - norm(squeeze(Y_all(sub,:,:)) - squeeze(L_all(sub,:,:))*s_mtw,'fro')^2/norm(squeeze(Y_all(sub,:,:)),'fro')^2;
        end
        fprintf('%f\t%f\t%f\t%f\t%f\n',mean(MTW_auc),mean(MTW_sd),mean(MTW_dle),mean(MTW_rmse),mean(MTW_se));
        if dx2<1e-5
            break;
        end
    end
    clearvars MTW_lamb MTW_mu m1 m2 thetaold thet1 thet2 dx2 s_real s_mtw MTW_theta;
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
    clear MTWVSSI_lambda1  MTWVSSI_lambda2  MTWVSSI_rho  MTWVSSI_mu eps eps s m dx2 s_old s_real Roc s_mtwvsi1;
    
%     load([save_path,'wMNEdata'],'wMNEdata');
%     wMNEdata.wMNE_auc(k,:)=wMNE_auc;
%     wMNEdata.wMNE_sd(k,:)=wMNE_sd;
%     wMNEdata.wMNE_dle(k,:)=wMNE_dle;
%     wMNEdata.wMNE_rmse(k,:)=wMNE_rmse;
%     wMNEdata.wMNE_s(k,:,:,:)=wMNE_s_group;
%     wMNEdata.wMNE_ev(k,:)=wMNE_ev;
%     save([save_path,'wMNEdata'],'wMNEdata');
%     clear wMNE_s_group wMNEdata;
%     
%     load([save_path,'LORETAdata'],'LORETAdata');
%     LORETAdata.LORETA_auc(k,:)=LORETA_auc;
%     LORETAdata.LORETA_sd(k,:)=LORETA_sd;
%     LORETAdata.LORETA_dle(k,:)=LORETA_dle;
%     LORETAdata.LORETA_rmse(k,:)=LORETA_rmse;
%     LORETAdata.LORETA_ev(k,:)=LORETA_ev;
%     LORETAdata.LORETA_s(k,:,:,:)=LORETA_s_group;
%     save([save_path,'LORETAdata'],'LORETAdata');
%     clear LORETA_s_group LORETAdata;
%     
%     load([save_path,'L2L2L2data'],'L2L2L2data');
%     L2L2L2data.L2L2L2_s(k,:,:,:)=L2L2L2_s_group;
%     L2L2L2data.L2L2L2_auc(k,:)=L2L2L2_auc;
%     L2L2L2data.L2L2L2_sd(k,:)=L2L2L2_sd;
%     L2L2L2data.L2L2L2_dle(k,:)=L2L2L2_dle;
%     L2L2L2data.L2L2L2_rmse(k,:)=L2L2L2_rmse;
%     L2L2L2data.L2L2L2_ev(k,:)=L2L2L2_ev;
%     save([save_path,'L2L2L2data'],'L2L2L2data');
%     clear L2L2L2_s_group L2L2L2data;
    
    load([save_path,'GLdata'],'GLdata');
    GLdata.GL_auc(k,:)=GL_auc;
    GLdata.GL_sd(k,:)=GL_sd;
    GLdata.GL_dle(k,:)=GL_dle;
    GLdata.GL_rmse(k,:)=GL_rmse;
    %GLdata.GL_s(k,:,:,:)=GL_s_group;
    GLdata.GL_se(k,:)=GL_se;
    save([save_path,'GLdata'],'GLdata');
    %clear GL_s_group GLdata;
    
    load([save_path,'Dirtydata'],'Dirtydata');
    Dirtydata.Dirty_auc(k,:)=Dirty_auc;
    Dirtydata.Dirty_sd(k,:)=Dirty_sd;
    Dirtydata.Dirty_dle(k,:)=Dirty_dle;
    Dirtydata.Dirty_rmse(k,:)=Dirty_rmse;
    %Dirtydata.Dirty_s(k,:,:,:)=Dirty_s_group;
    Dirtydata.Dirty_se(k,:)=Dirty_se;
    save([save_path,'Dirtydata'],'Dirtydata');
    %clear Dirty_s_group Dirtydata;
    
    load([save_path,'MTWdata'],'MTWdata');
    MTWdata.MTW_auc(k,:)=MTW_auc;
    MTWdata.MTW_sd(k,:)=MTW_sd;
    MTWdata.MTW_dle(k,:)=MTW_dle;
    MTWdata.MTW_rmse(k,:)=MTW_rmse;
    %MTWdata.MTW_s(k,:,:,:)=MTW_s_group;
    MTWdata.MTW_se(k,:)=MTW_se;
    save([save_path,'MTWdata'],'MTWdata');
    %clear MTW_s_group MTWdata;
    
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
