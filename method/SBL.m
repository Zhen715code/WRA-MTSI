function [ImagingKernel,par]=SBL(B,L,varargin)
%==========================================================================
% Discriptions: this function is to compute the unknown source s.The model
%    is B=Ls+noise. 
% Input: L:   the leadfield matrix
%        B:   the observed data

%        varargin:%!!!!就是'flags',0,'prune',[1 1e-6],'epsilon',1e-4, 'max_iter',150,'print',1,'nlf',1！！！！！
%       flags:  choose which algorithm to compute the unknown s
%                      flags=0:  EM algorithm
%                      flags=1:  Mackay updatas
%                      flags=2:  Convexity_based approach
%                      flags=[3,p]: SMAP (MFOCUSS)
%      prune (1x2 vector):   whether prune the smallest gamma during iterations (0 or 1) and the prune threshold ( typically 1e-6)
%      epsilon:  stop condition (typically 1e-4)
%      max_iter: the largest numbers of iteration  
%      print :   whether print the iteration process
%      nlf   :   whether normalize the lead field matrix
          
% Output: ImagingKernel
%        gamma : estimated hyperparameters
%        evidence : cost funcitons

% Reference: David Wipf and Srikantan Nagarajan
%         A unified Bayesian framework for MEG/EEG source imaging. 2009

% data: 2012/8/30 (version 1.1)
%       2013/12/2 (version 1.2)
% Author: Liu Ke
%==========================================================================
 
% Dimension of the inverse problem
 nSource=size(L,2);                 % number of total sources
 [nSensor,nSnap]=size(B);           % number of sensors and time points    B矩阵是db*n
 % Default Control Parameters 如果没有给初始值，就使用默认值
epsilon        = 1e-6;        % threshold for stopping iteration. 
MAX_ITER       = 300;        % maximum iterations
print          = 1;           % show progress information
flags          = 2;           % iteration algorithm (0:EM; 1:Makay; 2:Convexity; 3: MFOCUSS)
prune          = [0 1e-6];    % whether prune the smallest gamma第一个值表示是否删掉最小的gamma（0或者1）,第二个值表示删除阀值
NLFM           = 0;           % whether normalize the lead field matrix
beta           = 1;           % variance of sensor noise
Cov_n          = eye(nSensor);%db*db的单位矩阵
% get input argument values
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})%把大写字母转换成小写字母
            case 'nlf'
                NLFM = varargin{i+1}; 
            case 'epsilon'   
                epsilon = varargin{i+1}; 
            case 'print'    
                print = varargin{i+1}; 
            case 'max_iter'
                MAX_ITER = varargin{i+1}; 
            case 'cov_n'
                Cov_n = varargin{i+1}; 
            case 'flags'
                flags = varargin{i+1};  
            case 'prune'
                prune = varargin{i+1};  
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end

fprintf('\nRunning Sparse Bayesian Learining...\n'); 
 %---------------initial states--------------%
 gamma = (ones(nSource,1) + randn(nSource,1)*1e-4)*nSensor/trace(L*L');% ？？？？？？ds*1的1矩阵加上ds*1的均匀分布随机矩阵，乘以ds,除以L*L'迹??????initial values of hyperparameters
 cost_old = 1;
 keep_list = (1:nSource)';%[1;2;3;4;....;ds]列向量
 evidence = zeros(MAX_ITER,1);%每次迭代的成本函数
 
%======================================%
%         Leadfield Matrix normalization
%=====================================%
if (NLFM)
% iG=sqrt(sum(L.^2,1))';
% iG(iG == 0) = min(iG(iG ~= 0));
% iG = spdiags((1 ./ iG) .^ 0.5, 0, length(iG), length(iG)) ;
% L =L* iG ;
% clear iG
    LfvW = (mean(L.^2,1)).^0.5;%将L每个元素平方，对列求均值，再每个元素开根，与每列直接求均值近似,1*ds
    L = L.*kron(ones(nSensor,1),1./LfvW);%？？？？？ones(nSensor,1)是db*1的1矩阵，1./LfvW是1*ds的矩阵，kron(ones(nSensor,1),1./LfvW)是db*ds矩阵
end
%=========================================================================%
%               iteration
%=========================================================================%
 for iter=1:MAX_ITER
% ******** Prune things as hyperparameters go to zero **********%当超参数变为0时删除
if (prune(1))
    index = find(gamma > max(gamma)* prune(2));
    gamma = gamma(index);
    L = L(:,index);
    keep_list = keep_list(index);
    if (isempty(gamma))   
        break; 
    end
end
% ========== Update hyperparameters ===========%
Ns = numel(keep_list);%ds 的个数
Cov_b = beta*Cov_n + L.*repmat(gamma',nSensor,1)*L';%(16)gamma'是1*ds,扩展到db行，列不扩展，repmat(gamma',nSensor,1)就是cov_s,db*ds,注意，前面是点乘L.*cov_s,后面是乘*L'，cov_b是db*db
  
for k = 1:Ns%把1到ds个源位置的信息循环,然后求出gama（1），gama（2），...,gama（ds），也就是求出了gama向量，ds*1
                          
    M = gamma(k)*L(:,k)'/Cov_b;%第k个源位置的gama乘以L的所有行第k列的转置，即L的第k个源位置的所有传感器信息                                                                         
    if (flags(1) == 0)
        %-------------EM algorithm---------------%
        gamma(k) = (norm(M*B,'fro'))^2/nSnap+gamma(k)-trace(M*L(:,k)*gamma(k));%(27)
     elseif (flags(1) == 1)
              % -----------MacKay updatas--------------%
           
        gamma(k) = (norm(M*B,'fro'))^2/(nSnap*trace(M*L(:,k))); %(29)                                                                                                                                                             
    elseif (flags(1) == 2)
        %--------Convexity_based approach--------%
        gamma(k) = gamma(k)*norm(L(:,k)'/Cov_b*B,'fro')/sqrt(nSnap*trace((L(:,k)'/Cov_b)*L(:,k)));%(30)
    elseif (flags(1) == 3)
        %------------S-MAP-----------------------%
        p=flags(2);
        gamma(k) = ((norm(M*B,'fro'))^2/nSnap)^((2-p)/2);
    end
end
%  beta = norm(beta*Cov_b\B,'fro')^2/(nSnap*beta*trace(inv(Cov_b)))

cost = trace(B*B'/Cov_b)+nSnap*log(det(Cov_b))+nSnap*nSensor*log(2*pi);%(18)？？？？第三部分？？？？
cost = -0.5*cost;%？？？？？
MSE = (cost-cost_old)/cost;%？？？？？
cost_old = cost;
if (print)%展示过程
    disp(['SBLiters: ',num2str(iter),'   num voxels: ',num2str(Ns),'  MSE: ',num2str(MSE)])%输出这一次迭代的次数，ds未知源的数量，成本函数MSE
end
if (abs(MSE) < epsilon)  break; end%成本函数要足够小，都小于阀值了，就可以停止
evidence(iter)=cost;
 end
%% compute ImagingKernel and s
evidence=evidence(1:iter-1,:);

% Cov = diag(gamma)-diag(gamma)*L'/(eye(nSensor)+L*diag(gamma)*L')*L*diag(gamma);
% temp = zeros(nSource,1);
% temp(keep_list) = diag(Cov);
% Cov = temp;

temp=zeros(nSource,nSensor);%ds*db
ImagingKernel = (repmat(gamma,1,nSensor).*L')/Cov_b;%（13）repmat(gamma,1,nSensor)是cov_s,ds*db,repmat(gamma,1,nSensor).*L'是点乘，得到ds*db,ImagingKernel还是ds*db

if (Ns>0)
    temp(keep_list,:) = ImagingKernel;% 把ImagingKernel赋值给temp的第1到ds行，所有列
end
 ImagingKernel = temp;
 s=ImagingKernel*B;
 par.gamma = gamma;%把gama赋值给par.gama,把par构造成一个结构体
 par.evidence = evidence;
 par.beta = beta;
 par.keeplist = keep_list;
fprintf('\nSparse Bayesian Learining Finished!\n'); 

