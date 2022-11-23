function [Kernel,par] = MNE(B,Gain3,L,Cortex,InverseMethod,varargin)
% Discription: This script is to compute the L2-norm inverse problem,
% including MNE, wMNE, sloreta, dspm and LORETA. The basis code is from
% Brainstorm.
% S = min_S {||B - LS||_F^2 + alpha*||WS||_2^2}
% alpha is estimate via maximizing the evidence p(B|alpha).


% Input: 
%       B:               Whitened M/EEG data
%       Gain3:            Three orientation leadfield matrix (is used to computing the weight of each source for wMNE and LORETA)
%       L:               Whitened single orientation leadfield matrix
%       VertConn:        Connective condition of each source (used to computing Laplacian Matrix in LORETA)
%       InverseMethod:   optional methods are 'MNE', 'wMNE', 'sloreta', 'LORETA', 'dspm'

% Output: 
%      Kernel:           Imaging Kernel of the inverse method. The estimate
%      of source is    S = Kernel*B;

% Author: Liu Ke
% Date:   2014/5/16
%% ===== DEFINE DEFAULT OPTIONS =====
[nSensor,nSource] = size(L);
nSnap = size(B,2);
SNR = 3;
pQ = [];
OPTIONS.fMRI        = [];
OPTIONS.fMRIthresh  = [];
OPTIONS.fMRI_prior  = [];
OPTIONS.depth = 1;
OPTIONS.weightlimit = 10;
OPTIONS.weightexp = 0.5;
Reg = 1;

for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'reg'   
                Reg = varargin{i+1}; 
            case 'fmri'   
                pQ = varargin{i+1}; 
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
end

if ~isempty(pQ)%pQÊÇ£¿
    supra_threshold = [];
    for i = 1:numel(pQ)
        supra_threshold = union(supra_threshold,find(pQ{i}~=0));
    end
    OPTIONS.fMRI_prior = 0.1*ones(nSource,1);
    OPTIONS.fMRI_prior(supra_threshold) = 1;
end
%% ===== Depth compensation =====
% Depth compensation
if any(strcmpi(InverseMethod, {'wmne','wMNE'}))
    if size(Gain3) == size(L)
        % if there is no three orientation lead filed matrix
        if isempty(OPTIONS.fMRI_prior)
            w = ones(nSource,1);
            %         w = 1./sum(Gain.^2,1)';
        else
            w = OPTIONS.fMRI_prior;
        end
    else
        if isempty(OPTIONS.fMRI_prior)
            szl = size(Gain3);
            w = squeeze(sum(sum((reshape(Gain3,[szl(1) 3 szl(2)/3])) .^2,1),2));
            w = 1 ./ w;
        else
            w = OPTIONS.fMRI_prior;
        end
        %% ===== APPLY DEPTH WEIGHTHING =====
        if OPTIONS.depth
            % ===== APPLY WEIGHT LIMIT =====
            % Applying weight limit.
            display('Applying weight limit.')
            weightlimit2 = OPTIONS.weightlimit .^ 2;
            %limit=min(w(w>min(w)*weightlimit2));  % This is the Matti way.
            limit = min(w) * weightlimit2;  % This is the Rey way (robust to possible weight discontinuity).
            w(w>limit) = limit;
            
            % ===== APPLY WEIGHT EXPONENT =====
            % Applying weight exponent.
            display('Applying weight exponent.')
            w = w .^ OPTIONS.weightexp;
            clear limit weightlimit2
            % elseif  ~isempty(OPTIONS.fMRI)
            %     display('wMNE> Applying fMRI priors.')
            %     ifmri = (OPTIONS.fMRI < OPTIONS.fMRIthresh);
            %     w(ifmri) = w(ifmri) * OPTIONS.fMRIoff;
            % end
        end
    end
elseif any(strcmpi(InverseMethod, {'mne','MNE','dspm','dSPM','sloreta','sLORETA','LORETA','loreta'}))
    w = ones(nSource,1);
end
%% ===== ADJUSTING SOURCE COVARIANCE MATRIX =====
% Adjusting Source Covariance matrix to make trace of L*C_J*L' equal to number of sensors.
display('Adjusting source covariance matrix.')
if any(strcmpi(InverseMethod, {'wmne','wMNE','dspm','dSPM','sloreta','sLORETA','mne','MNE'}))
    sspl = size(L,2);
    C_J = speye(sspl, sspl);
    C_J = spdiags(w, 0, C_J);
    trclcl = trace(L * C_J * L');
    C_J = C_J * (size(L,1)/ trclcl);
    Rc = chol(C_J,'lower');
    LW = L * Rc;
elseif  any(strcmpi(InverseMethod, {'loreta','LORETA'}))
    if nSource == 15028
        load 'LORETA15028.mat'
        trclcl = trace(L * C_J * L');
        C_J = C_J * (nSensor / trclcl);
        Rc = chol(C_J,'lower');
%     elseif nSource ==6003
%         load 'LORETA6003Rc.mat'
%     elseif nSource ==15002
%         load 'LORETA15002Rc.mat'
    else %||nSource~=15002
        if isempty(Cortex.VertConn)
            Cortex.VertConn = tess_vertconn(Cortex.Vertices, Cortex.Faces);
        end
        M = spdiags(ones(nSource,1),0,nSource,nSource) - sparse(double(Cortex.VertConn))./(sum(Cortex.VertConn,2)*ones(1,nSource));
        %==============Employ Depth Compensation for LORETA================%    
        sspl = size(L,2);
        C_J = speye(sspl, sspl);
        M = spdiags(w, 0, C_J)*M;
        C_J = M'*M+0.0001*trace(M'*M)/nSource*speye(nSource);
        %==============Not Employ Depth Compensation for LORETA=============%
%         C_J = M'*M+0.0001*trace(M'*M)/nSource*speye(nSource);

        C_J = inv(C_J);
        
        trclcl = trace(L * C_J * L');
        C_J = C_J * (nSensor / trclcl);
        Rc = chol(C_J,'lower');
    end
    LW = L * Rc;
end
clear C_J trclcl sspl itangential rnkC_noise
%% Estimate the Regularization Parameter via Data-Driven Process 
if Reg
    MAX_iter = 100;
    evidence = zeros(MAX_iter,1);
    gamma = 1;
    cost_old = 0;
    for iter = 1 : MAX_iter
        Cov_B = eye(nSensor) + gamma*LW*LW';
        gamma = norm(gamma*LW'/Cov_B*B,'fro')^2/(nSnap*trace(gamma*LW'/Cov_B*LW));
        cost = -(nSnap*log(det(Cov_B)) + trace(B*B'/Cov_B) + nSnap*nSensor*log(2*pi))/2;
        MSE = (cost - cost_old)/cost;
        cost_old = cost;
        evidence(iter) = cost;
        disp(['Iteration:   ',num2str(iter),'    MSE:  ',num2str(MSE),'    gamma: ',num2str(gamma),'    cost: ',num2str(cost)]);
        if abs(MSE) < 1e-5
            break;
        end
    end
    lambda2 = 1/gamma;
else
    lambda2 = SNR^(-2);
    Cov_B = eye(nSensor) + LW*LW'/lambda2;
    cost = -(nSnap*log(det(Cov_B)) + trace(B*B'/Cov_B) + nSnap*nSensor*log(2*pi))/2;
end
par.Cost = cost;
%% Compute SVD.
display('Computing SVD of whitened and weighted lead field matrix.')
[U,S,V] = svd(LW,'econ');
s = diag(S);
%lambda2 = 0.1*max(s);
% lambda2 = SNR^(-2);
ss = s ./ (s.^2 + lambda2);
Kernel = Rc * V * diag(ss) * U';
%% ===== WHITENED dSPM IMAGING KERNEL =====
% Compute dSPM operator.
if any(strcmpi(InverseMethod, {'dspm','dSPM'}))
display('Computing dSPM inverse operator.')
dspmdiag = sum(Kernel.^2, 2);
dspmdiag = sqrt(dspmdiag);
Kernel = bst_bsxfun(@rdivide, Kernel, dspmdiag);
%% ===== WHITENED sLORETA IMAGING KERNEL =====
% Compute sLORETA operator.
elseif any(strcmpi(InverseMethod, {'sloreta','sLORETA'}))
    display('Computing sLORETA inverse operator.')
    sloretadiag = sqrt(sum(Kernel .* L', 2));
    Kernel = bst_bsxfun(@rdivide, Kernel, sloretadiag);
end
disp(' ');
end