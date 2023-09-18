function [W,H,Hc,dc_mvnmf_obj] = DC_MVNMF_L1_Norm(X,options)
% modified from "DC_MVNMF", L1-normalization is imposed on columns of U, 
% not L2-normalization like "DC_MVNMF".
% sum_v |Xv-WvHv|_F^2 + alpha*sum_v |Hv-[Y,H*]|_F^2
% + beta*sum_v tr(Hv*Lv*Hv')
% s.t. Wv,Hv >=0 
% refer to 2020-KBS
n_view = length(X);
for i = 1:n_view
    X{i} = X{i}';
    [mFea{i},nSmp] = size(X{i});
end

nSubSpace = options.nSubSpace;
dim = options.nClas*nSubSpace;
alpha = options.alpha;
beta = options.beta;

gnd_known = options.gnd_known;


nSmp_label = length(gnd_known);
nSmp_unlabel = nSmp - nSmp_label;
%% construct I_block
nClas = options.nClas;
if nSmp_label > 0
    I_label = zeros(nClas*nSubSpace,nSmp_label);
    for i = 1:nClas
        idx = i == options.gnd_known;
        I_label((i-1)*nSubSpace+1:i*nSubSpace,idx) = 1;
    end
    Y = I_label;
elseif nSmp_label == 0 && nSmp_unlabel > 0 
    Y = [];
end
Y = sparse(Y);
%% load graph laplacian of all samples
if beta > 0
    for i = 1:n_view
        Slu{i} = sparse(beta*options.Slu{i});
%         Du{i} = sum(Su{i},2);
        Dlu{i} = spdiags(full(sum(Slu{i},2)),0,nSmp,nSmp);
%         Du{i} = diag(Du{i});
        Llu{i} = Dlu{i} - Slu{i};
    end
else
    for i = 1:n_view
        Llu{i} = zeros(nSmp,nSmp);
    end
end
%% initialize W{i} H{i} Hc
init_nmf = 0;
if ~init_nmf
    W = cell(1,n_view);
    H = cell(1,n_view);
    for i = 1:n_view
        W{i} = rand(mFea{i},dim);
        H{i} = rand(dim,nSmp);
    end
    Hc = rand(dim,nSmp_unlabel);
else
    opts_nmf = [];
    opts_nmf.maxIter = 100;
    opts_nmf.error = 1e-6;
    opts_nmf.nRepeat = 30;
    opts_nmf.minIter = 50;
    opts_nmf.meanFitRatio = 0.1;
    for i = 1:n_view
        [W{i},H{i}] = NMF(X{i},dim, opts_nmf, [], []);H{i} = H{i}';
        Hc = 1/n_view*H{i};
    end
%     Hc = rand(dim,nSmp_unlabel);
end
%% normalize W{i} H{i}
for v = 1:n_view
    Norm = 1;
    NormV = 0;
    [W{v},H{v}] = NormalizeUV(W{v}, H{v}, NormV, Norm);
end
%% ------------------- normalization ------------------------%
% normalize the column vectors of W and consequently convey the
% norm to the coefficient matrix H
% for v = 1:n_view
%     normW = max(1e-15,sqrt(sum(W{v}.*W{v},1)));normW = normW';
%     W{v} = W{v}*spdiags(normW.^-1,0,dim*nSubSpace,dim*nSubSpace);
%     H{v} = spdiags(normW,0,dim*nSubSpace,dim*nSubSpace)*H{v};
% end
%% calculate obj at step 0
dc_mvnmf_obj = [];
% obj = CalculateObj(X,W,H,Hc,Y,alpha,beta,p,Llu);
% dc_mvnmf_obj = [dc_mvnmf_obj,obj];
%% start optimization
max_iter = options.max_iter;
iter = 0;
while  iter<=max_iter 
    iter = iter + 1;
    for v = 1:n_view
    %################# updata U1 and V1 ###################%
    %--------------------- update U1 ----------------------%
        XDH = X{v}*H{v}';
        WHDH = W{v}*H{v}*H{v}';
        if alpha > 0 % consistent
            HcBar = [Y,Hc];
            HH = diag((HcBar*H{v}').*eye(dim));% Y3 in paper.
            XDH = XDH + alpha*repmat(HH',mFea{v},1);
            
            normW = max(1e-15,sum(W{v},1));normW = normW';
            Qv = spdiags(normW,0,dim,dim);
            Y1 = diag(Qv*((H{v}*H{v}').*eye(dim)));
            WHDH = WHDH + alpha*repmat(Y1',mFea{v},1);
        end
        if beta > 0 % laplacian
            normW = max(1e-15,sum(W{v},1));normW = normW';
            Qv = spdiags(normW,0,dim,dim);
            Y2n = diag(Qv*((H{v}*Slu{v}*H{v}').*eye(dim)));
            XDH = XDH + repmat(Y2n',mFea{v},1);
            Y2p = diag(Qv*((H{v}*Dlu{v}*H{v}').*eye(dim)));
            WHDH = WHDH + repmat(Y2p',mFea{v},1);
        end

        W{v} = W{v}.*(XDH./max(WHDH,1e-10)); % 3mk
    %------------------- normalization ------------------------%
    % normalize the column vectors of W and consequently convey the
    % norm to the coefficient matrix H
        normW = max(1e-15,sum(W{v},1));normW = normW';
        W{v} = W{v}*spdiags(normW.^-1,0,dim,dim);
        H{v} = spdiags(normW,0,dim,dim)*H{v};
    %--------------------- update V1 ----------------------%
        WXD = W{v}'*X{v}; % mnk or pk (p<<mn)
        WWHD = W{v}'*W{v}*H{v}; % mk^2
        if alpha > 0
            HcBar = [Y,Hc];
            WXD = WXD + alpha*HcBar;
            WWHD = WWHD + alpha*H{v};
        end
        if beta > 0 
            WXD = WXD + H{v}*Slu{v};
            WWHD = WWHD + H{v}*Dlu{v};
        end

        H{v} = H{v}.*(WXD./max(WWHD,1e-10));
    end
    %--------------------- update V* ----------------------%
    if alpha > 0
        sum_H = 0;
        for i = 1:n_view
            sum_H = sum_H + H{i}(:,nSmp_label+1:end);
        end
        Hc = sum_H/n_view;
    end

%         newobj = CalculateObj(X,W,H,Hc,Y,alpha,beta,p,Llu);
% %     %     differror = abs(newobj - objhistory(end))/abs(objhistory(end));
%         dc_mvnmf_obj = [dc_mvnmf_obj newobj]; %#ok<AGROW>
% %         disp(num2str(iter))
end

% for v = 1:n_view
%     Norm = 1;
%     NormV = 0;
%     [W{v},H{v}] = NormalizeUV(W{v}, H{v}, NormV, Norm);
% end

end
%==========================================================================

function obj = CalculateObj(X,W,H,Hc,Y,alpha,beta,p,Llu)
    n_view = length(X);
    
    obj_NMF = 0;
    for v = 1:n_view
        dX = W{v}*H{v}-X{v};
        obj_NMF = obj_NMF + sum(sqrt(sum(dX.^2,1)).^p);
    end
    
    obj_consis = 0;
    if alpha > 0
        for v = 1:n_view
            obj_consis = obj_consis + sum(sum((H{v}-[Y,Hc]).^2));
        end
    end
    obj_consis = alpha*obj_consis;

    obj_lap = 0;
    if beta > 0
        for v = 1:n_view
            obj_lap = obj_lap + sum(sum(H{v}.*(H{v}*Llu{v})));
        end
    end

    
    obj = obj_NMF + obj_consis + obj_lap;
end

function [U, V] = NormalizeUV(U, V, NormV, Norm)
    K = size(U,2);
    if Norm == 2
        if NormV
            norms = max(1e-15,sqrt(sum(V.^2,1)));
            V = spdiags(norms.^-1,0,K,K)*V;
            U = U*spdiags(norms,0,K,K);
        else
            norms = max(1e-15,sqrt(sum(U.^2,1)))';
            U = U*spdiags(norms.^-1,0,K,K);
            V = spdiags(norms,0,K,K)*V;
        end
    else
        if NormV
            norms = max(1e-15,sum(abs(V),1));
            V = spdiags(norms.^-1,0,K,K)*V;
            U = U*spdiags(norms,0,K,K);
        else
            norms = max(1e-15,sum(abs(U),1))';
            U = U*spdiags(norms.^-1,0,K,K);
            V = spdiags(norms,0,K,K)*V;
        end
    end
end
