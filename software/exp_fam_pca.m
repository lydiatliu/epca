function [S_recolored,white_shr_eval,white_V,white_U,white_eval,recolor_eval,recolor_V,D_n, ...
    white_shr_eval_scaled,estim_SNR_improvement] ...
    = exp_fam_pca(Y,exp_fam_type,r,sigma_sq,bino_n)

%ePCA: High Dimensional Exponential Family PCA
%computes a Principal Component decomposition of exponential-family valued data
%see the ePCA paper by Liu, Dobriban and Singer for details of the method

%Inputs
%Y - nxp data matrix with n rows containing p-dimensional noisy signal
%          vectors
%exp_fam_type - distribution of the entries of Y, a type of exponential family
%         possibilities: 'normal', 'poisson','binomial'
%r - rank estimate used in PCA and denoising
%sigma_sq,bino_n - parameters for respective exponential families

%Outputs
%S_recolored - whitened, shrunk and recolored covariance matrix
%white_shr_eval - whitened, shrunk eigenvalue (no recoloring)
%white_V - top r eigenvectors of whitened covariance matrix (also right
%                  singular vectors of centered and whitened data matrix)
%white_U - top r left singular vectors of centered and whitened data matrix
%white_eval - whitened eigenvalues (no shrinkage, or recoloring)
%                  (standard Marchenko-Pastur)
%recolor_eval,recolor_V - eigendecomposition of S_w_op
%D_n        - diagonal debiasing matrix
%white_shr_eval_scaled - scaling to white_shr_eval, to remove bias
%                   for estimating true spike
%estim_SNR_improvement - estimated SNR improvement due to whitening

[n,p] = size(Y);
gamma = p/n;

%impute missing data to column means
impute_missing=1;
if impute_missing
    m = nanmean(Y);
    for i=1:p
        ind_nan= isnan(Y(:,i));
        Y(ind_nan,i) = m(i);
    end
    vars = var(Y);
    Y = Y(:,vars>0); %remove columns with zero variance
end

if ~exist('r','var')
    r = 1;
end
if ~exist('exp_fam_type','var')
    exp_fam_type = 'normal';
    if ~exist('sigma_sq','var')
        sigma_sq = 1;
    end
end
%mean-variance map
switch exp_fam_type
    case 'normal'
        V = @(x) sigma_sq;
    case 'poisson'
        V = @(x) x;
    case 'binomial'
        if ~exist('bino_n','var')
            bino_n = 1;
        end
        V = @(x) x.*(1-x./bino_n);
end

Y_bar = mean(Y);
D_n = V(Y_bar);

%center data
Y_c = Y - ones(n,1)*Y_bar;
%standardize data
Y_w = n^(-1/2)*Y_c*diag(D_n.^(-1/2));
%note: Y_w'*Y_w = S_w + I_d

[U,sval,V] = svd(Y_w,'econ');
sval = diag(sval);
white_eval = sval.^2;
E = sval(1:r).^2;
%white_shr_eval = op_norm_shrink( E, gamma)-1;
white_shr_eval = op_norm_shrink2( E, gamma)-1;

white_V = V(:,1:r);
white_covar = white_V*diag(white_shr_eval)*white_V'; %rank r
S_recolored = diag(D_n.^(1/2))*white_covar*(diag(D_n.^(1/2)));
white_U = U(:,1:r);
D_n = D_n';

[recolor_V,recolor_eval] = eigs(S_recolored,r);
recolor_eval = diag(recolor_eval);
[recolor_eval,ind] = sort(recolor_eval,'descend');
recolor_V = recolor_V(:, ind);

%scaled estimator of the true spike
[~,c2] = standard_spiked_forward(white_shr_eval,gamma);
s2 = 1-c2;
tau = (sum(D_n)*white_shr_eval)./(p*recolor_eval);
alpha = zeros(r,1);
for i=1:r
    if c2(i)>0
        alpha(i) = (1-s2(i)*tau(i))/c2(i);
    else
        alpha(i)=1;
    end
end
white_shr_eval_scaled =alpha.*recolor_eval;
estim_SNR_improvement =tau./alpha;
