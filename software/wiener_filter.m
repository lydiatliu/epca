function [ X ] = wiener_filter( mean_hat, sigma_hat, Y,exp_fam_type, sigma_sq,bino_n )
% given the noisy data, mean estimate and covariance estimate,
% computes the estimated clean data, using the best linear estimator.
% This is described in Sec. 5. "Denoising" in the ePCA paper. Wiener filter
% = BLP

%Inputs:
% Y must be p x n (each sample is a column vector).
% mean_hat is p x 1. sigma_hat is p x p.
%exp_fam_type - distribution of the entries of Y, a type of exponential family
%         possibilities: 'normal', 'poisson','binomial'
%sigma_sq,bino_n - parameters for respective exponential families

if ~exist('exp_fam_type','var')
    exp_fam_type = 'poisson';
end
%mean-variance map
switch exp_fam_type
    case 'normal'
      if ~exist('sigma_sq','var')
        sigma_sq = 1;
    end
        V = @(x) sigma_sq;
    case 'poisson'
        V = @(x) x;
    case 'binomial'
        if ~exist('bino_n','var')
            bino_n = 1;
        end
        V = @(x) x.*(1-x./bino_n);
end

Y_bar = diag(mean_hat);
D_n = V(Y_bar);
n = size(Y,2);
%diagnostic: how many pixels have zero intensity?
%n_zero   = sum(mean_hat==0);

%at the moment, the regularization is 
%set to ridge by default
if ~exist('reg','var')
    reg = 'ridge';
end

switch reg
    case 'ridge'
        %approach 1: ridge regularization
        eps = 0.1;
        m=mean(diag(D_n));
        Sigma_reg = sigma_hat+(1-eps)*D_n+eps*m*eye(length(mean_hat));
        S_r_inv = Sigma_reg^(-1);
        X = sigma_hat*S_r_inv*Y + D_n*S_r_inv*repmat(mean_hat,1,n);
                
    case 'none'
        X = sigma_hat/(D_n+sigma_hat)*Y + D_n/(D_n+sigma_hat)*repmat(mean_hat,1,n);
        
end
end

