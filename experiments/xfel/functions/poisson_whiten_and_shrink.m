function [ covar,eigvec,eigval ] = poisson_whiten_and_shrink( debias_cov, sample_mean, gamma, rank_guess, bias_correction)
%% given a debiased covariance matrix, performs whitening and shrinkage
%% returns the resulting covariance matrix, eigenvector, eigenvalues (as a diagonal matrix), after de-whitening.
%% also implements heuristics to prevent eigenvalue reordering due to shrinkage

if ~exist('bias_correction','var')
    bias_correction = false;
end

addpath('../core')
p = length(sample_mean);

% whitening
D = diag(sample_mean.^(-1/2));
D(isinf(D)) = 0;
white_cov = D*debias_cov*D;
n_sig = floor(1.5*rank_guess);  % keep more than just the top rank_guess e-vectors
[V,E] = eigs(white_cov,n_sig);
E = diag(E); E = E(1:n_sig); V = V(:,1:n_sig);
[E,ind] = sort(E,'descend'); % sort the e-values in descending order after each call to eigs
V0 = V(:, ind);


%% number of eigenvalues to shrink and keep = rank_guess

%shrink eigenvalues and recompute eigenvectors
white_shr_eval = op_norm_shrink2( E+1, gamma)-1;
D = diag(sample_mean.^(1/2));
recolored_covar = (D*V0*diag(sqrt(white_shr_eval)))*(D*V0*diag(sqrt(white_shr_eval)))';  % M*M' to avoid loss of symmetry for numerical reasons
[V_col,E_col] = eigs(recolored_covar, n_sig);
E_col= diag(E_col);
[recolor_eval,ind] = sort(E_col,'descend');
V_col = V_col(:, ind); 

%% corrected estimator of the true spike
if bias_correction
    [~,c2] = standard_spiked_forward(white_shr_eval,gamma);
    s2 = 1-c2;
    tau = (sum(sample_mean)*white_shr_eval)./(p*recolor_eval);
    alpha = ((1-s2.*tau)./c2).*(c2>0)+ (c2<=0);
    white_shr_eval_corrected = alpha.*recolor_eval;
    eigval = diag(white_shr_eval_corrected(1:rank_guess));

else
    eigval = diag(recolor_eval(1:rank_guess));
end

eigvec = V_col(:,1:rank_guess); % keep just the top rank_guess e_vectors
temp = sqrt(eigval);
covar = (eigvec*temp)*(eigvec*temp)';

end

