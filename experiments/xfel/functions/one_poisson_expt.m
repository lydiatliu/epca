function [ est_shrunk_bc, est_shrunk, est, sample ] = one_poisson_expt( clean_matrix, rank_guess, true_cov)
%one_poisson_expt: takes in clean_data, a pxn matrix, and
%   performs one Poisson PCA experiment. For each method, returns cell
%   containing
%   1) covariance matrix estimation error (Spectral Norm of Difference)
%   2) covariance matrix estimation error (Frobenius Norm of Difference)
%   3) eigenvalue estimation percentage error - vector of length rank_guess
%   Assumes that clean_data lies on a low rank subspace of rank, rank_guess
%%%%
%%%% Initialize the return values
%%%%
est_shrunk_bc = cell(1,3);
est_shrunk = cell(1,6);
est = cell(1,6);
sample = cell(1,6);


%%get dimension information from clean_matrix
p = size(clean_matrix,1);
n = size(clean_matrix,2);
gamma = p/n;

%%Add poisson noise
noisy_matrix = poissrnd(clean_matrix,p,n);

[~, ~, sample_mean, sample_cov] = estimate_stats2(clean_matrix,noisy_matrix);
estim_cov = sample_cov-diag(sample_mean);
[estim_eigvec,estim_eigval] = eig(estim_cov);
[sample_eigvec,sample_eigval] = eig(sample_cov);
[true_eigvec,true_eigval] = eig(true_cov);

[ white_covar,white_est_eigvec,white_est_eigval ] = poisson_whiten_and_shrink( estim_cov, sample_mean, gamma, rank_guess);
[ white_covar2,~,white_est_eigval2 ] = poisson_whiten_and_shrink( estim_cov, sample_mean, gamma, rank_guess, true);


%   1) covariance matrix estimation error (Spectral Norm of Difference)
est_shrunk{1} = norm(true_cov-white_covar);
est{1}  =norm(true_cov-estim_cov) ;
sample{1}  =norm(true_cov-sample_cov);
est_shrunk_bc{1}= norm(true_cov-white_covar2);

%   2) covariance matrix estimation error (Frobenius Norm of Difference)
est_shrunk_bc{2}= norm(true_cov-white_covar2, 'fro') ;
est_shrunk{2} = norm(true_cov-white_covar, 'fro') ;
est{2} =norm(true_cov-estim_cov, 'fro') ;
sample{2}  =norm(true_cov-sample_cov, 'fro') ;



%   3) eigenvalue estimation percentage error - vector of length rank_guess
temp = sort(diag(true_eigval), 'descend');
top_true_vals = temp(1:rank_guess);
temp0 = sort(diag(white_est_eigval2), 'descend');
est_shrunk_bc{3} = (temp0(1:rank_guess) - top_true_vals)./top_true_vals;
temp1 = sort(diag(white_est_eigval), 'descend');
est_shrunk{3} = (temp1(1:rank_guess) - top_true_vals)./top_true_vals;
temp2 = sort(diag(estim_eigval), 'descend');
est{3} = (temp2(1:rank_guess) - top_true_vals)./top_true_vals;
temp3 = sort(diag(sample_eigval), 'descend');
sample{3} = (temp3(1:rank_guess) - top_true_vals)./top_true_vals;

%   4) eigenvector sq correlation - vector of length rank_guess
truevec = fliplr(true_eigvec(:,p-rank_guess+1:p));
est_shrunk{4} = (diag(truevec'*white_est_eigvec)).^2;
debiased_eigenvectors = fliplr(estim_eigvec(:,(p-rank_guess+1):p));
est{4} = (diag(truevec'*debiased_eigenvectors)).^2;
sample_eigenvectors = fliplr(sample_eigvec(:,(p-rank_guess+1):p));
sample{4} = (diag(truevec'*sample_eigenvectors)).^2;


%   5) projection matrix estimation error (Spectral Norm of Difference)
true_proj = truevec*truevec';
est_shrunk{5} = norm(true_proj-white_est_eigvec*white_est_eigvec');
est{5}  =norm(true_proj-debiased_eigenvectors*debiased_eigenvectors') ;
sample{5}  =norm(true_proj-sample_eigenvectors*sample_eigenvectors');

%   6) projection matrix estimation error (Frobenius Norm of Difference)
est_shrunk{6} = norm(true_proj-white_est_eigvec*white_est_eigvec', 'fro') ;
est{6} =norm(true_proj-debiased_eigenvectors*debiased_eigenvectors', 'fro') ;
sample{6}  =norm(true_proj-sample_eigenvectors*sample_eigenvectors', 'fro') ;




end

