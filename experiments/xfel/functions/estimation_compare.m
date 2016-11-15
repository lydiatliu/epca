
function [ ] = estimation_compare(image_data,gamma, rank_guess, true_cov , num_trials)

%%this experiment shows the advantage of the est+shrink method in
%%covariance estimation, eigenvalue estimation and eigenvector(subspace) estimation
%%over vanilla methods.

%Inputs
%image_data - pxn data matrix with n columns containing p-dimensional clean
%               data vectors
%gamma - vector of values of gamma to try
%rank_guess - rank estimate used in PCA and denoising
%true_cov - true covariance matrix of underlying distribution
%num_trials - number of Monte Carlo trials for each data point

% set seed for reproducibility
seed = 1;
rng(seed)
p = size(image_data,1);

n = p./gamma;

%%% experiment results
%   1) covariance matrix estimation error (Spectral Norm of Difference)
spec_norm = zeros(4, num_trials, length(gamma));
%   2) covariance matrix estimation error (Frobenius Norm of Difference)
fro_norm = zeros(4, num_trials, length(gamma));
%   3) eigenvalue estimation percentage error - vector of length rank_guess
val_est = zeros(4,num_trials, length(gamma), rank_guess);
%   4) eigenvector absolute correlation - vector of length rank_guess
vect_est = zeros(3,num_trials, length(gamma), rank_guess);

%   5) projection matrix error (Spectral Norm of Difference)
pspec_norm = zeros(4, num_trials, length(gamma));
%   6) projection matrix estimation error (Frobenius Norm of Difference)
pfro_norm = zeros(4, num_trials, length(gamma));

for k = 1:length(gamma)
    for i = 1:num_trials
        clean_matrix = datasample(image_data', n(k))';
        [ est_shrunk_bc,est_shrunk, est, sample ] = one_poisson_expt( clean_matrix, rank_guess, true_cov);
        spec_norm(1,i, k) = est_shrunk_bc{1};
        spec_norm(2,i, k) = est_shrunk{1};
        spec_norm(3,i, k) = est{1};
        spec_norm(4,i, k) = sample{1};
        
        fro_norm(1,i, k) = est_shrunk_bc{2};
        fro_norm(2,i, k) = est_shrunk{2};
        fro_norm(3,i, k) = est{2};
        fro_norm(4,i, k) = sample{2};
       
        pspec_norm(2,i, k) = est_shrunk{5};
        pspec_norm(3,i, k) = est{5};
        pspec_norm(4,i, k) = sample{5};
        
        pfro_norm(2,i, k) = est_shrunk{6};
        pfro_norm(3,i, k) = est{6};
        pfro_norm(4,i, k) = sample{6};
        
        val_est(1,i, k,:) = est_shrunk_bc{3};
        val_est(2,i, k,:) = est_shrunk{3};
        val_est(3,i, k,:) = est{3};
        val_est(4,i, k,:) =sample{3};
        
        vect_est(1,i, k,:) = est_shrunk{4};
        vect_est(2,i, k,:) = est{4};
        vect_est(3,i, k,:) = sample{4};
       
    end
end

% save results
save(['results_estimation_',num2str(num_trials), '.mat'], 'spec_norm', 'fro_norm', 'pspec_norm', 'pfro_norm', 'val_est', 'vect_est');

%%% plots
%   1) covariance matrix estimation error (Spectral Norm of Difference)
%   against gamma, with error bars

figure;
temp = reshape(spec_norm(4,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'r--o', 'linewidth',3);

hold on;
temp = reshape(spec_norm(3,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'b:o', 'linewidth',3);


xlabel('log(n)');
ylabel('Difference from true covariance matrix');
hold on;
temp = reshape(spec_norm(2,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'g-.o', 'linewidth',3);

hold on;
temp = reshape(spec_norm(1,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'k--o', 'linewidth',3);


legend( {'Sample', 'Debiased', 'Recolored', 'Scaled'},'Location', 'best');
set(gca,'fontsize',15)
hold off;
title(['Error of covariance estimation (Spectral norm, ', num2str(num_trials), ' trials)']);
print(['images/bc_err_cov_est_spec_rank',num2str(rank_guess)],'-dpng');



%   2) covariance matrix estimation error (Frobenius Norm of Difference)
%   against gamma, with error bars
figure;
temp = reshape(fro_norm(4,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'r--o', 'linewidth',3);

hold on;
temp = reshape(fro_norm(3,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'b:o', 'linewidth',3);

xlabel('log(n)');
ylabel('Difference from true covariance matrix');
hold on;
temp = reshape(fro_norm(2,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'g-.o', 'linewidth',3);
hold on;
temp = reshape(fro_norm(1,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'k--o', 'linewidth',3);
legend( {'Sample', 'Debiased', 'Recolored', 'Scaled'},'Location', 'best');
set(gca,'fontsize',15)
hold off;
title(['Error of covariance estimation (Frobenius norm, ', num2str(num_trials), ' trials)']);
print(['images/err_cov_est_fro_rank',num2str(rank_guess)],'-dpng');


%   3) eigenvalue estimation percentage error - vector of length rank_guess
%   against gamma, with error bars: rank_guess subplots

counter = 1;
for j = 1:rank_guess
    subplot(2,floor(rank_guess/2),counter);
    temp = 100*reshape(val_est(4,:,:,j), [num_trials, length(gamma)]);
    errorbar(log10(n), mean(temp,1),std(temp),'r--o', 'linewidth',3);
    hold on;
    temp = 100*reshape(val_est(3,:,:,j), [num_trials, length(gamma)]);
    errorbar(log10(n), mean(temp,1),std(temp),'b:o', 'linewidth',3);
    xlabel('log(n)');
    ylabel('% estimation error');
    hold on;
    temp = 100*reshape(val_est(2,:,:,j), [num_trials, length(gamma)]);
    errorbar(log10(n), mean(temp,1),std(temp),'g-.o', 'linewidth',3);
    hold on;
    temp = 100*reshape(val_est(1,:,:,j), [num_trials, length(gamma)]);
    errorbar(log10(n), mean(temp,1),std(temp),'k--o', 'linewidth',3);
    if j==rank_guess
        legend( {'Sample', 'Debiased', 'Recolored', 'Scaled'},'Location', 'best');
    end
    set(gca,'fontsize',16)
    hold off;
    title(['Eigenvalue ', num2str(j),', (', num2str(num_trials), ' trials)'],'fontsize',16);
    counter=counter+1;
end
%fig.PaperPositionMode = 'auto';

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10*rank_guess/4 10];

print(['images/err_eval_est_rank',num2str(rank_guess)],'-dpng', '-r0');

%   4) eigenvector absolute correlation - vector of length rank_guess
%   against gamma, with error bars: rank_guess subplots
counter = 1;
for j = 1:rank_guess
    subplot(2,floor(rank_guess/2),counter);
    temp = reshape(vect_est(3,:,:,j), [num_trials, length(gamma)]);
    errorbar(log10(n), mean(temp,1),std(temp),'r--o', 'linewidth',3);
    %mean(temp,1);
    
    xlabel('log(n)');
    ylabel('Corr^2');
    hold on;
    temp = reshape(vect_est(2,:,:,j), [num_trials, length(gamma)]);
    errorbar(log10(n), mean(temp,1),std(temp),'b:o', 'linewidth',3);
    %mean(temp,1);
    
    hold on;
    temp = reshape(vect_est(1,:,:,j), [num_trials, length(gamma)]);
    errorbar(log10(n), mean(temp,1),std(temp),'g-.o', 'linewidth',3);
    %mean(temp,1);
    if j==rank_guess
        legend( {'Sample',  'Debiased', 'Recolored'},'Location', 'best');
    end
    set(gca,'fontsize',16)
    hold off;
    title(['Eigenvector ', num2str(j),' (', num2str(num_trials), ' trials)'],'fontsize',16);
    counter=counter+1;
end
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10*rank_guess/4 10];

print(['images/err_evec_est_rank',num2str(rank_guess)],'-dpng', '-r0');



%   5) projection matrix estimation error (Spectral Norm of Difference)
%   against gamma, with error bars
figure;
temp = reshape(pspec_norm(4,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'r--o', 'linewidth',3);
hold on;
temp = reshape(pspec_norm(3,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'b:o', 'linewidth',3);

xlabel('log(n)');
ylabel('Difference from true projection matrix');
hold on;
temp = reshape(pspec_norm(2,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'g-.o', 'linewidth',3);

legend( {'Sample', 'Debiased', 'Recolored'},'Location', 'best');
set(gca,'fontsize',15)
hold off;
title(['Error of subspace estimation (Spectral norm, ', num2str(num_trials), ' trials)']);
print(['images/err_proj_est_spec_rank',num2str(rank_guess)],'-dpng');



%   6) projection matrix estimation error (Frobenius Norm of Difference)
%   against gamma, with error bars
figure;

temp = reshape(pfro_norm(4,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'r--o', 'linewidth',3);
hold on;
temp = reshape(pfro_norm(3,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'b:o', 'linewidth',3);

xlabel('log(n)');
ylabel('Difference from true projection matrix');
hold on;
temp = reshape(pfro_norm(2,:, :), [num_trials, length(gamma)]);
errorbar(log10(n), mean(temp,1),std(temp),'g-.o', 'linewidth',3);
legend( {'Sample', 'Debiased', 'Recolored'},'Location', 'best');
set(gca,'fontsize',15)
hold off;
title(['Error of subspace estimation (Frobenius norm, ', num2str(num_trials), ' trials)']);
print(['images/err_proj_est_fro_rank',num2str(rank_guess)],'-dpng');


end