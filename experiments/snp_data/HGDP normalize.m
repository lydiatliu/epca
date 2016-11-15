%SNP data analysis: compare HWE normalization and standardization on HGDP
addpath('../../software')
addpath('./process_HGDP/')
%% load data
clear, clc
data_DIR='.\process_HGDP\';
load([data_DIR 'HGDP_small.mat'],'X','Y','pops','pop_code');
X = [X;Y];
clear Y
%% look at a subset of one population
sub =1;
if sub==1
    sample_inds = [1:10;31:40];
    X = X(sample_inds,:);
end
%%
[n,p] = size(X);
%%
x = tic;
cor =1.1;
r = 2;
remove_const = 1;

normalization = 'SD';
[X_sd,ind,~] = corr_prune(X,cor,normalization, remove_const); %ctr+standardize
% compute SVD of original matrix
[U,S,V] = svd(X_sd,'econ');
E = diag(S.^2); %all eigenvalues
PC_sd = X_sd*V(:,1:r);

normalization = 'HW';
[X_hw,ind,~] = corr_prune(X,cor,normalization, remove_const); %ctr+standardize
% compute SVD of original matrix
[U,S,V] = svd(X_hw,'econ');
E = diag(S.^2); %all eigenvalues
PC_hw = X_hw*V(:,1:r);

inds = [ones(floor(n/2),1);ones(floor(n/2),1)];
figure, hold on
h1= gscatter([PC_sd(:,1);PC_hw(:,1)], [PC_sd(:,2); PC_hw(:,2)],[inds;2*inds])
xlabel('PC 1')
ylabel('PC 2')
set(gca,'fontsize',14)
str = sprintf('HGDP_PC_scores_r=%d_norm=%s_n=%d',r,normalization,n)
sf(str)

figure, hold on
subplot(1,2,1)
h1= gscatter(PC_sd(:,1), PC_sd(:,2))
xlabel('PC 1')
ylabel('PC 2')
legend('SD')
set(gca,'fontsize',14)
subplot(1,2,2)
h2= gscatter(PC_hw(:,1), PC_hw(:,2))
xlabel('PC 1')
ylabel('PC 2')
legend('HW')
set(gca,'fontsize',14)
str = sprintf('HGDP_PC_scores_subplots_r=%d_norm=%s_n=%d',r,normalization,n)
sf(str)
toc(x);
