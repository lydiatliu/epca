%% runs XFEL experiments with dataset
addpath('../../core');
addpath('./functions');

%%load the true covariance matrix - assumes 'true_cov.mat' file exists
load('true_cov.mat');

% get p x n data matrix - assumes 'data.cxi' file exists
data = cxi2datamatrix( 'data.cxi' );

% choose average pixel intensity to scale to
average_intensity = 0.04;
data = average_intensity*data/mean(data(:)); %rescale data

gamma = [4 2 1 1/2 1/4 1/16];
rank_guess = 10;
num_trials = 50;

%estimation
estimation_compare(data, gamma, rank_guess, true_cov , num_trials);

%denoising
denoising_compare(data, gamma, rank_guess);

%eigenvectors plot
experiment3_xfel(data, gamma, rank_guess, true_cov);



