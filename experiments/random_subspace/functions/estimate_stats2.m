function[true_mean, true_cov, noisy_mean, noisy_cov] = estimate_stats2(clean, noisy)
%assumes clean is p x n

X = clean;
R = noisy;
n= size(X,2);

% clean covariance matrix;
true_cov = cov(X'); %pxp matrix
true_mean = sum(X,2)/n;

%%poisson noise
noisy_mean = sum(R,2)/n;
noisy_cov = cov(R');


end
