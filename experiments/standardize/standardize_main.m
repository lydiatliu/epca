%% Compare standardization to whitening
addpath('../../software')
%% one MC trial of bulk
rng(2);
p = 500;
gamma = 1/2;
n = floor(p/gamma);
%modulate signal strength
M = 1;
m = M*ones(p,1);

%the choice of v should be important
v = M*linspace(-1,1,p)';

X = ones(n,1)*m'+(2*rand(n,1)-1)*v';
mini = min(min(X));
Y = poissrnd(X);
exp_fam_type = 'poisson';
[~,~,~,~,whitened_evals] = exp_fam_pca(Y,exp_fam_type);
null_e = whitened_evals(2:min(n,p));
null_mean = mean(null_e);
%% exp fam
figure, hold on
histogram(null_e, floor(2*sqrt(p)),'Normalization', 'pdf');
a = @(gamma) ((1-sqrt(gamma))^2);
b = @(gamma) ((1+sqrt(gamma))^2);
%MP
x = linspace(a(gamma), b(gamma), 100);
y = mp(x, gamma);
if(gamma>1)
    y = y*gamma;
end
plot(x,y,'r','LineWidth',2);
%%
sf('exp_fam_whitening')
%% usual
variances = var(Y)';
Y_n = n^(-1/2)*(Y-ones(n,1)*mean(Y))*diag(variances.^(-1/2));
norm_e = svd(Y_n).^2;
norm_e = norm_e(2:min(n,p));
norm_mean = mean(norm_e);
figure, hold on
histogram(norm_e, floor(2*sqrt(p)),'Normalization', 'pdf');
%histogram(norm_e/norm_mean, floor(2*sqrt(p)),'Normalization', 'pdf');
plot(x,y,'r','LineWidth',2);
%%
sf('usual_standardization')

%dividing by the mean is crucial
%because otherwise mean eigenvalue is 0.9, which is <1
%because we are "overcorrecting"
%after standardization, mean of all eigs is 1
%but there are spikes
%if instead we standardize using the exponential family method
%we can use the MP law directly

%% MC histogram of upper edge
rng(2);
p = 500;
gamma = 1/2;
n = floor(p/gamma);
M = 1;
m = M*ones(p,1);
v = M*linspace(-1,1,p)';
n_M = 100;
MP_edge_whiten = zeros(n_M,1);
MP_edge_sd = zeros(n_M,1);
exp_fam_type = 'poisson';
x = tic;
print_iter=1;
for i=1:n_M
    if print_iter==1
        str = sprintf('trial %d out of %d. time: %.2f sec \n',i,n_M, toc(x));
        fprintf(str);
    end
    X = ones(n,1)*m'+(2*rand(n,1)-1)*v';
    mini = min(min(X));
    Y = poissrnd(X);
    [~,~,~,~,whitened_evals] = exp_fam_pca(Y,exp_fam_type);
    MP_edge_whiten(i) = whitened_evals(2);
    
    variances = var(Y)';
    Y_n = n^(-1/2)*(Y-ones(n,1)*mean(Y))*diag(variances.^(-1/2));
    norm_e = svd(Y_n).^2;
    MP_edge_sd(i) = norm_e(2);
end

%% usual
figure, hold on
h1 = histogram(MP_edge_sd, floor(2*sqrt(p)));
h2 = histogram(MP_edge_whiten, floor(2*sqrt(p)));
%plot upper edge
SP=(1+sqrt(gamma))^2; %your point goes here
x=[SP,SP];
y=get(gca,'Ylim');
h3  = plot(x,y,'linewidth',2);
legend([h1,h2],{'Standardize','Whiten'})
set(gca,'FontSize',20)
%%
sf('upper_edge_hist')
