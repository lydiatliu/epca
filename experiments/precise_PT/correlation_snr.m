%% Plot Correlation
%cd('C:\Git\poissonpca\experiments\precise_PT')
addpath('../../software')
addpath('../random_subspace/functions');
%% MC simu: spike on grid
gamma = 1/2;
rng(2);
p = 500;
n = floor(p/gamma);
v = linspace(-1,1,p)'; v = v/norm(v);
u = linspace(1,3,p)';
PT  = sqrt(gamma)/(v'*diag(u)^(-1)*v);
c = linspace(0,10,20)';
ells = c.*(v'*v); %the true spikes we want to estimate

n_monte = 100;
correlation2_recolor = zeros(length(c),n_monte);
correlation2_debias = zeros(length(c),n_monte);
correlation2_sample = zeros(length(c),n_monte);
SNR_improvement =zeros(length(c),n_monte);
X_pos = zeros(length(c),n_monte);

x = tic;
print_iter=1;
for i=1:length(c)
    if print_iter==1
        str = sprintf('spike %d out of %d. time: %.2f sec \n',i,length(c), toc(x));
        fprintf(str);
    end
    for j=1:n_monte
        X = ones(n,1)*u'+c(i)^(1/2)*sqrt(3)*(2*rand(n,1)-1)*v';
        mini = min(min(X)); %need this >0
        if mini>0
            X_pos(i,j)=1;
        end
        %%% TRUE EIGENVECTOR
        true_cov = cov(X);
        [true_V,~] = eigs(true_cov,1); % top eigenvector
        
        Y = poissrnd(X);
        
        %%% SAMPLE EIGENVECTOR
        sample_cov = cov(Y);
        [sample_V,~] = eigs(sample_cov,1); % top eigenvector
        correlation2_sample(i,j) = (true_V'*sample_V)^2;
        
        %%% DEBIASED EIGENVECTOR
        debias_cov = cov(Y) - diag(mean(Y));
        [debias_V,~] = eigs(debias_cov,1); % top eigenvector
        correlation2_debias(i,j) = (true_V'*debias_V)^2;
        
        %%% RECOLORED EIGENVECTOR
        [~,white_shr_eval,white_vec,~,~,recolor_eval,recolor_vec,D_n]=...
            exp_fam_pca(Y, 'poisson');
        correlation2_recolor(i,j) = (true_V'*recolor_vec(:,1))^2;
    end
end
%% plot
% require that means are positive
pos = mean(X_pos,2);
max_pos =find(pos,1,'last');
%% spike estim
mean_correlation_s = mean(correlation2_sample, 2);
mean_correlation_d = mean(correlation2_debias, 2);
mean_correlation_r = mean(correlation2_recolor, 2);
rng(2); a = {'--',':','-.','-'};
figure, hold on
h1 = plot(ells(1:max_pos) ,mean_correlation_s(1:max_pos),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(ells(1:max_pos), mean_correlation_d(1:max_pos),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
h3 = plot(ells(1:max_pos) ,mean_correlation_r(1:max_pos),'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{4});
ylabel('Corr^2')
xlabel('Pop Spike')
m = max(1,max(mean_correlation_r(1:max_pos)));
ylim([0,m])
%plot PT location
SP=PT; %your point goes here
x=[SP,SP];
y=get(gca,'Ylim');
pt = plot(x,y,'linewidth',2,'color',rand(1,3));
set(gca,'fontsize',20)
xlim([min(ells(1:max_pos)),max(ells(1:max_pos))])
legend([h1,h2, h3,pt],{'Sample', 'Debiased', 'Recolored','Phase Transition'},'location','Best')
plot([0 1],[5 1], 'k-')
filename = sprintf( 'poisson_corr_snr_gamma=%.2f_PT=%.2f_n=%d_n_monte=%d',...
    gamma ,PT,n,n_monte);
%%
sf(filename)