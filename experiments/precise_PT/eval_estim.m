%% Eigenvalue estimation: experiments with whitening
%cd('C:\Git\poissonpca\experiments\precise_PT')
addpath('../../software')
%% MC simu: spike on grid
gamma = 1/2;
rng(2);
p = 500;
n = floor(p/gamma);
v = linspace(-1,1,p)'; v = v/norm(v);
%u = ones(p,1);
u = linspace(1,3,p)';
PT_white  = sqrt(gamma)/(v'*diag(u)^(-1)*v);
c = linspace(0,3,20)';
ells = c.*(v'*v); %the true spikes we want to estimate

n_monte = 100;
shr_spike_recolor =zeros(length(c),n_monte);
shr_spike =zeros(length(c),n_monte);
shr_spike_debias =zeros(length(c),n_monte);
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
        Y = poissrnd(X);
        [~,white_shr_eval,~,~,~,recolor_eval,~,D_n]=...
            exp_fam_pca(Y, 'poisson');
        shr_spike_recolor(i,j) =recolor_eval; %recolor estimation of spike
        [lambda,c2] = standard_spiked_forward(white_shr_eval,gamma);
        s2 = 1-c2;
        tau = (sum(D_n)*white_shr_eval)/(p*recolor_eval);
        if c2>0
            alpha = (1-s2*tau)/c2;
        else
            alpha=1;
        end
        shr_spike(i,j) =alpha*recolor_eval;
        Y_bar = mean(Y);
        Y_c = Y - ones(n,1)*Y_bar;
        S_d = 1/n*(Y_c'*Y_c)-diag(Y_bar);
        shr_spike_debias(i,j) = eigs(S_d,1);
        SNR_improvement(i,j) =tau/alpha;
    end
end
%% plot
% require that means are positive
pos = mean(X_pos,2);
max_pos =find(pos,1,'last');
%% spike estim
mean_shr_spike_recolor =mean(shr_spike_recolor,2);
mean_shr_spike_debias =mean(shr_spike_debias,2);
mean_shr_spike =mean(shr_spike,2);
rng(2); a = {'-','--','-.',':'};
figure, hold on
h1 = plot(ells(1:max_pos) ,ells(1:max_pos),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(ells(1:max_pos), mean_shr_spike_debias(1:max_pos),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
h3 = plot(ells(1:max_pos), mean_shr_spike_recolor(1:max_pos),'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{3});
h4 = plot(ells(1:max_pos) ,mean_shr_spike(1:max_pos),'linewidth',4,'color',rand(1,3));
set(h4,'LineStyle',a{4});
xlabel('Pop Spike')
ylabel('Est Spike')
m = max(ells(1:max_pos));
m = max(m,max(mean_shr_spike_recolor(1:max_pos)));
m = max(m, max(mean_shr_spike(1:max_pos)));
m = max(m, max(mean_shr_spike_debias(1:max_pos)));
ylim([0,m])
legend([h1,h2,h3,h4],{'True','Debiased','Recolored','Scaled'},'location','Best')
%plot PT location
SP=PT_white; %your point goes here
x=[SP,SP];
y=get(gca,'Ylim');
plot(x,y,'linewidth',2,'color',rand(1,3))
set(gca,'fontsize',20)
xlim([min(c(1:max_pos)),max(c(1:max_pos))])
filename = sprintf( 'poisson_spike_estim_gamma=%.2f_PT=%.2f_n=%d_n_monte=%d',...
    gamma ,PT_white,n,n_monte);
%%
sf(filename)
%% SNR improvement
mean_SNR_improvement =mean(SNR_improvement,2);
rng(2); a = {'-','--','-.',':'};
figure, hold on
h1 = plot(ells(1:max_pos) ,mean_SNR_improvement(1:max_pos),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
xlabel('Pop Spike')
ylabel('Est SNR Improment')
%plot PT location
SP=PT_white; %your point goes here
x=[SP,SP];
ylim([0,max(mean_SNR_improvement(1:max_pos))])
y=get(gca,'Ylim');
plot(x,y,'linewidth',2,'color',rand(1,3))
set(gca,'fontsize',20)
xlim([min(c(1:max_pos)),max(c(1:max_pos))])
ylim([min(mean_SNR_improvement(1:max_pos)),max(mean_SNR_improvement(1:max_pos))])
filename = sprintf( 'poisson_SNR_improvement_estim_gamma=%.2f_PT=%.2f_n=%d_n_monte=%d',...
    gamma ,PT_white,n,n_monte);
%%
sf(filename)