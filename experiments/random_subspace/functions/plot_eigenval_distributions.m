function [ X, Y ] = plot_eigenval_distributions( r , save_data)
% plots the eigenvalue distributions (Figure 1 and 2 in paper)
% set r=1 (constant vector case) or r=2 (spike model)

gamma = [1/4 1 4 ];
p= [100 1000 ];
n = floor((1./gamma)'*p);  %n is length(gamma) by length(p)

X = cell(length(p), length(gamma)); %all the clean images
Y = cell(length(p), length(gamma)); %all the noisy images
true_mean = cell(length(p), length(gamma)); %all the means
true_cov = cell(length(p), length(gamma)); %all the sample covariances
sample_cov = cell(length(p), length(gamma)); %all the sample covariances
sample_mean = cell(length(p), length(gamma)); %all the sample means
estim_cov = cell(length(p), length(gamma)); %all the estimated covariance

%%whitened stuff
white_sam_cov = cell(length(p), length(gamma));
white_est_cov = cell(length(p), length(gamma));
white_true_cov = cell(length(p), length(gamma));

white_est_eigval = cell(length(p), length(gamma));
white_est_eigvec = cell(length(p), length(gamma));
white_sam_eigval = cell(length(p), length(gamma));
white_sam_eigvec = cell(length(p), length(gamma));
white_true_eigval = cell(length(p), length(gamma));

shrunk_eigval = cell(length(p), length(gamma));
estim_eigval = cell(length(p), length(gamma));
estim_eigvec = cell(length(p), length(gamma));
sample_eigval = cell(length(p), length(gamma));
sample_eigvec = cell(length(p), length(gamma));
true_eigval = cell(length(p), length(gamma));
true_eigvec = cell(length(p), length(gamma));

for i = 1:length(p)
    for j = 1:length(gamma)
        if save_data
            [X{i,j},Y{i,j}] = simulate_poisson(n(j,i), p(i),r);
            [true_mean{i,j}, true_cov{i,j}, sample_mean{i,j}, sample_cov{i,j}] = estimate_stats2(X{i,j},Y{i,j});
        else
            [X, Y] = simulate_poisson(n(j,i), p(i),r);
            [true_mean{i,j}, true_cov{i,j}, sample_mean{i,j}, sample_cov{i,j}] = estimate_stats2(X,Y);
        end
        estim_cov{i,j} = sample_cov{i,j}-diag(sample_mean{i,j});
        [estim_eigvec{i,j},estim_eigval{i,j}] = eig(estim_cov{i,j});
        [sample_eigvec{i,j},sample_eigval{i,j}] = eig(sample_cov{i,j});
        [true_eigvec{i,j},true_eigval{i,j}] = eig(true_cov{i,j} );
    end
end

for i = 1:length(p)
    for j = 1:length(gamma)
        %%'whitening'
        D = diag(sample_mean{i,j}.^(-1/2));
        D(isinf(D)) = 0;
        white_sam_cov{i,j} = D*sample_cov{i,j}*D;
        white_est_cov{i,j} = D*estim_cov{i,j}*D;
        white_true_cov{i,j}= D*true_cov{i,j}*D;
        
        white_est_eigval{i,j} = eig(white_est_cov{i,j});
        white_sam_eigval{i,j} = eig(white_sam_cov{i,j});
        white_true_eigval{i,j}=eig(white_true_cov{i,j});
    end
end


%% utility
a = @(gamma) ((1-sqrt(gamma))^2);
b = @(gamma) ((1+sqrt(gamma))^2);

%%whitened eigenvalues plots

%whitened estimated eigenvalues.
figure;
counter=1;
for j = 1:length(gamma)
    for i = 1:length(p)
        subplot(length(gamma), length(p),counter);
        %histogram with pdf normalization - bar areas sum to 1
        eigenvalues = real(white_est_eigval{i,j});
        eigenvalues(eigenvalues < -1+1e-6) = []; % remove point mass at -1
        nbins = floor(2*(p(i))^(1/2));
        histogram(eigenvalues,nbins, 'Normalization', 'pdf');
        set(gca,'fontsize',12)
        
        hold on
        
        if r>= 2
            %%plot the true eigenvalue of the (transformed) covariance matrix
            temp = white_true_eigval{i,j};
            scatter([max(real(temp))] ,zeros(r-1,1),'MarkerEdgeColor','r', 'LineWidth',0.3);
            
            hold on;
            %plot max est eigenval
            scatter([max(eigenvalues)] ,zeros(r-1,1),30,'MarkerEdgeColor','none','MarkerFaceColor','r');
            hold on;
             if counter == length(gamma)*length(p)
                 legend({'Debiased EV dist.', 'True EV', 'Top Debiased EV', 'Top Debiased Shrunken EV'}, 'Location','Best');
             end
        end
        
        %MP
        x = linspace(a(gamma(j)), b(gamma(j)), 100);
        y = mp(x, gamma(j));
        if(gamma(j)>1)
            y = y*gamma(j);
        end
        plot(x-1,y,'r','LineWidth',2);
        hold off;
        
        title(['$\gamma$ = ',num2str(gamma(j)), ', $p$ = ', num2str(p(i))], 'Interpreter', 'latex');
        
        
        counter = counter+1;
    end
end

print(['images/whitened_estim_eigval_rank', num2str(r-1)],'-dpng');

end

