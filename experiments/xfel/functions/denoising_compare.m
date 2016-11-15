function [ ] = denoising_compare(image_data, gamma, rank_guess)

%%compares denoising methods (ePCA + EBLP vs Vanilla PCA)
% for different values of gamma

%Inputs
%image_data - pxn data matrix with n columns containing p-dimensional clean
%data vectors
%gamma - vector of values of gamma to try
%rank_guess - rank estimate used in PCA and denoising


mse1 = 1:length(gamma);
mse2 = 1:length(gamma);
mse3 = 1:length(gamma);
mse4 = 1:length(gamma);
mse0 = 1:length(gamma);
n = 1:length(gamma);

p = size(image_data,1);
edge = sqrt(p);


for u = 1:length(gamma)
    n(u) = p/gamma(u);
    clean_matrix = datasample(image_data', n(u))';
    %%noisy matrix
    noisy_matrix = poissrnd(clean_matrix,p,n(u));

    [~, ~, sample_mean, sample_cov] = estimate_stats2(clean_matrix,noisy_matrix);
    estim_cov = sample_cov-diag(sample_mean);
    [sample_eigvec,~] = eig(sample_cov);
    
    % Scaled covariance matrix
    [ white_covar,white_est_eigvec,~] = poisson_whiten_and_shrink( estim_cov, sample_mean, gamma(u), rank_guess);
    [ white_covar2,~,~] = poisson_whiten_and_shrink( estim_cov, sample_mean, gamma(u), rank_guess,true);

    %% Get PCA projections
    
    %colormap low high
    high = max(max(clean_matrix));
    low = 0;
    
    figure;
    subplot(1,5,1);
    projected_images = projections(n(u), noisy_matrix, sample_eigvec(:,(p-rank_guess+1):p));
    mse0(u) = sum(sum((clean_matrix-projected_images).^2))/n(u)/p;
    plotfaces(projected_images(:,1:9),3, 3, false, edge);
    caxis manual
    caxis([low high])
    title({'Vanilla projection', ['(Sample evec), ','$\gamma = $', num2str(gamma(u))]}, 'Interpreter', 'latex','fontsize',9 );

    
    subplot(1,5,2);
    projected_images = projections(n(u), noisy_matrix, white_est_eigvec);
    mse1(u) = sum(sum((clean_matrix-projected_images).^2))/n(u)/p;
    plotfaces(projected_images(:,1:9),3, 3, false, edge);
    caxis manual
    caxis([low high])
    title({'Vanilla projection','(Recolored evec)'}, 'Interpreter', 'latex','fontsize',9 );

    %% get EBLP projections
    projected_images = wiener_filter( sample_mean, sample_cov, noisy_matrix );
    mse3(u) = sum(sum((clean_matrix-projected_images).^2))/n(u)/p;
    subplot(1,5,3);
    plotfaces(projected_images(:,1:9),3, 3, false, edge);
    caxis manual
    caxis([low high])
    title({'EB linear predictor','(Sample evec)'}, 'Interpreter', 'latex','fontsize',9 );
    
    
    %%recolored cov
    projected_images = wiener_filter( sample_mean, white_covar, noisy_matrix );
    mse2(u) = sum(sum((clean_matrix-projected_images).^2))/n(u)/p;
    
    
    %%scaled cov
    projected_images = wiener_filter( sample_mean, white_covar2, noisy_matrix );
    mse4(u) = sum(sum((clean_matrix-projected_images).^2))/n(u)/p;
    subplot(1,5,4);
    plotfaces(projected_images(:,1:9),3, 3, false, edge);
    caxis manual
    caxis([low high])
    title({'EB linear predictor','(Scaled Covariance)'}, 'Interpreter', 'latex','fontsize',9 );
    
    subplot(1,5,5);
    plotfaces(clean_matrix,3, 3, false, edge);
    caxis manual
    caxis([low high])
    title('Clean images', 'Interpreter', 'latex','fontsize',9);
    
    colormap jet
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 9 4];
    
    
    print(['images/denoise_viz_', num2str(u), '_rank', num2str(rank_guess)],'-dpng', '-r0');

 
end
   
rng(2); a = {'--',':','-.','-'};
figure;
plot(log10(n), mse0, 'Color',rand(1,3),'LineStyle',a{1},'linewidth',4);
hold on;
plot(log10(n), mse1, 'Color',rand(1,3),'LineStyle',a{2},'linewidth',4);
hold on;
plot(log10(n), mse4,'Color',rand(1,3),'LineStyle',a{4},'linewidth',4);
xlabel('log(n)')
ylabel('MSE');
legend({'Vanilla projection (Sample evec)', 'Vanilla projection (Recolored evec)','EBLP (Scaled covariance)'} ,'Location', 'best')
title('MSE of denoising methods')
set(gca,'fontsize',16)
print(['images/denoising_compare_rank',num2str(rank_guess)],'-dpng', '-r0');


end