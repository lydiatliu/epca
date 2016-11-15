%% experiment 3: adds poisson noise and performs whitening and estimation
%%on n samples of xfel to estimate the eigenimages

function [ ] = experiment3_xfel(image_data, gamma, rank_guess, true_cov)

% set seed for reproducibility
seed = 1;
rng(seed)
%% extract pxn matrix of clean images. 
p = size(image_data,1);
edge = sqrt(p);

%%takes in pxN matrix image_data and returns pxn matrix, where each col is uniformly
%%sampled from the columns of image_data

for u = 1:length(gamma)
    n = p/gamma(u);
    clean_matrix = datasample(image_data', n)';
    
    %%noisy matrix
    noisy_matrix = poissrnd(clean_matrix,p,n);
    
    %% goals 1) estimate the eigenvalues (use MP to determine the rank of the
    %%low rank subspace) 2) get the estimated covariance matrix
    %%3) estimate the low rank subspace (the eigenvectors) 4) visualize the
    %%eigenvectors as eigenfaces.
    
    [~, ~, sample_mean, sample_cov] = estimate_stats2(clean_matrix,noisy_matrix);
    estim_cov = sample_cov-diag(sample_mean);
    [estim_eigvec,estim_eigval] = eig(estim_cov);
    [sample_eigvec,sample_eigval] = eig(sample_cov);
    [true_eigvec,true_eigval] = eig(true_cov);
    
    % with eigenvalue bias correction (related to improvement in SNR)
    [ ~,recolor_eigvec,recolor_eigval ] = poisson_whiten_and_shrink( estim_cov, sample_mean, gamma(u), rank_guess, true);
    
  
    %[ ~,recolor_eigvec,recolor_eigval ] = poisson_whiten_and_shrink(
    %estim_cov, sample_mean, gamma(u), rank_guess); %% no eigenvalue bias
    %correction
    
%     % whitening
%     D = diag(sample_mean.^(-1/2));
%     D(isinf(D)) = 0;
%     white_est_cov = D*estim_cov*D;
%     
%     % visualize whitening matrix:
%     %figure;
%     %histogram(diag(D));
%     %print(['images/whitening_matrix_viz_', num2str(person_label), '_', num2str(u)],'-dpng');
%     
%     %shrink eigenvalues and recompute eigenvectors
%     D = diag(sample_mean.^(1/2));
%     [V,E] = eigs(white_est_cov,rank_guess);
%     temp = zeros(rank_guess);
%     for j = 1:rank_guess
%         temp(j,j) = sqrt(ONshrink(E(j,j), gamma(u)));
%     end
%     white_covar = (D*V*temp)*(D*V*temp)';
%     [white_est_eigvec,white_est_eigval] = eigs(white_covar,rank_guess);

    %%eigenfaces
    figure;
    truevec = fliplr(true_eigvec(:,p-rank_guess+1:p));
    estim_eigvec = fliplr(estim_eigvec(:,(p-rank_guess+1):p));
    sample_eigvec = fliplr(sample_eigvec(:,(p-rank_guess+1):p));
    sign1 = diag(sign(diag(truevec'*recolor_eigvec)));
    sign2 = diag(sign(diag(truevec'*estim_eigvec)));
    sign3 = diag(sign(diag(truevec'*sample_eigvec)));
    
    h=ceil(rank_guess/3);
    w=rank_guess/h;
    
    % compute caxis high low
    %temp1 = max(max([truevec(:) estim_eigvec(:) sample_eigvec(:) recolor_eigvec(:)]));
    %temp2 = min(min([truevec(:) estim_eigvec(:) sample_eigvec(:) recolor_eigvec(:)]));
    high = max(max(truevec(:)), -min(truevec(:)));
    low = -high;
    
    subplot(1, 4, 1);
    vals = diag(true_eigval);
    vals = vals(p-rank_guess+1:p);
    plot_eigenvectors(truevec,w, h, false, edge, flipud(vals));
    title(['True Eigenimages'], 'Interpreter', 'latex', 'Fontsize', 14);
    caxis manual
    caxis([low high])
    axis off
    
    
    subplot(1, 4, 2);
    plot_eigenvectors(recolor_eigvec*sign1,w, h, false, edge,diag(recolor_eigval));
    title(['Recolored Eigenimages, ', '$\gamma = $', num2str(gamma(u))], 'Interpreter', 'latex', 'Fontsize', 14);
    caxis manual
    caxis([low high])
    axis off
    
    subplot(1, 4, 3);
    vals = diag(estim_eigval);
    vals = vals(p-rank_guess+1:p);
    plot_eigenvectors(estim_eigvec*sign2,w, h, false, edge, flipud(vals));
    title('Debiased Eigenimages', 'Interpreter', 'latex', 'Fontsize', 14);
    caxis manual
    caxis([low high])
    axis off
    
    subplot(1, 4, 4);
    vals = diag(sample_eigval);
    vals = vals(p-rank_guess+1:p);
    plot_eigenvectors(sample_eigvec*sign3,w, h, false, edge, flipud(vals));
    title('Sample Eigenimages', 'Interpreter', 'latex', 'Fontsize', 14);
    caxis manual
    caxis([low high])
    axis off
    
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 17/3*w 6];
    
    print(['images/eigenimages_',  num2str(u), '_rank', num2str(rank_guess)],'-dpng');
    
    abs(diag(truevec'*recolor_eigvec));
    abs(diag(truevec'*estim_eigvec));
    abs(diag(truevec'*sample_eigvec));
    
end

end
