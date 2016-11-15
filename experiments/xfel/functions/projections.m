function[projected_images] = projections(n, noisy, eigenvectors)

U = eigenvectors;
Y = noisy(:,1:n);
rowmean = mean(Y,2);
Y = Y-repmat(rowmean, 1,n);
projected_images = U*U'*Y;
projected_images = projected_images + repmat(rowmean, 1,n);

end
