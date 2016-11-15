function[Y] = plot_eigenvectors(projected_images, numPerLine, ShowLine, createFig, faceL, eig_val)
% eig_val - vector of eigenvalues (length numPerLine*ShowLine)

if createFig 
    figure;
end

faceW = faceL;
faceH = faceL; 
Y = zeros(faceH*ShowLine,faceW*numPerLine); 
temp = projected_images';
for i=0:ShowLine-1 
  	for j=0:numPerLine-1 
    	Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(temp(i*numPerLine+j+1,:),[faceH,faceW]); 
  	end 
end 

imagesc(Y); colormap(gray);

for i=0:ShowLine-1 
  	for j=0:numPerLine-1 
        index = i*numPerLine+j+1;
        text((j+1)*faceW - 17,i*faceH+1+6, num2str(eig_val(index), '%.2f'),'BackgroundColor', 'w', 'fontsize', 10);
  	end 
end 

axis image

% print(['images/images_16x16_noisy_p', num2str(person_label)],'-dpng');
end