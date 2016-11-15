function[Y] = plotfaces(projected_images, numPerLine, ShowLine, createFig, faceL)
if nargin == 1
  numPerLine = 16; 
  ShowLine = 4;
  createFig = true;
  faceL = 32;
end

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
axis image

set(gca,'XTickLabel','','YTickLabel','')
ax=gca;
ax.GridColor = 'w';
ax.GridAlpha = 1;
ax.XTick = linspace(0, numPerLine*faceL, numPerLine+1);
ax.YTick = linspace(0, ShowLine*faceL, ShowLine+1);
grid on


% print(['images/images_16x16_noisy_p', num2str(person_label)],'-dpng');
end