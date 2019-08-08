%%% Vertical& Horizontal velocity plots as quivers

load rwb_colormap.mat
figure(1);clf
% subplot(211)
pimag = squeeze(U(:,:)');
pimag = exp(0.2*t).*(pimag./max(abs(pimag(:))));
quiver(xaxis-source.x,zaxis,U',W',1.4)
title(['Horizontal displacment |  t = ',num2str(t),' s'])
xlabel('Distance, km')
ylabel('Depth, km')
axis image; 
axis ij
axis([x1-source.x x2-source.x z1 z2])
colormap(rwb)
caxis([-.5 .5])
drawnow
