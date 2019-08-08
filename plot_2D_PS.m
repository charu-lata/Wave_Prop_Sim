%%% P and S wave plot
figure(1); clf
load rwb_colormap.mat
pp  = hypot(Up,Wp); ss = hypot((U-Up),(W-Wp));
maxp = max(abs(pp(:)));
maxs = max(abs(ss(:)));

Up_plot = Up;
Wp_plot = Wp;
U_plot = U;
W_plot = W;
    
tempp1 = U_plot-Up_plot; tempp2 = W_plot-Wp_plot;
maxm1 = max(3,max([Up_plot(:);Wp_plot(:);tempp1(:);tempp2(:)]));

subplot(311)
fill([-13 13 13 -13],[water_depth water_depth 0 0]-3,[0.85 0.85 1],'edgecolor','none'); hold on
% if maxs>0.01*maxp
if maxp>0
    h1=quiver(xaxis-source.x,zaxis-3, (U_plot-Up_plot)'./maxm1, (W_plot-Wp_plot)'./maxm1,0,'color','r'); %exp(0.1*t).*
    % end
    hold on
    
    % if
    h2=quiver(xaxis-source.x,zaxis-3, (Up_plot)'./maxm1, (Wp_plot)'./maxm1,0,'color','b') ; %exp(0.1*t).*
    scale= t/4;
    hU1 = get(h1,'UData');
    hV1 = get(h1,'VData');
    set(h1,'UData',scale*hU1,'VData',scale*hV1)
    hU2 = get(h2,'UData');
    hV2 = get(h2,'VData');
    set(h2,'UData',scale*hU2,'VData',scale*hV2)
end

contour(xaxis-source.x,zaxis-3,l2m');
plot(trace.range,trace.z-3,'ok','Markersize',2.5)
plot(0,source.z-3,'ob','MarkerFaceColor','r','MarkerSize',3)
title(['Time = ',num2str(t),' s'])
xlabel('Range, km')
ylabel('Depth, km')
colormap(rwb)
axis ij
axis image
axis([0 13 0 6])
% axis([0 19 3 9])
% axis([-1.5 1.5 2 2.8])
% axis([x1-source.x x2-source.x z1 4])
caxis([-.5 .5])
set(gca,'TickDir','out','fontsize',14);
drawnow
hold off


%%%%%%%%%%%%%%%%% Plot x-t traces %%%%%%%%%%%%%%%%%%%%%%%%%%%

tvec = 0:dt:total_time;
v_red = 2;  range_x = [trace.range];
ax1 = subplot(3,1,2); xlabel('Range, km'); ylabel(strcat(['Time - Range/',num2str(v_red),', s'])); title('Vertical'); axis tight; set(ax1,'fontsize',14,'Ydir','Normal','tickDir','out');  box on; grid on; hold(ax1,'on')
ax2 = subplot(3,1,3); xlabel('Range, km'); ylabel(strcat(['Time - Range/',num2str(v_red),', s'])); title('Radial'); axis tight; set(ax2,'fontsize',14,'Ydir','Normal','tickDir','out');  box on; grid on; hold(ax2,'on')
seisx = Ut';
seisz = Wt';
[~,seisx1] = filter_butter([],seisx,10,35,24,1,1,dt);
[~,seisz1] = filter_butter([],seisz,10,35,24,1,1,dt);
if max(seisz(:))>0
    for i=1:numel(range_x); maxx = max([seisx(:,i);seisz(:,i)]/0.5);  %maxx = max(1,max([seisx1(:,i);seisz1(:,i)]));
        if range_x(i)>0
            seisx1 = -seisx1;
        end
        tt = tvec-abs(range_x(i))/v_red;
        fill(ax1,seisz1(:,i)./maxx+range_x(i),tt,'k');
        fill(ax2,seisx1(:,i)./maxx+range_x(i),tt,'k');
    end
    set(ax1,'ylim',[-3 4.2],'xlim',[0 13],'ytick',-3:1:4) ;
    set(ax2,'ylim',[-3 4.2],'xlim',[0 13],'ytick',-3:1:4) ;
end
return
