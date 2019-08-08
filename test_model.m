
dr = 0.01; %in kms
x1 = 0; x2 = 30;
z1 = 0; z2 = 15;
x = x1:dr:x2;
z = z1:dr:z2;
nx = numel(x);
nz = numel(z);
water_depth = 5.4; %%% added 3km extra padding (source depth  = 3 km)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xaxis = linspace(x1,x2,nx+2);
zaxis = linspace(z1,z2,nz+2);
[z2d, x2d] = meshgrid(zaxis, xaxis);
Vp = zeros(nz+2,1);
Vs = zeros(nz+2,1);
r = zeros(nz+2,1);

%Water velocity
dz = zaxis(2)-zaxis(1);
nzw = floor(water_depth/dz)+1;
Vp(1:nzw,1) = 1.5;
Vs(1:nzw,1) = 0;
r(1:nzw,1) = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make up the 1-D earth model
l1 =  0;
l2 = 0.5;
nz1 = nzw + floor(l1/dz)+1;
nz2 = nz1 + floor(l2/dz)+1;

Vp(nzw+1:nz1,1) = 2.2; %linspace(1.5,4,numel(nzw+1:nz1)); %
Vs(nzw+1:nz1,1) =  0.6  ; %linspace(1.5/2.2,4/2.2,numel(nzw+1:nz1));%
r(nzw+1:nz1,1) = 2;     %linspace(1,2.0,numel(nzw+1:nz1)); %
Vp(nz1+1:nz2,1) = 2.8; %linspace(1.5,4,numel(nzw+1:nz1)); %
Vs(nz1+1:nz2,1) =  1  ; %linspace(1.5/2.2,4/2.2,numel(nzw+1:nz1));%
r(nz1+1:nz2,1) = 2.2; %linspace(1,2.0,numel(nzw+1:nz1)); %
Vp(nz2+1:end,1) = 4.5; %linspace(1.5,4,numel(nzw+1:nz1)); %
Vs(nz2+1:end,1) = 2.6; %linspace(1.5/2.2,4/2.2,numel(nzw+1:nz1));%
r(nz2+1:end,1) = 2.5; %linspace(1,2.0,numel(nzw+1:nz1)); %

velp = repmat(Vp(:)',[nx+2 1]);
vels = repmat(Vs(:)',[nx+2 1]);
rhog = repmat(r(:)',[nx+2 1]);
[x2d, z2d] = meshgrid(xaxis,zaxis);


% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Velocity plot

% figure(21);clf
% plot(Vp,zaxis-3,'.-b');axis ij
% ylim([0 7]); xlim([0 7]); hold on
% figure(21);
% plot(Vs,zaxis-3,'.-r');axis ij
% ylim([0 7])
% grid on
% legend('Vp','Vs')
% xlabel('Velocity(km/s)')
% ylabel('Depth(km)')
% set(gca,'fontsize',13,'tickDir','out','XAxisLocation','bottom')
% 
