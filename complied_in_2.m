clearvars
% close all

%%% input model x1 x2 z1 z2 nx nz lam mu b
test_model;

dt = 0.001; %in seconds grad:0.0025
total_time = 7;
S = dt/dr;
c1 = 9/8;
c2 = -1/24;

%%%stability check
stab = 1/sqrt(2);
for ix = 1:nx
    for iz = 1:nz
        if velp(ix,iz)*dt/dr > stab
            fprintf('reduce dt or increase dx \n')
            keyboard    
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SOURCE STRUCTURE
% Set the time sample step size
source.dt = 0.005; 
% Choose a time vector length
source.tlen = 0.1;
% Set the source position and radiusl
source.x = x2/2;
source.z = 3;
source.radius = 2*dr;


source.sample_rate = 1/source.dt;
fnyq = source.sample_rate/2;
source.freq = 10;
source.t = 0:source.dt:0.1; %0:source.dt:source.tlen;
source.nsamp = length(source.t);
pf  = pi^2*source.freq^2;
d = (1-2*pf*source.t.^2).*exp(-pf*source.t.^2);
% d = sin(2*pi*source.freq.*source.t);
d=d./max(d);
% figure(2);clf
% plot(source.t,d,'.-k')
source.sig = d;



%%interpolate source to new time sampling
% Time length of source:
tso = (source.nsamp-1)*source.dt;
tveco = linspace(0,tso,source.nsamp);
ns  = round(tso/dt)+1;		% length of resampled source time series
ts = (ns-1)*dt;            %Time length of reamapled source signature
tvec_s = linspace(0,ts,ns);
source.src = interp1(tveco,source.sig,tvec_s,'pchip');
source.nsamp = length(tvec_s);
source.t = tvec_s;
source.sample_rate = 1/source.dt;
source.ix = round((source.x-x1)/dr);
source.iz = round((source.z-z1)/dr);
source.irad = round(source.radius/dr);
amp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TRACE HEADER STRUCTURE
trace.num_traces = x2/0.25;
trace.x = linspace(x1,x2,trace.num_traces);
trace.z = ones(trace.num_traces,1).*water_depth+2*dr;
trace.dt = dt;
trace.tlen = total_time;
trace.nt = round(trace.tlen/dt) + 1;
% Compute nearest node of receiver in model space
trace.ix = round((trace.x-x1)/dr)+1;
trace.iz = round((trace.z-z1)/dr)+1;
trace.range =  trace.x-source.x;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

get_staggered_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialise FD and timestep
U = zeros(nx+2,nz+2);
W = zeros(nx+2,nz+2);
Txx = zeros(nx+2,nz+2);
Tzz = zeros(nx+2,nz+2);
Txz = zeros(nx+2,nz+2);
Sp = zeros(nx+2,nz+2);
Up = zeros(nx+2,nz+2);
Wp = zeros(nx+2,nz+2);
k = 0; t = 0; nfill = 0;
Ut = zeros(trace.num_traces,trace.nt); PPt = zeros(trace.num_traces,trace.nt);
Wt = zeros(trace.num_traces,trace.nt); SSt = zeros(trace.num_traces,trace.nt);

%%% video file
snapshot_times_filename = '../OUT/snaptimes_step_ss3_wide.avi';
vidfile = VideoWriter(snapshot_times_filename);
open(vidfile);


while k < total_time/dt
    if (k>=1 && k <= source.nsamp)
        for ix = max(2,source.ix-source.irad):min(nx+2,source.ix+source.irad)
            for iz = max(2,source.iz-source.irad):min(nz+2,source.iz+source.irad)
                x_temp = x1 + (ix)*dr;
                dxo = source.x-x_temp;
                z_temp = z1 + (iz)*dr;
                dzo = source.z-z_temp;
                rso = sqrt(dxo*dxo + dzo*dzo);
                arg = rso/source.radius;
                if (arg <= 1.0)
                    fac = cos(0.5*pi*arg);
                    fac = fac*fac;                    
                    Txx(ix+ioPx,iz+ioPz) = Txx(ix+ioPx,iz+ioPz) + amp*fac*source.src(k);
                    Tzz(ix+ioPx,iz+ioPz) = Tzz(ix+ioPx,iz+ioPz) + amp*fac*source.src(k);
                    Sp(ix+ioPx,iz+ioPz) = 1./(1 + (lam(ix+ioPx,iz+ioPz)/l2m(ix+ioPx,iz+ioPz))) * (Txx(ix+ioPx,iz+ioPz) + Tzz(ix+ioPx,iz+ioPz)) ; 
                
                end
            end
        end
    end
    
  
% % update velocities from stresses
for ix=3:nx+1
    for iz = 2:nz
        U(ix,iz) = U(ix,iz)- S*rox(ix,iz)*...
                (c1*(Txx(ix,iz) - Txx(ix-1,iz)+Txz(ix,iz+1) - Txz(ix,iz)) +...
                c2*(Txx(ix+1,iz) - Txx(ix-2,iz)+Txz(ix,iz+2) - Txz(ix,iz-1)));
            
            Up(ix,iz) = Up(ix,iz)- S*rox(ix,iz)*...
                (c1*(Sp(ix,iz)-Sp(ix-1,iz))+c2*(Sp(ix+1,iz)-Sp(ix-2,iz)));
    end
end

for ix=2:nx
    for iz = 3:nz+1
        W(ix,iz) = W(ix,iz) - S*roz(ix,iz)*...
            (c1*(Tzz(ix,iz) - Tzz(ix,iz-1)  +Txz(ix+1,iz) - Txz(ix,iz)) +...
            c2*(Tzz(ix,iz+1) - Tzz(ix,iz-2)+Txz(ix+2,iz) - Txz(ix-1,iz)));
        
        Wp(ix,iz) = Wp(ix,iz) - S*roz(ix,iz)*...
            (c1*(Sp(ix,iz)-Sp(ix,iz-1))+c2*(Sp(ix,iz+1)-Sp(ix,iz-2)));
        
    end
end

%%% set_bc
[U, W, ~] = set_bc_tension(U,W);


% % %update stresses from velocities
dvvx = zeros(nz,1); dvvz = zeros(nz,1); 
for ix = 2:nx
        for iz = 2:nz
            dvvx(iz) = c1*(U(ix+1,iz) - U(ix,iz)) + ...
                c2*(U(ix+2,iz) - U(ix-1,iz));
        end
        
        for iz = 2:nz
            dvvz(iz) = c1*(W(ix,iz+1) - W(ix,iz)) + ...
                c2*(W(ix,iz+2) - W(ix,iz-1));
        end
        
        for iz = 2:nz
            Txx(ix,iz) = Txx(ix,iz) - S*l2m(ix,iz)*dvvx(iz) - S*lam(ix,iz)*dvvz(iz);
            Tzz(ix,iz) = Tzz(ix,iz) - S*l2m(ix,iz)*dvvz(iz) - S*lam(ix,iz)*dvvx(iz);
        end   
end

for ix = 3:nx+1
    for iz = 3:nz+1
        Txz(ix,iz) = Txz(ix,iz) - S*muu(ix,iz)*(c1*(U(ix,iz) - U(ix,iz-1) + ...
                W(ix,iz) - W(ix-1,iz)) + c2*(U(ix,iz+1) - U(ix,iz-2) + ...
                W(ix+1,iz) - W(ix-2,iz)));
    end
end

Sp = (1./(1 + (lam./l2m))).* (Txx + Tzz) ;

%%% set_bc
[Txx, Tzz, Txz] = set_bc_tension(Txx,Tzz,Txz);
[Up, Wp, Sp] = set_bc_tension(Up,Wp,Sp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update traces

for i=1:trace.num_traces
    Ut(i,k+1) = U(trace.ix(i),trace.iz(i));
    Wt(i,k+1) = W(trace.ix(i),trace.iz(i));
end
pp  = hypot(Up,Wp); ss = hypot((U-Up),(W-Wp));
for i=1:trace.num_traces
    PPt(i,k+1) = hypot(Up(trace.ix(i),trace.iz(i)),Wp(trace.ix(i),trace.iz(i))) ;
    SSt(i,k+1) = hypot(U(trace.ix(i),trace.iz(i)) - Up(trace.ix(i),trace.iz(i)),W(trace.ix(i),trace.iz(i)) - Wp(trace.ix(i),trace.iz(i))) ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(figure(1),'position',[112    30   790   892])
% plot
if (mod(k,40)==0)% && t>4.5)
    tic
    plot_2D_PS;
%     plot_2D_UW;
%     pause
    fprintf(['Time propogated: ',num2str(t), 's \n'])
    F = getframe(gcf);
    writeVideo(vidfile,F);
    toc
 end
 k=k+1;
 t = t+dt;
%  keyboard
end
% % 
F = getframe(gcf);
writeVideo(vidfile,F);
close(vidfile)
save_seis

