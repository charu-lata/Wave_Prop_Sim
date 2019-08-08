function [At, Bt, Ct] = set_bc_tension(A,B,C,w)
%%% just a 2-D taper
%%w =  fraction of sides to be tapered wrt complete length

if nargin==2
    C = zeros(size(A));
    w = 0.1;
elseif nargin==3
    w = 0.1;
end
%%Initialise
At = zeros(size(A)); Bt = zeros(size(B)); Ct = zeros(size(C)); 

[n,m,p] =size(A);
taper_matrix = repmat(tukeywin(n,w),1,m).*rot90(repmat(tukeywin(m,w),1,n));

for i = 1:p
    At(:,:,i) = A(:,:,i).*taper_matrix;
    Bt(:,:,i) = B(:,:,i).*taper_matrix;
    Ct(:,:,i) = C(:,:,i).*taper_matrix;
    
end


end



% % side wall-1
% [nx,nz,~] = size(Txx);
% pp1 = Txx(1,2:nz,2)+Txx(2,2:nz,2) - Txx(2,2:nz,1);
% pp2 = U(1,2:nz);		% note that vel() is scaled with dr/dt in input_model
% pp3 = Txx(2,2:nz,2) - Txx(1,2:nz,2) - Txx(3,2:nz,1) + Txx(2,2:nz,1);
% Txx(1,2:nz,3) = pp1+pp2.*pp3;
% ppp1 = Tzz(1,2:nz,2)+Tzz(2,2:nz,2) - Tzz(2,2:nz,1);
% ppp2 = W(1,2:nz);		% note that vel() is scaled with dr/dt in input_model
% ppp3 = Tzz(2,2:nz,2) - Tzz(1,2:nz,2) - Tzz(3,2:nz,1) + Tzz(2,2:nz,1);
% Tzz(1,2:nz,3) = ppp1-ppp2.*ppp3;
% 
% % side wall-2:
% pp1 = Txx(nx,2:nz,2)+Txx(nx-1,2:nz,2) - Txx(nx-1,2:nz,1);
% pp2 = U(nx,2:nz);		% note that vel() is scaled with dr/dt in input_model
% pp3 = Txx(nx,2:nz,2) - Txx(nx-1,2:nz,2) - Txx(nx-2,2:nz,1) + Txx(nx-1,2:nz,1);
% Txx(nx,2:nz,3) = pp1-pp2.*pp3;
% ppp1 = Tzz(1,2:nz,2)+Tzz(2,2:nz,2) - Tzz(2,2:nz,1);
% ppp2 = W(1,2:nz);		% note that vel() is scaled with dr/dt in input_model
% ppp3 = Tzz(2,2:nz,2) - Tzz(1,2:nz,2) - Tzz(3,2:nz,1) + Tzz(2,2:nz,1);
% Tzz(nx,2:nz,3) = ppp1-ppp2.*ppp3;
% 
% %top
% Txx(1:nx,1,3) = 0.0;
% Txz(1:nx,1,3) = 0.0;
% Tzz(1:nx,1,3) = 0.0;
% 
% % bottom
% pp1 = Txx(1:nx,nz,2) + Txx(1:nx,nz-1,2) - Txx(1:nx,nz-1,1);
% pp2 = U(1:nx,nz);
% pp3 = Txx(1:nx,nz,2) - Txx(1:nx,nz-1,2) - Txx(1:nx,nz-1,1) + Txx(1:nx,nz-2,1);
% Txx(1:nx,nz,3) = pp1-pp2.*pp3;
% ppp1 = Tzz(1:nx,nz,2) + Tzz(1:nx,nz-1,2) - Tzz(1:nx,nz-1,1);
% ppp2 = W(1:nx,nz);
% ppp3 = Tzz(1:nx,nz,2) - Tzz(1:nx,nz-1,2) - Tzz(1:nx,nz-1,1) + Tzz(1:nx,nz-2,1);
% Tzz(1:nx,nz,3) = ppp1-ppp2.*ppp3;

