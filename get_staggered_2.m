%%% GIven: cp(nx+4,nz+4) cs(nx+4,nz+4) rho(nx+4,nz+4)
cp = velp; cs = vels; ro = rhog;
fac_temp = 1;

%%% for 4th order
ioXx = 3; ioXz = 2;
ioZz = 3; ioZx = 2;
ioPx = 2; ioPz = 2;
ioTx = 3; ioTz = 3;

rox = zeros(size(ro)); roz = rox;
l2m = zeros(size(cp)); lam = zeros(size(cp));
muu = zeros(size(cp));

%%Rightmost border
for iz=nz-1
    for ix = 1:nx-2
        cp2  = cp(ix,iz)*cp(ix,iz);
        cs2  = cs(ix,iz)*cs(ix,iz);
        cs2a = cs(ix+1,iz)*cs(ix+1,iz);
        cs11 = cs2*ro(ix,iz);
        cs12 = cs2*ro(ix,iz);
        cs21 = cs2a*ro(ix+1,iz);
        cs22 = cs2a*ro(ix+1,iz);
        
        if (cs11 > 0.0)
            mul  = 4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);
        else
            mul = 0.0;
        end
        mu   = cs2*ro(ix,iz);
        lamda2mu = cp2*ro(ix,iz);
        lamda    = lamda2mu - 2*mu;
        
        bx = 0.5*(ro(ix,iz)+ro(ix+1,iz));
        bz = ro(ix,iz);
        rox(ix+ioXx,iz+ioXz)=fac_temp/bx;
        roz(ix+ioZx,iz+ioZz)=fac_temp/bz;
        l2m(ix+ioPx,iz+ioPz)=fac_temp*lamda2mu;
        lam(ix+ioPx,iz+ioPz)=fac_temp*lamda;
        muu(ix+ioTx,iz+ioTz)=fac_temp*mul;
    end
end


%%% Bottom boundry
for ix = nx
    for iz=1:nz-2
        cp2  = cp(ix,iz)*cp(ix,iz);
        cs2  = cs(ix,iz)*cs(ix,iz);
        cs2b = cs(ix,iz+1)*cs(ix,iz+1);
        cs11 = cs2*ro(ix,iz);
        cs12 = cs2b*ro(ix,iz+1);
        cs21 = cs2*ro(ix,iz);
        cs22 = cs2b*ro(ix,iz+1);
        
        if (cs11 > 0.0)
            mul  = 4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);
            
        else
            mul = 0.0;
        end
        
        mu   = cs2*ro(ix,iz);
        lamda2mu = cp2*ro(ix,iz);
        lamda    = lamda2mu - 2*mu;
        
        bx = ro(ix,iz);
        bz = 0.5*(ro(ix,iz)+ro(ix,iz+1));
        rox(ix+ioXx,iz+ioXz)=fac_temp/bx;
        roz(ix+ioZx,iz+ioZz)=fac_temp/bz;
        l2m(ix+ioPx,iz+ioPz)=fac_temp*lamda2mu;
        lam(ix+ioPx,iz+ioPz)=fac_temp*lamda;
        muu(ix+ioTx,iz+ioTz)=fac_temp*mul;
    end
end

%%% corner node
for ix = nx-1
    for iz = nz-1

		cp2  = cp(ix,iz)*cp(ix,iz);
		cs2  = cs(ix,iz)*cs(ix,iz);
		mu   = cs2*ro(ix,iz);
		lamda2mu = cp2*ro(ix,iz);
		lamda    = lamda2mu - 2*mu;
		bx = ro(ix,iz);
		bz = ro(ix,iz);
		rox(ix+ioXx,iz+ioXz)=fac_temp/bx;
		roz(ix+ioZx,iz+ioZz)=fac_temp/bz;
		l2m(ix+ioPx,iz+ioPz)=fac_temp*lamda2mu;
		lam(ix+ioPx,iz+ioPz)=fac_temp*lamda;
		muu(ix+ioTx,iz+ioTz)=fac_temp*mu;
    end
end


%%% Rest
for ix = 1:nx-2
    for iz = 1:nz-2
        cp2  = cp(ix,iz)*cp(ix,iz);
        cs2  = cs(ix,iz)*cs(ix,iz);
        cs2a = cs(ix+1,iz)*cs(ix+1,iz);
        cs2b = cs(ix,iz+1)*cs(ix,iz+1);
        cs2c = cs(ix+1,iz+1)*cs(ix+1,iz+1);
        
        cs11 = cs2*ro(ix,iz);
        cs12 = cs2b*ro(ix,iz+1);
        cs21 = cs2a*ro(ix,iz);
        cs22 = cs2c*ro(ix,iz+1);
        
        if (cs11 > 0.0)
            mul  = 4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);
        else
            mul = 0.0;
        end
        
        mu   = cs2*ro(ix,iz);
        lamda2mu = cp2*ro(ix,iz);
        lamda    = lamda2mu - 2*mu;
        
        bx = 0.5*(ro(ix,iz)+ro(ix+1,iz));
        bz = 0.5*(ro(ix,iz)+ro(ix,iz+1));
        rox(ix+ioXx,iz+ioXz)=fac_temp/bx;
        roz(ix+ioZx,iz+ioZz)=fac_temp/bz;
        l2m(ix+ioPx,iz+ioPz)=fac_temp*lamda2mu;
        lam(ix+ioPx,iz+ioPz)=fac_temp*lamda;
        muu(ix+ioTx,iz+ioTz)=fac_temp*mul;
        
    end
end

