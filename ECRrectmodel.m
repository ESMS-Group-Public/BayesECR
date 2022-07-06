function [y] = ECRrectmodel(k, D, ax, ay, az, t)
% Solves the ECR rectangle model as described in Yasuda and Hishinuma's 
% 1996 paper in J. Solid State Chem. 
% Returns a vector y = M(t)/M(inf)
% k = effective exchange coefficient
% D = Effective Diffusion coefficient
% ax, ay, az are the x, y, and z dimensions of the rectangle
% t = is the time vector


%D =  bulk diffusion coeff.
%k =  surface exchange coeff.
%ax =.05; % half of the sample length along x-axis.
%k = 5e-5;
%D = 5e-5;
%N = length(times);
%for n = 1:N;
   %t(n) = times(n);
%end   

%[t,z] = ECRDATAFUNCTION;

Lx = ax*k/D;

TERMx = 1;
SUMx = 1;

box = fzero(@(bn)bn*tan(bn)-Lx,[-.00000000001, (pi/2)-.000000001]);
TERMx = 2*Lx^2*exp(-(box^2)*D.*(t)/(ax^2))/(box^2*(box^2+Lx^2+Lx));
 
SUMx = TERMx;

ix = 2;
while(max(TERMx./SUMx) > 0.0001)
  bn = fzero(@(bn)bn*tan(bn)-Lx,[(ix-1)*pi-.00000000001, (2*(ix-1)+1)*(pi/2)-.000000001]);
  TERMx = 2*Lx^2*exp(-(bn^2)*D.*(t)/(ax^2))/(bn^2*(bn^2+Lx^2+Lx));
  SUMx = SUMx + TERMx;
  ix = ix+1;
end

%ay = .05; %Enter half of sample length along y-axis.  
Ly = ay*k/D;
TERMy = 1;
SUMy = 1;

boy = fzero(@(bn)bn*tan(bn)-Ly,[-.00000000001, (pi/2)-.000000001]);
TERMy = 2*Ly^2*exp(-(boy^2)*D.*(t)/(ay^2))/(boy^2*(boy^2+Ly^2+Ly));

SUMy = TERMy;

iy = 2;
while(max(TERMy./SUMy) > 0.0001)
  bny = fzero(@(bn)bn*tan(bn)-Ly,[(iy-1)*pi-.00000000001, (2*(iy-1)+1)*(pi/2)-.000000001]);
  TERMy = 2*Ly^2*exp(-(bny^2)*D.*(t)/(ay^2))/(bny^2*(bny^2+Ly^2+Ly));
  SUMy = SUMy + TERMy;
  iy = iy+1;
end


%az = .05; %Enter half of sample length along z-axis.  
Lz = az*k/D;

TERMz = 1;
SUMz = 1;

boz = fzero(@(bn)bn*tan(bn)-Lz,[-.00000000001, (pi/2)-.000000001]);
TERMz = 2*Lz^2*exp(-(boz^2)*D.*(t)/(az^2))/(boz^2*(boz^2+Lz^2+Lz));

SUMz = TERMz;

iz = 2;
while(max(TERMz./SUMz) > 0.0001)
  bnz = fzero(@(bn)bn*tan(bn)-Lz,[(iz-1)*pi-.00000000001, (2*(iz-1)+1)*(pi/2)-.000000001]);
  TERMz = 2*Lz^2*exp(-(bnz^2)*D.*(t)/(az^2))/(bnz^2*(bnz^2+Lz^2+Lz));
  SUMz = SUMz + TERMz;
  iz = iz+1;
end

SUM3D = SUMx.*SUMy.*SUMz;
y = 1-SUM3D;
end