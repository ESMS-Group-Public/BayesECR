%%Exract ECR data for REF 24
Matrix = readmatrix('Graphically_Extracted_dataREF_24.csv');
t = Matrix(:,1)';
z = Matrix(:,2)';

%Dimensions of cells, in cm
ax= 4;
ay= .15;
az= .15;

%Priors
kmin     =   1e-6;
kmax     =   .1;
k = 1e-5;
SIGMAk   =   .1;

Dmin     =   1e-6 ;
Dmax     =   .1;
D        =   1e-6;
SIGMAD   =   .1;

psimin   =   0.00001;
psimax   =   .1;
psi      =   0.01;
SIGMApsi =   0.1;

BT       = 1;
N        = 1000;
xbins    = 10;
ybins    = 10;
nu       = 1;
tau      = 1;
kbins    = 10;
Dbins    = 10;
psibins  = 10;
thinfact = .1;

[kpost, Dpost, psipost,SUCCESSk, SUCCESSD, SUCCESSpsi] = MCMCECRig(ax, ay, az, z, t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact);
covdraws = 50;
drawind = randi(N-100,1,covdraws);
plot(t,z, 'r.');
hold on
for i=1:covdraws
    y = ECRrectmodel(kpost(100+drawind(i)), Dpost(100+drawind(i)), ax, ay, az, t);
    plot(t, y, 'LineWidth', 0.4)
end