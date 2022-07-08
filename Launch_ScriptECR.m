%%Exract ECR data for REF 24
Matrix = readmatrix('Graphically_Extracted_dataREF_24.csv');
t = Matrix(:,1)';
z = Matrix(:,2)';


ax= 4;
ay= .15;
az= .15;

%Priors
kmin     =   1e-11;
kmax     =   1;
k        =   1e-3;
SIGMAk   =   0.1;

Dmin     =   1e-10;
Dmax     =   .01;
D        =   1e-5;
SIGMAD   =   1;

ps      = 0.02;
N        = 5000;
nu       = 1000;
tau      = 1010*ps^2;

thinfact = .15;

[kpost, Dpost, psipost,SUCCESSk, SUCCESSD, SUCCESSpsi] = MCMCECRig(ax, ay, az, z, t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact);

figure
covdraws = 50;
drawind = randi(N-1000,1,covdraws);
plot(t,z, 'r.');
hold on
for i=1:covdraws
    y = ECRrectmodel(kpost(1000+drawind(i)), Dpost(1000+drawind(i)), ax, ay, az, t);
    plot(t, y, 'LineWidth', 0.4)
end

figure
scatdraws = 1000;
scatind = randi(N-1000,1,scatdraws);
scattervecs = zeros(scatdraws, 2);
for i=1:scatdraws
    scattervecs(i,:) = [kpost(scatind(i)+1000) Dpost(scatind(i)+1000)];
end
scatter(scattervecs(:,1), scattervecs(:,2), 'MarkerFaceAlpha', 0.2)