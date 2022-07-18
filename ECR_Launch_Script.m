%% LAUNCHER FOR ECR
% Please see Readme and Example case of ECR for more inforamtion

% Reduce your data to two row vectors

t        =     ;     % Time [1xN]
z        =     ;     % Row vector as a function of t [1xN]

% Dimensions of rectangular speciman
% Use consistent units
ax       =     ;
ay       =     ;
az       =     ;

% Priors, or initial guesses at values
% k, surface reation rate
kmin     =     ;
kmax     =     ;
SIGMAk   =     ;
%D, bulk diffusion constant
Dmin     =     ;
Dmax     =     ;
D        =     ;
SIGMAD   =     ;

%SAMPLE OF VALUES FOR IG DISTRUBUTION; Experiment with your own values
ps       = 0.02;        % StD of Obs. Variance
N        = 5000;        % Number of cycles to run the calibration
nu       = 1000;        % Shape parameter of the ig prior
tau      = 1010*ps^2;   % Scale Parameter of the ig prior

thinfact =     ;        % Thinning Factor between [0,1]

[kpost, Dpost, psipost,SUCCESSk, SUCCESSD, SUCCESSpsi] = MCMCECRig(ax, ay, az, z, t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact);

% Plots
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