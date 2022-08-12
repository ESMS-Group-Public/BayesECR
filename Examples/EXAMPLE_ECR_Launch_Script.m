%% Example of ECR use cited as REF 25 
% Original paper that details experiement can be found at
% https://www.sciencedirect.com/science/article/abs/pii/S0167273800005543

% Data Extraction
% Data is formatted to two row vectors
Matrix = readmatrix('Graphically_Extracted_data[REF_24].csv');
t = Matrix(:,1)';       % Time Row Vector
z = Matrix(:,2)';       % Row Vector as a function of t

% Dimensions of rectangular speciman, in this case LSCF
% Use consistent units, in this case cm
ax= 4;
ay= .15;
az= .15;
% FOR ALL CASES, dimensions should be halved
ax = ax/2;
ay = ay/2;
az = az/2;
% Priors, or initial guesses at values
% since minimum values are small, guesses can be an order of magnitude
% k, surface reaction rate
kmin     =   1e-11;      % Smallest permissible k*
kmax     =   1;          % Largest permissible k*
k        =   1e-3;       % Initial value of k*
SIGMAk   =   0.1;        % Standard Deviation of the draws on k
%D, bulk diffusion constant
Dmin     =   1e-10;      % Smallest permissible D*
Dmax     =   .01;        % Largest permissible D*
D        =   1e-5;       % Initial value of D*
SIGMAD   =   1;          % Standard Deviation of the draws on D

ps       = 0.02;         % Standard Deviation of Observational error
N        = 5000;         % Number of cycles to run the calibration
nu       = 1000;         % Shape parameter of the inverse gamma prior
tau      = 1010*ps^2;    % Scale parameter of the inverse gamma prior

thinfact = .15;          % Thinning factor b/w (0,1) fraction of data to be used for calibration

[kpost, Dpost, psipost,SUCCESSk, SUCCESSD, SUCCESSpsi] = MCMCECRig(ax, ay, az, z, t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact);

% Plots
% Hardcoded subtraction of 1000 ( NOT to be confused with scatdraws ) is
% used to reduce burn-in
% A rule-of-thumb is that by eliminating at least 1000 samples, burn-in can
% be reduced
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