%% Example of SIMS use cited as REF 23
% Original paper that details experiement can be found at
% https://onlinelibrary.wiley.com/doi/abs/10.1002/fuce.201300087

% A dataset is included.  It is called set2Array1 =
% pH2O-dependent 1.txt.  This is one of the SIMS data sets produced by Dr.
% Roger De-Souza that has already been analyzed.
% This is from the pH20 = 220mbar experiment.

% Import data as row vectors
set2Array1 = dlmread('pH2O-dependent 1.txt');

% The first row here is the set of sample thicknesses in cm. 
% The second row is the set of O18 fractions.

H220x = set2Array1(:,1)'; 
H220z = set2Array1(:,2)';


% This is actually two experiments.  We wish to only analyze one at a time,
% so the last half of the elements in x and z need to be deleted.

H220x(585:1:1171) = [];
H220z(585:1:1171) = [];


H220t = 5640;


% Now all that is needed is to enter the values we want to use for the other
% inputs.  These can be entered either in the workspace or in the command
% window.


%Note the wide bounds on k and D.
% k, surface reaction rate
kmin     = 1e-15;      % Smallest permissible k*
kmax     = 1e-2;       % Largest permissible k*
k        = 1e-6;       % Initial value of k*
SIGMAk   = .1;         % Standard deviation of the draws on k
%D, bulk diffusion constant
Dmin     = 1e-15;      % Smallest permissible D*
Dmax     = 1e-2;       % Largest permissible D*
D        = 5e-6;       % Initial value of D*
SIGMAD   = .1;         % Standard deviation of the draws on D

ps       = 0.03;       % Standard Deviation of Observational error
N        = 5000;       % Number of cycles to run the calibration
nu       = 1000;       % Shape parameter of the inverse gamma prior
tau      = 1001*ps^2;  % Scale parameter of the inverse gamma prior

thinfact = 1;          % Thinning factor b/w [0,1] fraction of data to be used for calibration

[kpost, Dpost, psipost, weights, successk, successD] = MCMCSIMSig(H220x, H220z, H220t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact);
