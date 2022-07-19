%% LAUNCHER FOR SIMS
% Please see Readme and Example case of SIMS for more information

% Import data as row vectors [1xN]
x        =     ;      % depth in cm
z        =     ;      % tracer site fraction, distance corresponds to x

t        =     ;      % time associated with the SIMS dataset [1x1]


% Priors, or initial guesses at values

% k, surface reation rate
kmin     =     ; 
kmax     =     ;
k        =     ;
SIGMAk   =     ;
%D, bulk diffusion constant
Dmin     =     ;
Dmax     =     ;
D        =     ;
SIGMAD   =     ;

%SAMPLE OF VALUES FOR IG DISTRUBUTION; Experiment with your own values
ps       =     ;        % Best guess for Observational Error Standard Deviation
N        = 5000;        % Number of cycles to run the calibration
nu       = 4   ;        % Shape parameter of the ig prior
tau      = 5*ps^2;      % Scale Parameter of the ig prior

thinfact =     ;        % Thinning Factor between [0,1]

[kpost, Dpost, psipost, weights, successk, successD] = MCMCSIMSig(x, z, t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact);
