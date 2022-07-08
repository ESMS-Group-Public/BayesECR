%MCMCSIMS Tutorial

%I have included a data set for this tutorial.  It is called set2Array1 =
%pH2O-dependent 1.txt.  This is one of the SIMS data sets produced by Dr.
%Roger De-Souza that has already been analyzed.

%The first step is to get the data imported as row vectors.

set2Array1 = dlmread('pH2O-dependent 1.txt');

%The first column here is the set of sample thicknesses in cm. 
%The second column is the set of O18 fractions.

H220x = set2Array1(:,1)'; 
H220z = set2Array1(:,2)';

%This is from the pH20 = 220mbar experiment.

%This is actually two experiments.  We wish to only analyze one at a time,
%so the last half of the elements in x and z need to be deleted.

H220x(585:1:1171) = [];
H220z(585:1:1171) = [];


H220t = 5640;


%Now all that is needed is to enter the values we want to use for the other
%inputs.  These can be entered either in the workspace or in the command
%window.


%Note the wide bounds on k and D.
kmin = 1e-15;
kmax = 1e-2;
k = 1e-6;
SIGMAk = .1;

Dmin = 1e-15;
Dmax = 1e-2;
D = 5e-6;
SIGMAD = .1;

ps      = 0.03;
N        = 5000;
nu       = 1000;
tau      = 1001*ps^2;

thinfact = 1; 

[kpost, Dpost, psipost, weights, successk, successD] = MCMCSIMSig(H220x, H220z, H220t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact);
