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
N        = 50000;         % Number of cycles to run the calibration
nu       = 4;         % Shape parameter of the inverse gamma prior
tau      = 5*ps^2;    % Scale parameter of the inverse gamma prior

thinfact = .15;          % Thinning factor b/w (0,1) fraction of data to be used for calibration

[kpost, Dpost, psipost,weights, SUCCESSk, SUCCESSD] = MCMCECRig(ax, ay, az, z, t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact);

% Plots
% Hardcoded subtraction of 1000 ( NOT to be confused with scatdraws ) is
% used to reduce burn-in
% A rule-of-thumb is that by eliminating at least 1000 samples, burn-in can
% be reduced
figure % Coverage Plot
covdraws = 50;
drawind = randi(N-1000,1,covdraws);
plot(t,z, 'r.');
hold on
for i=1:covdraws
    y(i,:) = ECRrectmodel(kpost(1000+drawind(i)), Dpost(1000+drawind(i)), ax, ay, az, t);
end
cut = floor(length(t)*0.025);     %For 95% confidence interval
nn = covdraws;
for i=1:length(t)
    points_at_t = y(:,i);
    for j = 1:cut
        [~,ind] = max(points_at_t);
        points_at_t(ind) = NaN;
        [~,ind] = min(points_at_t);
        points_at_t(ind) = NaN;
    end
    points_at_t = [points_at_t ones(covdraws,1)*t(i)];
    y_cut = rmmissing(points_at_t);
    if i == 1
        all_y = y_cut(:,1);
        all_t = y_cut(:,2);
    else
        all_y = [all_y;y_cut(:,1)];
        all_t = [all_t;y_cut(:,2)];
    end
end
bounds = boundary(all_t,all_y);
plot(all_t(bounds),all_y(bounds))
fill(all_t(bounds),all_y(bounds),'y','FaceAlpha',0.3)
ylabel("Normalized Conductivity");
xlabel("Time (s)");

figure
histdraws = N; % Histogram pulls from all the posterior
histind = randi(N-1000,1,histdraws);
histvecs = zeros(histdraws, 2);
for i=1:histdraws
    histvecs(i,:) = [kpost(histind(i)+1000) Dpost(histind(i)+1000)];
end
posts = [kpost' Dpost'];
histogram2(histvecs(:,1), histvecs(:,2),100,'DisplayStyle','tile')%'MarkerFaceAlpha', 0.2)
ylabel("D");
xlabel("k");