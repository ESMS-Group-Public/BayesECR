%% LAUNCHER FOR ECR
% Please see Readme and Example case of ECR for more inforamtion

% Reduce your data to two row vectors

t        =     ;     % Time [1xN]
z        =     ;     % Row vector as a function of t [1xN]

% Dimensions of rectangular specimen
% Use consistent units
% Halve dimensions of specimen
ax       =     ;
ay       =     ;
az       =     ;

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