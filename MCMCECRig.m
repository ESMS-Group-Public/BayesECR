function[kpost, Dpost, psipost, weights, successk, successD] = MCMCECRig(ax, ay, az, z, t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact)

%The inputs of the function are:
%
% ax          = x dimension of the rectangle 
% ay          = y dimension of the rectangle
% az          = z dimension of the rectangle
% Note: The code uses dimensionless parameters internally, so the user must
%       only use consistent units. 
% z           = Data input as a function of t. It should be the ratio of total 
%               oxygen diffused at time t to that at time infinity. This can be
%               assumed to equal the ratio of the difference of conductivities 
%               (see 1996 Yasuda and Hishinuma J. of Solid State Chem.) 
%               [MUST BE A ROW VECTOR]
% t           = Time vector from data (independant value) [MUST BE A ROW VECTOR]
% kmin        = The smallest permissible k*
% kmax        = The largest permissible k*
% k           = Initial value of k*
% SIGMAk      = Standard deviation of the proposal draws on k
% Dmin        = The smallest permissible D*
% Dmax        = The largest permissible D*
% D           = Initial value of D*
% SIGMAD      = Standard deviation of the proposal draws on D
% N           = Number of cycles to run the calibration (suggested = 100000. 
%               There is a burn-in of around 5000, so the program returns
%               the results from the last N-5000)
% nu          = shape parameter of the inverse gamma prior for the
%               observation error variance (psi)
% tau         = scale parameter of the inverse gamma prior for the
%               observation error variance
% thinfact    = thinning factor, fraction of the data you want to use in
%               the calibration, number between 0 and 1
% 
% And it returns:
% 
% kpost       = the posterior k* distribution
% Dpost       = the posterior D* distribution
% psipost     = the posterior psi distribution
% SUCCESSk    = The fraction of accepted k*
% SUCCESSD    = The fraction of accepted D*

    disp(strcat(13,' ' ))
    disp(' MCMCECRig is free software: you can redistribute it and/or modify')
    disp(' it under the terms of the GNU General Public License as published by')
    disp(' the Free Software Foundation, either version 3 of the License, or')
    disp(' (at your option) any later version.')

    %initialization and preallocation

    kpost = zeros(1, N);
    Dpost = zeros(1, N);
    psipost = zeros(1, N);
    likelies = zeros(1,N);

    kpost(1) = k;
    Dpost(1) = D;
    psipost(1) = tau/(nu+1);


    acceptedk = 0;
    acceptedD = 0;

    totacck = 0;
    totaccD = 0;

    [y] = ECRrectmodel(k, D, ax, ay, az, t);

    [LOGLIKELY, EPSILON] = LikelihoodECRrect(y, z, t, psipost(1));

    likelies(1) = LOGLIKELY;


    % Monte Carlo one-at-a-time iterations
    for i = 2:N
        if mod(i,100)==0
            disp(strcat('currently on draw',32,int2str(i),' of',32,int2str(N),' total draws'))
            disp(strcat('acceptace rate k = '))
            disp(acceptedk/100)
            disp(strcat('acceptace rate D = '))
            disp(acceptedD/100)
            if acceptedk/100 < 0.01
                SIGMAk = SIGMAk/2;
            elseif acceptedk/100 > 0.1
                SIGMAk = 1.5*SIGMAk;
            end
            if acceptedD/100 < 0.01
                SIGMAD = SIGMAD/2;
            elseif acceptedD/100 > 0.1
                SIGMAD = 1.5*SIGMAD;
            end
            totacck = totacck + acceptedk;
            totaccD = totaccD + acceptedD;
            acceptedD = 0;
            acceptedk = 0;

        end

        %proposing k and testing
        proposedk = log10(kmin) - 1;
        while (proposedk < log10(kmin) || proposedk > log10(kmax))
            proposedk = normrnd(log10(kpost(i-1)),SIGMAk);
        end
        proposedk = 10^(proposedk);

        %[proposedy] =  SIMS(x, proposedk, Dpost(i-1), t);
        [proposedy] = ECRrectmodel(proposedk, Dpost(i-1), ax, ay, az, t);

        [proposedLikely,propsilon] = LikelihoodECRrect(proposedy, z, t, psipost(i-1));

        if METROPOLIS(proposedLikely, LOGLIKELY, thinfact)
            kpost(i) = proposedk;
            LOGLIKELY = proposedLikely;
            EPSILON = propsilon;
            acceptedk = acceptedk + 1;
        else
            kpost(i) = kpost(i - 1);
        end


        %proposing D and testing
        proposedD = log10(Dmin) - 1;
        while (proposedD < log10(Dmin) || proposedD > log10(Dmax))
            proposedD = normrnd(log10(Dpost(i-1)),SIGMAD);
        end
        proposedD = 10^(proposedD);


        [proposedy] = ECRrectmodel(kpost(i), proposedD, ax, ay, az, t);

        [proposedLikely,propsilon] = LikelihoodECRrect(proposedy, z, t, psipost(i-1));


        if METROPOLIS(proposedLikely, LOGLIKELY, thinfact)
            Dpost(i) = proposedD;
            LOGLIKELY = proposedLikely;
            EPSILON = propsilon;
            acceptedD = acceptedD + 1;
        else
            Dpost(i) = Dpost(i - 1);
        end

        psipost(i) = ig(nu,tau, EPSILON);
        likelies(i) = LOGLIKELY;
        

    end

    successk = totacck/N;
    successD = totaccD/N;

    weights = (1-thinfact)*likelies;
    weights = weights - max(weights);
    weights = exp(weights);
 
end