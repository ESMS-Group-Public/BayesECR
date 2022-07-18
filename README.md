# BayesECR Guide
A Matlab Bayesian Fitting of Parameters in Electrical Conductivity Relaxation (ECR) and Isotope Exchange/Secondary Ion Mass Spectrometry (SIMS) Experiments [^1].
This guide will help explain how the Bayesian functions perform their calculations and how to implement your own data to receive results.

## Getting started
Determine which analysis is appropriate for your experiment. There are currently two different analyses that can be run. These are the MCMCECRig and MCMCSIMSig. ECR and SIMS are two different ways to measure the same parameters. ECR measures conductivity as a function of time after a step change in the gas concentration, while SIMS measures isotope concentration as a function of distance at a certain exposure time. For more information and examples regarding each analysis refer to the corresponding section. After determining which analysis you are going to run open the related launch script. The accuracy of the analysis will vary drastically depending on what is inputted into the priors/initial guesses sections. We have supplied some generalized inputs for the IG distribution (Inverse Gamma). These inputs will work on most problems and feel free to change these values if needed.

## Definition of variables used
#### Data Reduction variables
```
t:   Time [1xN]
z:   Row vector as a function of t [1xN]
```

#### Priors/initial guesses variables
```
kmin:      Minimum of the surface reaction rate
kmax:      Maximum of the surface reaction rate
k:         Surface reaction rate
SIGMAk:    Standard Deviation of the surface reaction rate

Dmin:      Minimum of the bulk diffusion constant
Dmax:      Maximum of the bulk diffusion constant
D:         Bulk diffusion constant
SIGMAD:    Standard Deviation of the bulk diffusion constant
```

#### Values for Inverse Gamma (IG) distribution
```
ps:         Standard Deviation of Observation Variance
N:          Number of cycles to run
nu          Shape parameter of the IG prior
tau:        Scale Parameter of the IG prior
thinfact:   Thinning Factor reduces data must be between 0 and 1
```
#### Specific Variables for ECR analysis only
```
Dimensions of rectangular specimen:
ax:    x dimension of specimen
ay:    y dimension of specimen
az:    z dimension of specimen
t:   Time [1xN]
z:   Row vector as a function of t [1xN]
```

#### Specific Variables for SIMS analysis only
```
x:    Depth in cm
z:    Tracer site fraction, distance corresponds to x

t     Time associated with the SIMS dataset [1x1]
```


## ECR analysis
Inorder to conduct the ECR analysis ECR_Launch_Script.m should be ran. First thing to do is input your time data into the time variable, t. The z variable is the variable you are looking to test; this should be a function of the time variable.
```
t        =     ;     % Time [1xN]
z        =     ;     % Row vector as a function of t [1xN]
```
Below this you should enter your dimensions of the specimen.
```
ax       =     ;
ay       =     ;
az       =     ;
```
After this you need to enter the initial guess. There are two sets of initial guesses. One has to do with k the surface reaction rate and D the bulk diffusion constant.
```
kmin     =     ;
kmax     =     ;
k        =     ;
SIGMAk   =     ;

Dmin     =     ;
Dmax     =     ;
D        =     ;
SIGMAD   =     ;
```
We have supplied variables from our own testing that seem to work on most problems for the Inverse Gamma distribution. Please change these values if needed.
A good rule when testing to produce more accurate results would be first to increase the value of N then change the other variables to tweak the accuracy. The thinfact variable allows you to reduce the amount of data tested this will decrease the runtime depending on the value. The value must be between 0 and 1.
```
ps       = 0.02;        
N        = 5000;        
nu       = 1000;        
tau      = 1010*ps^2;   

thinfact =     ;
```

The launch script for ECR currently plots the ECRrectmodel vs time and a scatter plot for the posterior distribution of k*  and D*. 

## SIMS analysis
Inorder to conduct the SIMS analysis SIMS_Launch_Script.m should be ran. Input your depth in centimeters into the x variable, the z variable should be a distance that is function of x, and t should be the time associated with the dataset.
```
x        =     ;      % depth in cm
z        =     ;      % tracer site fraction, distance corresponds to x

t        =     ;      % time associated with the SIMS dataset [1x1]
```


After this you need to enter the initial guess. There are two sets of initial guesses. One has to do with k the surface reaction rate and D the bulk diffusion constant.
```
kmin     =     ;
kmax     =     ;
k        =     ;
SIGMAk   =     ;

Dmin     =     ;
Dmax     =     ;
D        =     ;
SIGMAD   =     ;
```
We have supplied variables from our own testing that seem to work on most problems for the Inverse Gamma distribution. Please change these values if needed.
A good rule when testing to produce more accurate results would be first to increase the value of N then change the other variables to tweak the accuracy. The thinfact variable allows you to reduce the amount of data tested this will decrease the runtime depending on the value. The value must be between 0 and 1.
```
ps       = 0.02;        
N        = 5000;        
nu       = 1000;        
tau      = 1010*ps^2;   

thinfact =     ;
```




## License

BSD 2-Clause License

Copyright (c) 2022, ESMS-Group-Public
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

[^1]: A Bayesian approach to electrical conductivity relaxation and isotope exchange/secondary ion mass spectrometry, https://www.sciencedirect.com/science/article/abs/pii/S0167273814005177.

