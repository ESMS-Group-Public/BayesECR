function[LOGLIKELY, epsilon] = LikelihoodECRrect(y, z, t, psi)
%This function calculates the likelihood given the
%results from ECRrectmodel
%LOGLIKELY = 
%epsilon   = 
%
%y         =
%z         =
%t         =
%psi       =


    p = length(t);
    epsilon = z-y;
    LOGLIKELY = -(p/2)*log(psi)-(1/(2*psi))*((epsilon)*transpose(epsilon)); 
    
    
end