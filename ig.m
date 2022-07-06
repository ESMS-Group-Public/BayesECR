function psi = ig(nu,tau, epsilon)
%ig calculates the inverse gamma and returns:
%psi     = 
%
%nu      = 
%tau     = 
%epsilon = 



  n   = length(epsilon);
  nu1 = nu + n/2;
  tau1 = tau + norm(epsilon)^2/2;
  psi = 1/gamrnd(nu1, 1/tau1);

  
end