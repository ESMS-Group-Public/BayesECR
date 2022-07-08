function [y]  = SIMS(x, k, D, t)

h = k/D;
y = erfc(x./(2*sqrt(D*t)))-exp(h*x+h^2*D*t).*erfc(x./(2*sqrt(D*t))+h*sqrt(D*t));

end
