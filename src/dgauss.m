% creates a matrix containing the first derivative of the gaussian
% funtion, using the specified sigma
function dG = dgauss(sig)

x = floor(-3*sig):ceil(3*sig);
G = exp(-0.5*x.^2/sig^2);
G = G/sum(G);
dG = -x.*G/sig^2;