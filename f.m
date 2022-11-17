function y = f(psi)
% N=5:1:10;
% N=5;
global N;
AF_norm = sin(N*(psi/2))/(N*sin(psi/2));
y = 20*log10(AF_norm)+11;