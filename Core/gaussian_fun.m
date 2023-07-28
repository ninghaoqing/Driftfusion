function y = gaussian_fun(coeff,t)
%UNTITLED4 此处提供此函数的摘要
%   此处提供详细说明
Amp=coeff(1);
tau=coeff(2);
t0=coeff(3);
% t=coeff(4);

y=Amp*exp(-(t-t0).^2./(2*tau^2));

end