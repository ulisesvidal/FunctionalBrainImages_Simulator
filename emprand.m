function xr = emprand(len,x,cdf)
%EMPRAND Generates a random vector of length LEN according to the
%statistical distribution specified in CDF. 

% Generate uniform random number between 0 and 1
ur =  rand(len,1);

% Interpolate ur from empirical cdf and extraplolate for out of range
% values.
xr = interp1(cdf,x,ur,'linear','extrap'); 


