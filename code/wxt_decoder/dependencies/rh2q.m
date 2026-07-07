function q=rh2q(t,p,hrh,whoknows)
% q=rh2q(atmp,bpr,hrh) ...  % from S. Bigorre with correction

% compute specific humidity in g/kg from relative humidity
% uses Teten's formula from Buck (1981, JAM 20, 1527-1532)
%  input:  t   air temperature in degC
%          p   barometric pressure in mbars
%  output: q   specific humidity in g/kg
% S. Bigorre 01/10/2008

%First, calculate saturated vapor pressure es in mb, using Buck:
es=6.1121.*exp(17.502.*t./(t+240.97)).*(1.0007+3.46e-6*p);

%Second, by definition of relative humidity, vapor pressure is:
e=es.*hrh/100; %where HRH is in % (typical value 70)

%Finally, compute specific humidity using:
q=.62197*e./(p-.378*e); % q in kg/kg
q=q*1000; % q in g/kg

