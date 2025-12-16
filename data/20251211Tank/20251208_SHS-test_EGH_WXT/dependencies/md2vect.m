function [vx,vy] = md2vect(mag,dir);

%function [vx,vy] = md2vect(mag,dir);
% converts speed and heading to vectors east and north (vx,vy)

jk=find(dir > 360);
dir(jk)=dir(jk) - 360;

 jk=find(dir<0);
 if(length(jk)>0)
       dir(jk)=(dir(jk))+360;
 end
%
% to radians
dir=dir*pi/180;

vx = mag .* sin(dir);
vy = mag .* cos(dir);

