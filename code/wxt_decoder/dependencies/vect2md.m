function [mag,dir] = vect2md(vx,vy);

%function [mag,dir] = vect2md(vx,vy);
% converts vectors east and north (vx,vy) to speed and heading

mag = sqrt(vx.^2 + vy.^2);

dir = atan2(vx,vy)*180/pi;

jk=find(dir<0);

if(length(jk)>0)
      dir(jk)=360+(dir(jk));
end
