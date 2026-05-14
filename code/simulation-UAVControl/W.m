function [res,in]=W(phi,theta,psi)
    res = [1 sin(phi)*tan(theta) cos(theta)*tan(theta);...
           0 cos(phi) -sin(phi);
           0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
    in = [1 0 -sin(theta);0 cos(phi) sin(phi)*cos(theta);0 -sin(phi) cos(phi)*cos(theta)];
end