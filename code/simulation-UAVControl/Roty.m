function res=Roty(theta)
    s_th=sin(theta);c_th=cos(theta);
    res=[c_th,0,s_th;0,1,0;-s_th,0,c_th];
end