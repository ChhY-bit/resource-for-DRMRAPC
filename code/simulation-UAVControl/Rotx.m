function res=Rotx(phi)
    s_th=sin(phi);c_th=cos(phi);
    res=[1,0,0;0,c_th,-s_th;0,s_th,c_th];
end