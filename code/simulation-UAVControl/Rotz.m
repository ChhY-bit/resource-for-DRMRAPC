function res = Rotz(psi)
    s_th=sin(psi);c_th=cos(psi);
    res=[c_th,-s_th,0;s_th,c_th,0;0,0,1];
end