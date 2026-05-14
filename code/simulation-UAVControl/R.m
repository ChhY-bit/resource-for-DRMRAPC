function res = R(phi,theta,psi)
    res = Rotz(psi)*Roty(theta)*Rotx(phi);
end