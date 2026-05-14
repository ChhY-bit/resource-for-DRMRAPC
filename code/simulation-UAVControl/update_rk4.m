function x_new = update_rk4(dynamic_fun,x,u,t,dt)
K1 = dynamic_fun(x,u,t);
K2 = dynamic_fun(x+K1*dt/2,u,t+dt/2);
K3 = dynamic_fun(x+K2*dt/2,u,t+dt/2);
K4 = dynamic_fun(x+K3*dt,u,t+dt);
x_new = x+dt*(K1+2*K2+2*K3+K4)/6;
end