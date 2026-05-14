function x_new = update_rk4(dynamic_fun,x,u,t,dt)
% update_rk4 - Perform one-step integration using 4th-order Runge-Kutta method
%
% Inputs:
%   dynamic_fun - Function handle for system dynamics: x_dot = dynamic_fun(x, u, t)
%   x           - Current state vector (n x 1)
%   u           - Input vector (p x 1)
%   t           - Current time
%   dt          - Time step for integration
%
% Outputs:
%   x_new       - Next state vector after one integration step (n x 1)
%
% NOTES: Implements classic RK4 integration scheme
K1 = dynamic_fun(x,u,t);
K2 = dynamic_fun(x+K1*dt/2,u,t+dt/2);
K3 = dynamic_fun(x+K2*dt/2,u,t+dt/2);
K4 = dynamic_fun(x+K3*dt,u,t+dt);
x_new = x+dt*(K1+2*K2+2*K3+K4)/6;
end