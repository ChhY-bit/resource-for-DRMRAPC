function [u,u_series]=MPCSim_solv(MPCSim_prob,err,xr,ur)
% MPCSim_solv - Solve MPC optimization problem to compute optimal control sequence
%
% Inputs:
%   MPCSim_prob - Struct containing MPC problem definition (output from MPCSim_init)
%       .Sx      - State prediction matrix (N*n x n)
%       .Su      - Input prediction matrix (N*n x N*p)
%       .Q       - Extended state weight matrix (N*n x N*n)
%       .R       - Extended input weight matrix (N*p x N*p)
%       .Bc      - Continuous-time input matrix (n x p)
%   err         - Current state error (n x 1), defined as x_ref - x_current
%   xr          - Reference state trajectory in future horizon (n x N)
%   ur          - Reference input trajectory in future horizon (p x N)
%
% Outputs:
%   u           - First optimal control action (p x 1)
%   u_series    - Complete optimal control sequence over horizon (N*p x 1)
%
% NOTES: Uses OSQP solver for quadratic programming
%        xr/ur is not the current reference vector,
%        but a series of vectors in predictive horizon.

%% Reading MPC Problem
Sx = MPCSim_prob.Sx;
Su = MPCSim_prob.Su;
Q = MPCSim_prob.Q;
R = MPCSim_prob.R;

%% QP Construction
% Hessian Matrices:
H = Su'*Q*Su+R;
H = (H'+H)/2;   % Symmetrization to avoid floating-point errors
% Coefficient Vector:
F = err'*Sx'*Q*Su;

%% Create Constraints
[A,Ub]=MPCSim_cnst(MPCSim_prob,err,xr,ur);

%% Solve QP
% using quadprog:
% options = optimoptions('quadprog','Display', 'off');
% u_series = quadprog(H,F',A,Ub,[],[],[],[],[],options);    % with constraints
% u_series = quadprog(H,F',[],[],[],[],[],[],[],options);   % without constraints

% using OSQP:
prob = osqp;
options = osqp.default_settings;
options.verbose = false;
prob.setup(H,F,A,[],Ub,options);    % with constraints
% prob.setup(H,F,[],[],[],options);    % without constraints
sol = prob.solve();
u_series = sol.x;

%% Return the Solution
p = size(MPCSim_prob.Bc,2);
if any(isnan(u_series(1:p)))
    u = zeros(p,1);
    warning("The optimization problem has no solution!")
else
    u = u_series(1:p);
end
end