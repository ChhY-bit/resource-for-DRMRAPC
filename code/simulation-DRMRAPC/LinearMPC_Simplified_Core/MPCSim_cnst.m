function [A,Ub]=MPCSim_cnst(MPCSim_prob,err,xr_part,ur_part)
% MPCSim_cnst - Construct constraint matrices for MPC optimization problem
%
% Inputs:
%   MPCSim_prob - Struct containing MPC problem definition (output from MPCSim_init)
%       .Sx      - State prediction matrix (N*n x n)
%       .Su      - Input prediction matrix (N*n x N*p)
%       .x_ub    - Upper bounds for states (n x 1)
%       .x_lb    - Lower bounds for states (n x 1)
%       .u_ub    - Upper bounds for inputs (p x 1)
%       .u_lb    - Lower bounds for inputs (p x 1)
%   err         - Current state error (n x 1)
%   xr_part     - Reference state trajectory in future horizon (n x N)
%   ur_part     - Reference input trajectory in future horizon (p x N)
%
% Outputs:
%   A   - Constraint matrix for QP problem (2*N*p + 2*N*n x N*p)
%   Ub  - Upper bound vector for QP problem (2*N*p + 2*N*n x 1)
%
% The constraints are in the form: A * U <= Ub, where U is the control sequence
% 
% NOTES: xr_part/ur_part is not the current reference vector, 
%        but a series of vectors in predictive horizon.

%% Reading MPC Problem
Sx = MPCSim_prob.Sx;
Su = MPCSim_prob.Su;
x_max = MPCSim_prob.x_ub;
x_min = MPCSim_prob.x_lb;
u_max = MPCSim_prob.u_ub;
u_min = MPCSim_prob.u_lb;

%% Conversion
% convert to error-representation:
xe_max = xr_part-x_min;
xe_min = xr_part-x_max;
ue_max = ur_part-u_min;  
ue_min = ur_part-u_max;

% convert to extended vectors:
Xe_max = xe_max(:);
Xe_min = xe_min(:);
Ue_max = ue_max(:);
Ue_min = ue_min(:);

% convert state constraints into control constraints：
Ux_max = Xe_max-Sx*err;
Ux_min = Xe_min-Sx*err;

%% Return Standard QP Matrices
I = eye(size(Su,2));
A = [I;-I;Su;-Su];
Ub = [Ue_max;-Ue_min;Ux_max;-Ux_min];
