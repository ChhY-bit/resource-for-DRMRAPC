function MPCSim_prob=MPCSim_init(Ac,Bc,Q_each,R_each,N,Ts,x_ub,x_lb,u_ub,u_lb)
% MPCSim_init - Initialize MPC problem and compute prediction matrices
%
% Inputs:
%   Ac      - Continuous-time state matrix (n x n)
%   Bc      - Continuous-time input matrix (n x p)
%   Q_each  - State weight sub-matrix for each time step (n x n)
%   R_each  - Input weight sub-matrix for each time step (p x p)
%   N       - Prediction horizon length
%   Ts      - Sampling time
%   x_ub    - Upper bounds for states (n x 1)
%   x_lb    - Lower bounds for states (n x 1)
%   u_ub    - Upper bounds for inputs (p x 1)
%   u_lb    - Lower bounds for inputs (p x 1)
%
% Outputs:
%   MPCSim_prob - Struct containing all MPC problem data
%       .Sx      - State prediction matrix (N*n x n)
%       .Su      - Input prediction matrix (N*n x N*p)
%       .Horizen - Prediction horizon length
%       .Ts      - Sampling time
%       .Ac      - Continuous-time state matrix
%       .Ad      - Discrete-time state matrix
%       .Bc      - Continuous-time input matrix
%       .Bd      - Discrete-time input matrix
%       .Q       - Extended state weight matrix (N*n x N*n)
%       .R       - Extended input weight matrix (N*p x N*p)
%       .x_ub    - Upper bounds for states (n x 1)
%       .x_lb    - Lower bounds for states (n x 1)
%       .u_ub    - Upper bounds for inputs (p x 1)
%       .u_lb    - Lower bounds for inputs (p x 1)
%
% NOTES: If no input arguments are provided, default values will be used

%% Default Parameter Settings
if nargin == 0
    % No input arguments - use all default values
    % System Model: described by Ac, Bc
    Ac = [0 1;-1 -2];
    Bc = [0;1];
    % Weight Sub-Matrices: for states
    Q_each = diag([10,1]);
    % Weight Sub-Matrices: for control
    R_each = diag(1);
    % Predictive Horizon
    N = 20;
    % Sampling Time
    Ts = 0.02;
    % Constraints: upper-bound of states
    x_ub = 10*[3;3];
    % Constraints: lower-bound of states
    x_lb = -x_ub;
    % Constraints: upper-bound of inputs
    u_ub = 2*5;
    % Constraints: lower-bound of inputs
    u_lb = -u_ub;
elseif nargin ~= 10
    error('MPCSim_init: Either provide no arguments to use defaults, or provide all 10 arguments (Ac, Bc, Q_each, R_each, N, Ts, x_ub, x_lb, u_ub, u_lb).');
end
%% Computing
% Discretization:
sys = c2d(ss(Ac,Bc,[],[]),Ts);
Ad = sys.A;
Bd = sys.B;

% Initialization:
n = size(Ad,1);     % dimension of state vector
p = size(Bd,2);     % dimension of input vector
Sx = zeros(N*n,n);      % predictive matrix for states
Su = zeros(N*n,N*p);    % predictive matrix for inputs
Q = zeros(N*n);     % state weight matrix
R = zeros(N*p);     % input weight matrix
for i = 1:N
    Sx((i-1)*n+1:i*n,:) = Ad^i;
    for j = 1:i
        Su((i-1)*n+1:i*n,(j-1)*p+1:j*p)=Ad^(i-j)*Bd;
    end
    Q((i-1)*n+1:i*n,(i-1)*n+1:i*n) = Q_each;%* ((i-N)/N+1);
    R((i-1)*p+1:i*p,(i-1)*p+1:i*p) = R_each;
end

%% Return an MPC Problem
MPCSim_prob.Sx = Sx;
MPCSim_prob.Su = Su;
MPCSim_prob.Horizen = N;
MPCSim_prob.Ts = Ts;
MPCSim_prob.Ac = Ac;
MPCSim_prob.Ad = Ad;
MPCSim_prob.Bc = Bc;
MPCSim_prob.Bd = Bd;
MPCSim_prob.Q = Q;
MPCSim_prob.R = R;
MPCSim_prob.Q_each = Q_each;
MPCSim_prob.R_each = R_each;
MPCSim_prob.x_ub = x_ub;
MPCSim_prob.x_lb = x_lb;
MPCSim_prob.u_ub = u_ub;
MPCSim_prob.u_lb = u_lb;
end
