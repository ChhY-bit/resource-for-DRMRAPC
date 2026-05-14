clear,clc
%% 仿真设置
dt = 1e-4;
T = 20;

A = [0 1;-2 -3];
B = [0 ; 2];

xr_fun = @(t) [5+sin(t);...
               0+cos(t)];
ur_fun = @(x,t) (2*x(1,:) + 3*x(2,:)- sin(t))/2;
x_ini = [10;10];
%% 初始化
prob = MPCSim_init(A,B,diag([100,1]),0.1,20,0.05,[10;10],-[10;10],10,0);
N = prob.Horizen;
Ts = prob.Ts;
tspan = 0:dt:T;
N_max = length(tspan);

tspan_r = 0:dt:(T+N*Ts);
xr = xr_fun(tspan_r);
ur = ur_fun(xr,tspan_r);
x = zeros(2,N_max);
x(:,1) = x_ini;
u = zeros(1,N_max);
%% 仿真计算
dynamic_fun = @(x,u,t)A*x + B*u;
for k = 1:N_max-1
    if mod(tspan(k),Ts) == 0    % 采样周期
    % =========================== 在此处编辑控制器 ========================
        x_measured = x(:,k);
        error = xr(:,k) - x_measured;
        % 注意：必须取预测时域长度的参考值送入MPC--------------------------
        ur_part = ur(:,k:Ts/dt:k+Ts/dt*(N-1));     % 取预测时域长度，用于约束
        xr_part = xr(:,k+Ts/dt:Ts/dt:k+Ts/dt*N);   % 取预测时域长度，用于约束
        % -----------------------------------------------------------------
        u(:,k) = ur(:,k) - MPCSim_solv(prob,error,xr_part,ur_part);
        % u(:,k) = ur(:,k);
    % =====================================================================
    else                        % 零阶保持
        u(:,k) = u(:,k-1);
    end

    % 动态更新
    x(:,k+1) = update_rk4(dynamic_fun,x(:,k),u(:,k),[],dt);
end
%% 结果
figure(1)
subplot(2,1,1);
plot(tspan,x(1,:))
hold on
plot(tspan,xr(1,1:N_max),'r--');
subplot(2,1,2);
plot(tspan,x(2,:))
hold on
plot(tspan,xr(2,1:N_max),'r--');

figure(2)
plot(tspan(1:N_max-1),u(1:N_max-1))
hold on
plot(tspan(1:N_max-1),ur(1:N_max-1),'r--');