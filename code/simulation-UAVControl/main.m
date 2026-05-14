clear,clc
close all
warning off
%% 仿真参数
dt = 1e-4;
T = 15;
tspan = 0:dt:T;
N = length(tspan);
%% 控制参数
freq_pos = 20;      % 位置环（加速度信号）控制频率 (Hz)
freq_att = 100;    % 姿态环控制频率 (Hz)
freq_spd = 500;    % 角速度环控制频率 (Hz)

% 期望输入 - 时变波形 (0-5s方波, 5-10s三角波, 10-15s正弦波)
T_wave = 5;           % 每个波形周期 (s)
amp = 1;              % 幅值 ±2

ax = zeros(N,1);
ay = zeros(N,1);
az = zeros(N,1);

% 向量化生成波形
t = tspan;

% 延迟相位：ax=0, ay=1/3周期, az=2/3周期
delay = [0, T_wave/3, 2*T_wave/3];

% 0-5s: 方波 (square wave)
idx_sq = (t >= 0) & (t <= 5);
ax(idx_sq) = amp .* sign(sin(2*pi*(t(idx_sq))/T_wave));
ay(idx_sq) = amp .* sign(sin(2*pi*(t(idx_sq) - delay(2))/T_wave));
az(idx_sq) = amp .* sign(sin(2*pi*(t(idx_sq) - delay(3))/T_wave));

% 5-10s: 三角波 (triangle wave)
idx_tr = (t > 5) & (t <= 10);
phase_tr_ax = mod(t(idx_tr) - 5 - delay(1), T_wave) / T_wave;
phase_tr_ay = mod(t(idx_tr) - 5 - delay(2), T_wave) / T_wave;
phase_tr_az = mod(t(idx_tr) - 5 - delay(3), T_wave) / T_wave;
triangle_ax = 2 * abs(2 * phase_tr_ax - 1) - 1;
triangle_ay = 2 * abs(2 * phase_tr_ay - 1) - 1;
triangle_az = 2 * abs(2 * phase_tr_az - 1) - 1;
ax(idx_tr) = amp .* triangle_ax;
ay(idx_tr) = amp .* triangle_ay;
az(idx_tr) = amp .* triangle_az;

% 10-15s: 正弦波 (sine wave)
idx_si = (t > 10) & (t <= 15);
ax(idx_si) = amp .* sin(2*pi*(t(idx_si) - delay(1))/T_wave);
ay(idx_si) = amp .* sin(2*pi*(t(idx_si) - delay(2))/T_wave);
az(idx_si) = amp .* sin(2*pi*(t(idx_si) - delay(3))/T_wave);
%psi_d = [zeros(1,(N-1)/2),0.5*ones(1,(N-1)/2+1)];
psi_d = 0*sin(tspan);
%% 初始化


UAV1 = QUAV();
state = zeros(12,N);

u = zeros(4,N);     % 期望输入
% u = [2*cos(tspan);2*sin(tspan+pi/4);0.5*ones(1,N);0.5*sin(tspan)];
eta_d = zeros(3,N); % 期望姿态
force_d = zeros(1,N);   % 期望升力
w_d = zeros(3,N);   % 期望角速度
T_d = zeros(3,N);   % 期望力矩
Omega = zeros(4,N); % 电机转速
last_eta_d = zeros(3,1);
%% 仿真计算
for k = 1:N
    %% 位置环（外层输入）
    if mod(tspan(k),1/freq_pos) == 0
        % 位置环的控制输入
        u(:,k) = [ax(k);ay(k);az(k);psi_d(k)];
        % 数据准备
        m = 0.52;
        % 控制量
        [phi_d,theta_d] = QUAVS_Controller_ComputeAttitude(u(1,k),u(2,k),u(3,k),u(4,k));
        eta_d(:,k) = [phi_d;theta_d;psi_d(k)];
        force_d(k) = QUAVS_Controller_ComputeForce(u(1,k),u(2,k),u(3,k),m);
        % 计算导数
        if k ~= 1
        eta_d_dot = (eta_d(:,k) - last_eta_d)*freq_pos;
        else
            eta_d_dot = 0;
        end
        last_eta_d = eta_d(:,k);
    else
        u(:,k) = u(:,k-1);
        eta_d(:,k) = eta_d(:,k-1);
        force_d(k) = force_d(k-1);
    end
    %% 姿态环
    if mod(tspan(k),1/freq_att) == 0
        % 数据准备
        % eta_d_dot = (eta_d(:,k)-eta_d(:,max(k-1/freq_att/dt,1)))*freq_att;
        eta = [UAV1.States.phi;UAV1.States.theta;UAV1.States.psi];  % 直接测量
        % 控制量
        w_d(:,k) = QUAVS_Controller_AttitudeLoop(eta_d(:,k),eta_d_dot,eta);
    else
        w_d(:,k) = w_d(:,k-1);
    end
    %% 角速度环
    if mod(tspan(k),1/freq_spd) == 0
        % 数据准备
        J = diag([15.50e-3,15.50e-3,3.30e-2]);
        w = [UAV1.States.w_x;UAV1.States.w_y;UAV1.States.w_z];  % 直接测量
        % 控制量
        [Tx,Ty,Tz] = QUAVS_Controller_SpeedLoop(w_d(:,k),w,J);
        T_d(:,k) = [Tx;Ty;Tz];
    else
        T_d(:,k) = T_d(:,k-1);
    end
    %% 状态更新
    % 控制分配
    [omg1,omg2,omg3,omg4] = QUAVS_Controller_Allocate(force_d(k),T_d(1,k),T_d(2,k),T_d(3,k));
    if ~isreal([omg1,omg2,omg3,omg4]')
        warning('出现负转速');
        % 出现负转速的原因：加速度给定太大；内环增益过大（收敛速度要求过快）
        % 解决：==> 引入约束机制
        omg1 = real(omg1);omg2 = real(omg2);omg3 = real(omg3);omg4 = real(omg4);
    end
    Omega(:,k) = [omg1;omg2;omg3;omg4];
    % 更新
    UAV1.state_update(omg1,omg2,omg3,omg4,k,dt);
    state(:,k) = struct2array(UAV1.States)';
end
%% 结果输出
% 设置LaTeX解释器和字体
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesFontSize',12);

% 计算实际加速度（速度差分）
a_x = diff(state(1,:))/dt;
a_y = diff(state(2,:))/dt;
a_z = diff(state(3,:))/dt;
a_actual = [a_x; a_y; a_z];

% 图1: 期望加速度 vs 实际加速度 (3×1子图)
figure('Name','期望与实际加速度对比');
ylim_a = [min([u(:); a_actual(:)]), max([u(:); a_actual(:)])];
labels = {'$a_x$ (m/s$^2$)', '$a_y$ (m/s$^2$)', '$a_z$ (m/s$^2$)'};
ax_labels = {'$a_{x,d}$', '$a_{y,d}$', '$a_{z,d}$'};
a_labels = {'$a_x$', '$a_y$', '$a_z$'};
for i = 1:3
    subplot(3,1,i);
    hold on;
    % 先画实际值（实线），再画期望值（虚线），虚线在上层
    plot(tspan(1:end-1), a_actual(i,:), 'r-', 'LineWidth', 1.5, 'DisplayName', a_labels{i});
    plot(tspan, u(i,:), 'b--', 'LineWidth', 1.5, 'DisplayName', ax_labels{i});
    xlabel('Time (s)');
    ylabel(labels{i});
    legend('Location','best', 'FontSize',11);
    grid on;
    xlim([0 T]);
    ylim(ylim_a);
end

% 图2-1: 姿态环 (3×1子图)
figure('Name','姿态环响应');
ylim_eta = [min([eta_d(1,:), eta_d(2,:), eta_d(3,:), state(10,:), state(11,:), state(12,:)]), ...
            max([eta_d(1,:), eta_d(2,:), eta_d(3,:), state(10,:), state(11,:), state(12,:)])];
for i = 1:3
    subplot(3,1,i);
    hold on;
    if i == 1
        plot(tspan, state(10,:), 'r-', 'LineWidth', 1.5, 'DisplayName', '$\phi$');
        plot(tspan, eta_d(i,:), 'b--', 'LineWidth', 1.5, 'DisplayName', '$\phi_d$');
        ylabel('$\phi$ (rad)');
    elseif i == 2
        plot(tspan, state(11,:), 'r-', 'LineWidth', 1.5, 'DisplayName', '$\theta$');
        plot(tspan, eta_d(i,:), 'b--', 'LineWidth', 1.5, 'DisplayName', '$\theta_d$');
        ylabel('$\theta$ (rad)');
    else
        plot(tspan, state(12,:), 'r-', 'LineWidth', 1.5, 'DisplayName', '$\psi$');
        plot(tspan, eta_d(i,:), 'b--', 'LineWidth', 1.5, 'DisplayName', '$\psi_d$');
        ylabel('$\psi$ (rad)');
    end
    xlabel('Time (s)');
    legend('Location','best', 'FontSize',11);
    grid on;
    xlim([0 T]);
    ylim(ylim_eta);
end

% 图2-2: 角速度环 (3×1子图)
figure('Name','角速度环响应');
ylim_w = [min([w_d(1,:), w_d(2,:), w_d(3,:), state(7,:), state(8,:), state(9,:)]), ...
          max([w_d(1,:), w_d(2,:), w_d(3,:), state(7,:), state(8,:), state(9,:)])];
for i = 1:3
    subplot(3,1,i);
    hold on;
    if i == 1
        plot(tspan, state(7,:), 'r-', 'LineWidth', 1.5, 'DisplayName', '$\omega_x$');
        plot(tspan, w_d(i,:), 'b--', 'LineWidth', 1.5, 'DisplayName', '$\omega_{x,d}$');
        ylabel('$\omega_x$ (rad/s)');
    elseif i == 2
        plot(tspan, state(8,:), 'r-', 'LineWidth', 1.5, 'DisplayName', '$\omega_y$');
        plot(tspan, w_d(i,:), 'b--', 'LineWidth', 1.5, 'DisplayName', '$\omega_{y,d}$');
        ylabel('$\omega_y$ (rad/s)');
    else
        plot(tspan, state(9,:), 'r-', 'LineWidth', 1.5, 'DisplayName', '$\omega_z$');
        plot(tspan, w_d(i,:), 'b--', 'LineWidth', 1.5, 'DisplayName', '$\omega_{z,d}$');
        ylabel('$\omega_z$ (rad/s)');
    end
    xlabel('Time (s)');
    legend('Location','best', 'FontSize',11);
    grid on;
    xlim([0 T]);
    ylim(ylim_w);
end

% 图3: 电机转速 (2×2子图)
figure('Name','电机转速');
ylim_omega = [0, max(Omega(:))];
omega_names = {'$\Omega_1$', '$\Omega_2$', '$\Omega_3$', '$\Omega_4$'};
omega_colors = {[0 0.447 0.741], [0.85 0.325 0.098], [0.929 0.694 0.125], [0.494 0.184 0.556]};
rotor_labels = {'Rotor1', 'Rotor2', 'Rotor3', 'Rotor4'};
for i = 1:4
    subplot(2,2,i);
    plot(tspan, Omega(i,:), 'Color', omega_colors{i}, 'LineWidth', 1.5);
    title(rotor_labels{i}, 'FontSize',12);
    xlabel('Time (s)');
    ylabel([omega_names{i}, ' (rad/s)']);
    grid on;
    xlim([0 T]);
    ylim(ylim_omega);
end
