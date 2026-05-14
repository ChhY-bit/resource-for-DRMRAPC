%% 基于扰动补偿的自适应模型预测控制
%
%  作者：杨宸涵
%  学校：北京理工大学
%  学号：1120222643
%  专业：自动化
%  学院：自动化学院
%  指导教师：孙中奇
%  日期：2026年4月

clear,clc,close all,warning off
addpath("LinearMPC_Simplified_Core")
addpath("Utils")

%% 1 仿真参数
dt = 1e-3;  % 仿真步长
T = 30;     % 仿真时长
expr_names = {'fixed_point','eight_shape','spiral_shape','square_shape'};
%% 2 控制参数
% 2.1 ============================== 控制模式 ==============================
% 0: PID
% 1: MPC
% 2: DR-MPC
% 3: MRAPC
% 4: DR-MRAPC
mode = 4;
mode_names = {'PID', 'MPC', 'DR-MPC', 'MRAPC', 'DR-MRAPC'};

% 2.2 ============================== 参数设置 ==============================
Ts = 0.02;      % 控制周期
umax = 10*[2;2;2]; % 最大控制量
umin = -umax;   % 最小控制量
xmax = 10*[2;1;1;1;2.2;1];           % 最大状态量
% xmin = [-3;-2;-3;-2;-0.01;-3];    % 最小状态量
xmin = -xmax;
x_ini = [0;0;0;0;0;0];  % 初始状态


    % 2.2.1 -------------------- PID参数 --------------------
    Kp = [0.5;0.5;0.5];           % 比例系数
    Ki = [0.0;0.0;0.0];     % 积分系数
    Kd = [5;5;1];           % 微分系数
    % 2.2.2 -------------------- MPC参数 --------------------
    Q = diag([10,1,10,1,10,1]);     % 误差权重
    R = diag([1,1,1]);           % 控制权重
    horizen = 50;                   % 预测时域
    % 2.2.3 -------------------- MRAC参数 --------------------
    Q_lyap = 10*eye(6);                   % Lyapunov矩阵
    Gamma = 1*diag([1,1,1,1,1,1]);    % Kx自适应率
    Lambda = 1*diag([1,1,1]);         % Ka自适应率
    poles = [-2;-3;-2;-3;-2;-3];    % 配置极点
    Ka_max = 10;              % Ka限幅
    Kx_max = 10;            % Kx限幅
    epsilon_Ka = 1;               % 保护阈值
    epsilon_Kx = 1;               % 保护阈值
    % 2.2.4 -------------------- DRC参数 --------------------
    r0 = 50;    % 微分器跟踪速率
    h = 5;      % 滤波倍数
    c1_x = 1;   c2_x = 1;           % x通道反步增益
    c1_y = 1;   c2_y = 1;           % y通道反步增益
    c1_z = 1;   c2_z = 1;           % z通道反步增益
    L = -diag([10,10,10,10,10,10]); % 扰动观测器增益

% 2.3 ============================== 系统设置 ==============================
I = [0 1;0 0];  E = [0;1];      % 双积分形式
    % 2.3.1 轨迹选择
    % 0:    定点
    % 1:    8字形
    % 2:    螺旋形
    % 3：   方形
    traj_select = 1;
    % 2.3.2 
%% 3 初始化与预分配
% 3.1 ============================== 基本数据 ==============================
tspan = 0:dt:T;     % 仿真时域
N = length(tspan);  % 仿真步数
fprintf('[INFO] 程序初始化完成\n');

% 3.2 ============================== 系统信息 ==============================
fprintf('[INFO] 控制模式: %s (mode=%d)\n', mode_names{mode+1}, mode);
fprintf('[INFO] 仿真参数: dt=%.1e, T=%d, Ts=%.2f, horizon=%d\n', dt, T, Ts, horizen);
Ac = blkdiag(I,I,I);
Bc = blkdiag(E,E,E);
% dynamic_fun = @(x,u,t) Ac*x + Bc*(u...
%               .*max([1-0.00*t;1-0.00*t;1-0.00*t],0.1))... % 模拟扰动
%               +disturbance(t);  % 扰动
% dynamic_fun = @(x,u,t) Ac*x + Bc*u;   % 理想动态
% dynamic_fun = @(x,u,t) Ac*x + 0.5*Bc*u;   % 模拟动力衰减
dynamic_fun = @(x,u,t) Ac*x + Bc*u + disturbance(t);   % 模拟动力衰减

% 3.2 ============================== 轨迹参数 ==============================
Gen = TrajGen(0:dt:T+horizen*Ts);   % 多给定horizon*Ts
switch traj_select
    case 0  % 定点
        [x,y,z,vx,vy,vz,ax,ay,az] = Gen.FixedPoint(1,1,2);
        fprintf('[INFO] 轨迹类型: 定点\n');
    case 1  % 8字形
        [x,y,z,vx,vy,vz,ax,ay,az] = Gen.EightTraj(0,0,2,2,1,15,0);  % xc,yc,zc,rx,ry,T,dir
        fprintf('[INFO] 轨迹类型: 8字形\n');
    case 2  % 螺旋形
        [x,y,z,vx,vy,vz,ax,ay,az] = Gen.SpiralTraj(0,0,0.5,0.05,15,1,0); % xc,yc,z0,v,T,r,dir
        fprintf('[INFO] 轨迹类型: 螺旋形\n');
    case 3  % 方形
        [x,y,z,vx,vy,vz,ax,ay,az] = Gen.SquareTraj(0,0,2,2,30,0);
        fprintf('[INFO] 轨迹类型: 方形\n');
end

% 3.3 ============================== 工具初始化 ==============================
cnst_fun_Ka = @(Ka) (Ka(1)^2+Ka(2)^2+Ka(3)^2+Ka(4)^2+Ka(5)^2+Ka(6)^2+Ka(7)^2+Ka(8)^2+Ka(9)^2 ...
                     - Ka_max^2)/(2*Ka_max*epsilon_Ka + epsilon_Ka^2);
cnst_fun_Kx = @(Kx) (Kx(1)^2+Kx(2)^2+Kx(3)^2+Kx(4)^2+Kx(5)^2+Kx(6)^2+Kx(7)^2+Kx(8)^2+Kx(9)^2 + ...
                     Kx(10)^2+Kx(11)^2+Kx(12)^2+Kx(13)^2+Kx(14)^2+Kx(15)^2++Kx(16)^2+Kx(17)^2+Kx(18)^2 ...
                     - Kx_max^2)/(2*Kx_max*epsilon_Kx + epsilon_Kx^2);
H = place(Ac,Bc,poles);
P = lyap((Ac-Bc*H)', Q_lyap);
LAW = MRAC_Law(Ts,Gamma,Lambda,Bc,P,cnst_fun_Ka,cnst_fun_Kx);
DRCx = DRC_2nd(Ts,c1_x,c2_x,1);
DRCy = DRC_2nd(Ts,c1_y,c2_y,1);
DRCz = DRC_2nd(Ts,c1_z,c2_z,1);
DOB = LinearDOB(Ts,Ac,Bc,L,[],[]);
TD_x = TrckDiff(Ts,r0,h*Ts);
TD_y = TrckDiff(Ts,r0,h*Ts);
TD_z = TrckDiff(Ts,r0,h*Ts);
% 3.4 ============================== 变量初始化 ==============================
xi = zeros(6,N);        % 实际状态量
xi(:,1) = x_ini;
xi_mea = zeros(6,N);    % 测量状态量
xi_e = zeros(6,N);      % 误差状态量
xi_r = [x;vx;y;vy;z;vz];% 参考状态量

xi_m = xi(:,:);     % 参考模型（用于MRAC）
xi_n = xi(:,:);     % 标称模型（用于DRC）

u = zeros(3,N);         % 实际控制量（最终执行）
u_r = [ax;ay;az];       % 参考控制量
u_D = u;                % 抗扰作用
u_M = u;                % 自适应作用

Ka = zeros(3,3,N);    % 自适应参数
Kx = zeros(3,6,N);    % 自适应参数
Ka(:,:,1) = eye(3);

d = zeros(6,N);         % 扰动实际值
d_hat = zeros(6,N);     % 扰动估计值
% 3.5 ============================== MPC初始化 ==============================
prob = MPCSim_init(Ac,Bc,Q,R,horizen,Ts,xmax,xmin,umax,umin);
Ad = prob.Ad;
Bd = prob.Bd;
fprintf('[INFO] 系统初始化完成, 开始仿真计算 (N=%d 步)...\n', N);

%% 4 仿真计算
for k = 1:N-1
    if mod(tspan(k),Ts) == 0
        %% 4.0 获取信息
        xi_mea(:,k) = xi(:,k) + 0.00*randn(6,1);    % 加噪
        xi_e(:,k) = xi_r(:,k) - xi_mea(:,k);        % 转换为误差调节
        %% 4.1 控制算法
            % -------------------- 4.1.0 PID对比 --------------------
            if mode == 0
                TD_x.update(xi_e(1,k));
                TD_y.update(xi_e(3,k));
                TD_z.update(xi_e(5,k));
                u_pid = Kp.*xi_e([1,3,5],k) + ...
                        Ki.*sum(xi_e([1,3,5],1:Ts/dt:k),2) + ...
                        Kd.*[TD_x.output;TD_y.output;TD_z.output];
            end

            % -------------------- 4.1.1 标称MPC --------------------
            if mode ~= 0
                % 注意：必须取预测时域长度的参考值送入MPC
                ur_part = u_r(:,k:Ts/dt:k+Ts/dt*(horizen-1));      % 取预测时域长度，用于约束
                xr_part = xi_r(:,k+Ts/dt:Ts/dt:k+Ts/dt*horizen);   % 取预测时域长度，用于约束
                if mode == 1
                    u_mpc = u_r(:,k) - MPCSim_solv(prob,xi_e(:,k),xr_part,ur_part);   % 其他情况不应该直接采用\xi_e作为反馈！
                elseif mode == 2 || mode == 4
                    u_mpc = u_r(:,k) - MPCSim_solv(prob,xi_r(:,k)-xi_n(:,k),xr_part,ur_part);   % 抗扰时采用\xi_n！
                elseif mode == 3
                    u_mpc = u_r(:,k) - MPCSim_solv(prob,xi_r(:,k)-xi_m(:,k),xr_part,ur_part);   % 自适应时采用\xi_m！
                end
            end

            % -------------------- 4.1.2 DR-MPC ----------------------
            if mode == 2 || mode == 4
                xi_n(:,k+1) = Ad*xi_n(:,k) + Bd*u_mpc;  % 更新标称模型（处理非匹配）
                % 抗扰计算：
                DRCx.update(xi_mea(1:2,k),xi_n(1:2,k),d_hat(1:2,k));
                DRCy.update(xi_mea(3:4,k),xi_n(3:4,k),d_hat(3:4,k));
                DRCz.update(xi_mea(5:6,k),xi_n(5:6,k),d_hat(5:6,k));
                % 实施抗扰补偿：
                comp = [DRCx.compensation;DRCy.compensation;DRCz.compensation];
                u_D(:,k) = u_mpc + comp;
            elseif mode == 0
                u_D(:,k) = u_pid;
            else
                u_D(:,k) = u_mpc;
            end
            % 扰动观测：
            DOB.update(u_D(:,k),xi_mea(:,k))
            d_hat(:,k+1) = DOB.d_hat;

            % -------------------- 4.1.4 MRAPC -------------------
            if mode == 3 || mode == 4
                xi_m(:,k+1) = Ad*xi_m(:,k) + Bd*u_D(:,k);
                e_mrac = xi_m(:,k) - xi_mea(:,k);
                LAW.update(u_D(:,k),xi_mea(:,k),xi_m(:,k));
                Ka(:,:,k) = LAW.Ka;
                Kx(:,:,k) = LAW.Kx;
                u_M(:,k) = Ka(:,:,k)*(Kx(:,:,k)*xi_mea(:,k) + u_D(:,k))+H*e_mrac;
            else
                u_M(:,k) = u_D(:,k);
            end

        %% 4.2 实际执行
        % u(:,k) = min(max(u_M(:,k),umin),umax);    % 物理限幅
        u(:,k) = u_M(:,k);
    else
        %% 4.3 零阶保持：
        xi_mea(:,k) = xi_mea(:,k-1);
        xi_e(:,k) = xi_e(:,k-1);
        u(:,k) = u(:,k-1);

        d_hat(:,k+1) = d_hat(:,k);
        xi_n(:,k+1) = xi_n(:,k);
        u_D(:,k) = u_D(:,k-1);

        xi_m(:,k+1) = xi_m(:,k);
        Ka(:,:,k) = Ka(:,:,k-1);
        Kx(:,:,k) = Kx(:,:,k-1);
        u_M(:,k) = u_M(:,k-1);
    end

    %% 4.4 更新
    xi(:,k+1) = update_rk4(dynamic_fun,xi(:,k),u(:,k),tspan(k),dt);
    d(:,k) = disturbance(tspan(k));

    % 进度显示
    if mod(k, floor(N/10)) == 0
        progress = round(k / N * 100);
        bar_str = repmat('=', 1, floor(progress/5));
        fprintf('\r[INFO] 仿真进度: [%s%-*s] %3d%%', bar_str, 20-length(bar_str), '', progress);
    end
end
fprintf('[INFO] 仿真计算完成!\n\n');

%% 预览
figure(1)
subplot(3,2,1)
plot(tspan,xi(1,:),tspan,xi_r(1,1:N))
title("x")
hold on
subplot(3,2,2)
plot(tspan,xi(2,:),tspan,xi_r(2,1:N))
title("v_x")
hold on
subplot(3,2,3)
plot(tspan,xi(3,:),tspan,xi_r(3,1:N))
title("y")
hold on
subplot(3,2,4)
plot(tspan,xi(4,:),tspan,xi_r(4,1:N))
title("v_y")
hold on
subplot(3,2,5)
plot(tspan,xi(5,:),tspan,xi_r(5,1:N))
title("z")
hold on
subplot(3,2,6)
plot(tspan,xi(6,:),tspan,xi_r(6,1:N))
title("v_z")
hold on


figure(2)
plot3(xi(1,:),xi(3,:),xi(5,:))
hold on
plot3(xi_r(1,:),xi_r(3,:),xi_r(5,:))
title("Trajectory")
grid on

figure(3)
subplot(3,1,1)
plot(tspan,u(1,:),tspan,u_r(1,1:N))
title("u_x")
hold on
subplot(3,1,2)
plot(tspan,u(2,:),tspan,u_r(2,1:N))
title("u_y")
subplot(3,1,3)
plot(tspan,u(3,:),tspan,u_r(3,1:N))
title("u_z")

figure(4)
subplot(3,2,1)
plot(tspan,d_hat(1,:),tspan,d(1,1:N))
title("d_1")
hold on
subplot(3,2,2)
plot(tspan,d_hat(2,:),tspan,d(2,1:N))
title("d_2")
subplot(3,2,3)
plot(tspan,d_hat(3,:),tspan,d(3,1:N))
title("d_3")
subplot(3,2,4)
plot(tspan,d_hat(4,:),tspan,d(4,1:N))
title("d_4")
subplot(3,2,5)
plot(tspan,d_hat(5,:),tspan,d(5,1:N))
title("d_5")
subplot(3,2,6)
plot(tspan,d_hat(6,:),tspan,d(6,1:N))
title("d_6")
%% 指标
fprintf('[INFO] ========== 跟踪误差指标 ==========\n');
kss = floor(N/2);
num = N-kss+1;
ex2 = (xi(1,kss:N)-xi_r(1,kss:N)).^2;
ey2 = (xi(3,kss:N)-xi_r(3,kss:N)).^2;
ez2 = (xi(5,kss:N)-xi_r(5,kss:N)).^2;
fprintf('[INFO]   位置 RMSE: X=%.4f(m)  Y=%.4f(m)  Z=%.4f(m)\n', sqrt(sum(ex2)/num), sqrt(sum(ey2)/num), sqrt(sum(ez2)/num));
fprintf('[INFO]   总位置 RMSE: %.4f (m)\n', sqrt(sum(ex2+ey2+ez2)/num));

evx2 = (xi(2,kss:N)-xi_r(2,kss:N)).^2;
evy2 = (xi(4,kss:N)-xi_r(4,kss:N)).^2;
evz2 = (xi(6,kss:N)-xi_r(6,kss:N)).^2;
fprintf('[INFO]   速度 RMSE: Vx=%.4f(m/s)  Vy=%.4f(m/s)  Vz=%.4f(m/s)\n', sqrt(sum(evx2)/num), sqrt(sum(evy2)/num), sqrt(sum(evz2)/num));
fprintf('[INFO]   总速度 RMSE: %.4f (m/s)\n', sqrt(sum(evx2+evy2+evz2)/num));
fprintf('[INFO] ==================================\n');
%% 保存数据
YN = input("\n[INFO] 是否保存实验数据? (Y/N): ", "s");
if YN == 'Y' || YN == 'y'
    xi_r = xi_r(:,1:N);     % 只保留有效部分
    u_r = u_r(:,1:N);     % 只保留有效部分
    str = ['Data\chap4\data_',num2str(mode),'.mat'];
    expr_name = expr_names{traj_select+1};
    save(str,"tspan","xi_r","u_r","xi","u","d","d_hat","expr_name");
    fprintf('[INFO] 数据已保存至: %s\n', str);
    
    % 询问是否保存 Kx 和 Ka
    YN_K = input("[INFO] 是否保存自适应参数 Kx 与 Ka? (Y/N): ", "s");
    if YN_K == 'Y' || YN_K == 'y'
        str_K = 'Data\chap4\adapt.mat';
        Ka_save = Ka(:,:,1:N);
        Kx_save = Kx(:,:,1:N);
        save(str_K, "Ka_save", "Kx_save");
        fprintf('[INFO] 自适应参数已保存至: %s\n', str_K);
    else
        fprintf('[INFO] 自适应参数未保存\n');
    end
else
    fprintf('[INFO] 已取消保存\n');
end
%% 保存Kx与Ka
