classdef QUAV < handle
    %QUAV 四旋翼无人机模型
    %   创建一个确定的四旋翼无人机模型
    %   参数说明：
    %   C_A         -   1.气动阻力常数
    %   C_R         -   2.反力矩系数
    %   L           -   3.旋翼臂长
    %   m           -   4.无人机质量
    %   Jxx         -   5.x主轴惯量
    %   Jyy         -   6.y主轴惯量
    %   Jzz         -   7.z主轴惯量
    %   Jr          -   8.翼桨惯量
    %   mu_x        -   9.x方向平移阻力系数
    %   mu_y        -   10.y方向平移阻力系数
    %   mu_z        -   11.z方向平移阻力系数
    %   lambda_x    -   12.x方向旋转阻力系数
    %   lambda_y    -   13.y方向旋转阻力系数
    %   lambda_z    -   14.z方向旋转阻力系数
    %   g           -   15.重力加速度
    properties
        Parameters
        States
    end

    methods
        function obj = QUAV(params,states_ini)
            % QUAV 构造此类的实例
            %   创建无人机示例，参数初始化
            if nargin == 0
                params = QUAVS_DefaultParams();
                states_ini = zeros(12,1);
            elseif length(params(:)) ~= 15
                error('Wrong Number of Parameters!')
            end
            obj.Parameters.C_A = params(1);
            obj.Parameters.C_R = params(2);
            obj.Parameters.L = params(3);
            obj.Parameters.m = params(4);
            obj.Parameters.Jxx = params(5);
            obj.Parameters.Jyy = params(6);
            obj.Parameters.Jzz = params(7);
            obj.Parameters.Jr = params(8);
            obj.Parameters.mu_x = params(9);
            obj.Parameters.mu_y = params(10);
            obj.Parameters.mu_z = params(11);
            obj.Parameters.lambda_x = params(12);
            obj.Parameters.lambda_y = params(13);
            obj.Parameters.lambda_z = params(14);
            obj.Parameters.g = params(15);

            obj.States.dx = states_ini(1);
            obj.States.dy = states_ini(2);
            obj.States.dz = states_ini(3);
            obj.States.x = states_ini(4);
            obj.States.y = states_ini(5);
            obj.States.z = states_ini(6);
            obj.States.w_x = states_ini(7);
            obj.States.w_y = states_ini(8);
            obj.States.w_z = states_ini(9);
            obj.States.phi = states_ini(10);
            obj.States.theta = states_ini(11);
            obj.States.psi = states_ini(12);
        end

        function state_update(obj,Omega_1,Omega_2,Omega_3,Omega_4,t,dt)
            %state_input 依据电机转速指令输入更新状态
            %   参数说明
            %   Omega_1     -   1号电机转速
            %   Omega_2     -   2号电机转速
            %   Omega_3     -   3号电机转速
            %   Omega_4     -   4号电机转速
            %   dt          -   仿真时间步长
            it = obj.Parameters;
            % 分配矩阵：
            B = diag([it.C_A,...
                      it.L*it.C_A/sqrt(2),...
                      it.L*it.C_A/sqrt(2),...
                      it.C_R])*...
                      [1 1 1 1;1 1 -1 -1;-1 1 1 -1;-1 1 -1 1];
            % 计算动力：
            tau = B*[Omega_1^2;Omega_2^2;Omega_3^2;Omega_4^2];
            %% 动力学方程
            % 数据准备
            eta = [obj.States.phi,obj.States.theta,obj.States.psi]';
            p = [obj.States.x,obj.States.y,obj.States.z]';
            dp = [obj.States.dx,obj.States.dy,obj.States.dz]';
            w = [obj.States.w_x,obj.States.w_y,obj.States.w_z]';
            mu = -diag([it.mu_x,it.mu_y,it.mu_z]);
            J = diag([it.Jxx,it.Jyy,it.Jzz]);
            lambda = diag([it.lambda_x,it.lambda_y,it.lambda_z]);
            % 动态方程
            Omega_r = Omega_1-Omega_2+Omega_3-Omega_4;
            Lambda = @(w)[0,-it.Jr*Omega_r-it.Jzz*w(3),it.Jyy*w(2);...
                          it.Jr*Omega_r+it.Jzz*w(3),0,-it.Jxx*w(1);...
                          -it.Jyy*w(2),it.Jxx*w(1),0]-lambda;
            dynaics = @(dp,w,eta,tau)...
                [dp;... % \dot{p}
                 (mu*dp-[0;0;it.m*it.g]+R(eta(1),eta(2),eta(3))*[0;0;tau(1)])/it.m;... % \ddot{p}
                 J\(Lambda(w)*w+tau(2:4));... % \dot{\omega}
                 W(eta(1),eta(2),eta(3))*w]; % \dot{\eta}
            %% 状态更新
            x = [p;dp;w;eta];
            u = tau;
            fun = @(x,u,t)dynaics(x(4:6),x(7:9),x(10:12),u);
            new_states=update_rk4(fun,x,u,t,dt);
            obj.States.x = new_states(1);
            obj.States.y = new_states(2);
            obj.States.z = new_states(3);
            obj.States.dx = new_states(4);
            obj.States.dy = new_states(5);
            obj.States.dz = new_states(6);
            obj.States.w_x = new_states(7);
            obj.States.w_y = new_states(8);
            obj.States.w_z = new_states(9);
            obj.States.phi = new_states(10);
            obj.States.theta = new_states(11);
            obj.States.psi = new_states(12);
        end
    end
end