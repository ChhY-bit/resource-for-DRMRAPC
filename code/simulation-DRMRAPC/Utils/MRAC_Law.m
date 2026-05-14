classdef MRAC_Law < handle
    %UNTITLED2 此处显示有关此类的摘要
    %   此处显示详细说明

    properties
        Ka      % 总增益       -   (p x p)
        Kx      % 状态反馈增益 -   (p x n)
        Lambda  % Ka自适应率   -   (p x p)
        Gamma   % Kx自适应率   -   (n x n)
        P       % Lyapunov方程矩阵  -   (n x n)
        Ts      % 采样时间
        Bm      % 参考输入矩阵
        dimensions % 维度信息

        ifproj  % 是否开启投影算子
        cnst_fun_Ka
        cnst_fun_Kx
        cnst_grd_Ka
        cnst_grd_Kx

        ifnorm % 是否开启归一化
    end

    methods
        function obj = MRAC_Law(Ts,Gamma,Lambda,Bm,P,cst_Ka,cst_Kx)
            %MRAC 构造此类的实例
            %   此处显示详细说明
            if nargin == 5 || isempty(cst_Ka) || isempty(cst_Kx)
                obj.ifproj = "off";
            else
                obj.ifproj = "on";
            end
            if nargin <= 7
                ifnorm = "on";
            end

            obj.dimensions.n = size(Gamma,1);
            obj.dimensions.p = size(Lambda,1);
            obj.Ts = Ts;
            obj.Bm = Bm;
            obj.Gamma = Gamma;
            obj.Lambda = Lambda;
            obj.P = P;
            obj.Ka = eye(obj.dimensions.p);
            obj.Kx = zeros(obj.dimensions.p,obj.dimensions.n);
            obj.ifnorm = ifnorm;
            if obj.ifproj == "on"
                obj.cnst_fun_Ka = cst_Ka;   % Ka的约束函数
                obj.cnst_fun_Kx = cst_Kx;   % Kx的约束函数
                Ka_num = obj.dimensions.p * obj.dimensions.p; %#ok<*NASGU>
                Kx_num = obj.dimensions.p * obj.dimensions.n;
                syms x_sym [Ka_num,1]
                grad = jacobian(obj.cnst_fun_Ka(x_sym),x_sym);
                obj.cnst_grd_Ka = matlabFunction(grad, 'Vars', {x_sym});    % Ka约束函数的梯度
                syms x_sym [Kx_num,1]
                grad = jacobian(obj.cnst_fun_Kx(x_sym),x_sym);
                obj.cnst_grd_Kx = matlabFunction(grad, 'Vars', {x_sym});    % Kx约束函数的梯度
            end
        end

        function update(obj,u,x,xm)
            %UPDATE 此处显示有关此方法的摘要
            % u     -   当前输入
            % x     -   当前实际状态
            % xm    -   当前参考状态
            e = xm - x;
            % 范数归一化：
            if obj.ifnorm == "on"
                Delta_Ka = obj.Ka * obj.Bm' * obj.P * e * (u + obj.Kx*x)' *...
                          obj.Ka' * obj.Lambda * obj.Ka...
                          ./(norm(u + obj.Kx * x)^2+1);
                Delta_Kx = obj.Bm' * obj.P * e * x' * obj.Gamma...
                          ./(norm(x)^2+1);
            end
            % 投影算子：
            if obj.ifnorm == "on"
                Ka_list = obj.Ka + obj.Ts*Delta_Ka;
                Ka_list = Ka_list(:);
                Delta_Ka = Proj(Ka_list,Delta_Ka(:),obj.cnst_fun_Ka,obj.cnst_grd_Ka);
                Delta_Ka = reshape(Delta_Ka,obj.dimensions.p,obj.dimensions.p);

                Kx_list = obj.Kx + obj.Ts*Delta_Kx;
                Kx_list = Kx_list(:);
                Delta_Kx = Proj(Kx_list,Delta_Kx(:),obj.cnst_fun_Kx,obj.cnst_grd_Kx);
                Delta_Kx = reshape(Delta_Kx,obj.dimensions.p,obj.dimensions.n);
            end
            % 参数更新：
            obj.Ka = obj.Ka + obj.Ts*Delta_Ka;
            obj.Kx = obj.Kx + obj.Ts*Delta_Kx; 
        end
    end
end