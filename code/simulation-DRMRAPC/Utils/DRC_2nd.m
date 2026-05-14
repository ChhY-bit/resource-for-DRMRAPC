classdef DRC_2nd < handle
    %DRC 二阶扰动补偿器
    %   适用于双积分系统
    %   针对匹配扰动，直接在输出通道实施补偿
    %   针对非匹配扰动，使用反步法实施补偿

    properties
        compensation        % 补偿控制量
        TD_e                % 对跟踪误差的跟踪微分器
        TD_d1               % 对非匹配扰动的跟踪微分器
        c1                  % 增益1
        c2                  % 增益2
        Ts                  % 采样时间
        if_nonmatch         % 是否加入非匹配补偿
    end

    methods
        function obj = DRC_2nd(Ts,c1,c2,if_nonmatch)
            %DRC_2nd 构造此类的实例
            %   Ts  -   采样时间
            %   c1  -   增益1
            %   c2  -   增益2
            if nargin == 3
                if_nonmatch = false;   % 默认不开启
            end
            obj.Ts = Ts;
            obj.c1 = c1;
            obj.c2 = c2;
            obj.TD_d1 = TrckDiff(Ts,50,2*Ts);
            obj.TD_e = TrckDiff(Ts,50,2*Ts);
            obj.compensation = 0;
            obj.if_nonmatch = if_nonmatch;
        end

        function set_TD_d1(obj,r,h0)
            obj.TD_d1.r = r;
            obj.TD_d1.h0 = h0;
        end

        function set_TD_e(obj,r,h0)
            obj.TD_e.r = r;
            obj.TD_e.h0 = h0;
        end

        function update(obj,x,x_n,d)
        %UPDATE 更新补偿量
        %   x   -   实际状态量 (2 x 1)
        %   x_n -   标称状态量 (2 x 1)
        %   d   -   扰动值 (2 x 1)
        e = x-x_n;  ep=e(1);    ev=e(2);

        obj.TD_d1.update(d(1));
        d1_diff = obj.TD_d1.output;

        obj.TD_e.update(ep);
        ep_diff = obj.TD_e.output;
        
        w = -(obj.c1*obj.c2+1)*ep...
            - obj.c2*ev - obj.c1*ep_diff- obj.c2*d(1) - d1_diff;

        obj.compensation = obj.if_nonmatch*w - d(2);
        end
    end
end