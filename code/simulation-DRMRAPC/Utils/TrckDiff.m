classdef TrckDiff   < handle
    %TRCKDIFF 跟踪微分器
    %   TrckDiff 利用MATLAB类提供了跟踪微分器的实现
    %   初始化：TrckDiff(h,r,h0)
    %   更新：TrckDiff.update(input)
    %   输出：[diff,trck] = TrckDiff.output()
    %   设置状态：TrckDiff.set_states(x1,x2)
    %
    %   作者：<C.Yang> Github: <https://github.com/ChhY-bit>
    %   联系邮箱：<yangchenhan.work@gmail.com>
    %            <ych_0872@126.com>
    %   版本：v1.0
    %   更新日期：2026-4-04-07
    
    properties
        h   % 积分步长（通常与采样周期相同）
        r   % 跟踪速率
        h0  % 滤波因子
        x1  % 状态变量1（信号输入）
        x2  % 状态变量2（微分输出）
    end
    
    methods
        function obj = TrckDiff(h,r,h0)
            %TRCKDIFF 跟踪微分器的初始化
            %   输入参数：
            %   h   - 积分步长（即采样周期）
            %   r   - 跟踪速率
            %   h0  - 滤波因子
            %
            %   输出参数：
            %   obj - 跟踪微分器对象
            %
            %   注意：
            %       h0默认与h相同
            %       x1,x2初始默认为0，若需更改请使用set_states()方法
            if nargin == 2
                h0 = h;
            end
            obj.h = h;
            obj.r = r;
            obj.h0 = h0;
            obj.x1 = 0;
            obj.x2 = 0;
        end
        
        function set_states(obj,x1,x2)
            %SET_STATES 设置跟踪微分器状态值
            %   输入参数：
            %   x1  - 状态变量1（信号输入）
            %   x2  - 状态变量2（微分输出）
            obj.x1 = x1;
            obj.x2 = x2;
        end

        function set_params(obj,h,r,h0)
            %SET_PARAMS 设置跟踪微分器参数值
            %   输入参数：
            %   h   - 积分步长（即采样周期）
            %   r   - 跟踪速率
            %   h0  - 滤波因子
            obj.h = h;
            obj.r = r;
            obj.h0 = h0;
        end

        function update(obj,input)
            %UPDATE 跟踪微分器更新
            %   输入参数：
            %   input - 信号输入
            %
            % 注意：应当在每一次采样周期内、获取微分结果之前调用此方法
            err = obj.x1 - input;
            obj.x1 = obj.x1 + obj.h*obj.x2;
            obj.x2 = obj.x2 + obj.h*obj.fhan(err,obj.x2,obj.r,obj.h0);
        end

        function [diff,trck] = output(obj)
            %OUTPUT 跟踪微分器输出
            %   输出参数：
            %   diff - 微分结果
            %   trck - 跟踪信号
            diff = obj.x2;
            trck = obj.x1;
        end

        function res = fhan(obj,x1,x2,r,h0)
            %FHAN 辅助函数
            d = r * h0^2;
            a0 = h0 * x2;
            y = x1 + a0;
            a1 = sqrt(d*(d+8*abs(y)));
            a2 = a0 + (sign(y)*(a1-d))/2;
            a = (a0+y)*obj.fsg(y,d) + a2*(1-obj.fsg(y,d));
            res = -r*(a/d*obj.fsg(a,d)+sign(a)*(1-obj.fsg(a,d)));
        end

        function res = fsg(~,x,d)
            %FSG 辅助函数
            res = (sign(x+d)-sign(x-d))/2;
        end
    end
end

