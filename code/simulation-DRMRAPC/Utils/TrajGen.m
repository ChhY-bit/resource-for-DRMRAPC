classdef TrajGen < handle
    properties
        tspan
    end
    
    methods
        function obj = TrajGen(tspan)
            obj.tspan = tspan;
        end

        function [x,y,z,dx,dy,dz,ddx,ddy,ddz] = SquareTraj(obj,xc,yc,zc,L,T,dir)
            % 正方形轨迹：从(xc+L,yc+L,zc)开始
            % xc    -   中心x坐标(m)
            % yc    -   中心y坐标(m)
            % zc    -   中心z坐标(m)
            % L     -   每边边长(m)
            % T     -   完成一个轨迹周期的时间(s)
            % dir   -   运动方向    0:逆时针   1:顺时针
            v = L/(T/4);    % 速度
            fun = @(t) (t>=0 & t<T/4).* (-v*(t-T/8)) + ...
                       (t>=T/4 & t<T/2).* (-L/2) + ...
                       (t>=T/2 & t<3*T/4).* (v*(t-5*T/8)) + ...
                       (t>=3*T/4).* L/2;
            fun_d = @(t) (t>=0 & t<T/4).* (-v) + ...
                       (t>=T/4 & t<T/2).* 0 + ...
                       (t>=T/2 & t<3*T/4).* v + ...
                       (t>=3*T/4).* 0;
            fun_dd = @(t) 0.*t;
            if dir == 1
                x = xc + fun(mod(obj.tspan-T/4,T));
                y = yc + fun(mod(obj.tspan,T));
                dx = fun_d(mod(obj.tspan-T/4,T));
                dy = fun_d(mod(obj.tspan,T));
                ddx = fun_dd(mod(obj.tspan-T/4,T));
                ddy = fun_dd(mod(obj.tspan,T));
            else
                x = xc + fun(mod(obj.tspan,T));
                y = yc + fun(mod(obj.tspan-T/4,T));
                dx = xc + fun_d(mod(obj.tspan,T));
                dy = yc + fun_d(mod(obj.tspan-T/4,T));
                ddx = xc + fun_dd(mod(obj.tspan,T));
                ddy = yc + fun_dd(mod(obj.tspan-T/4,T));
            end
            z = zc*ones(1,length(obj.tspan));
            dz = zeros(1,length(obj.tspan));
            ddz = zeros(1,length(obj.tspan));
        end

        function [x,y,z,dx,dy,dz,ddx,ddy,ddz] = SpiralTraj(obj,xc,yc,z0,v,T,r,dir)
            % 生成螺旋轨迹（从第一象限开始）
            % xc    -   中心x坐标(m)
            % yc    -   中心y坐标(m)
            % z0    -   起始z坐标(m)
            % v     -   上升速度(m/s)
            % T     -   周期时间(s)
            % r     -   半径(m)
            % dir   -   运动方向    0:逆时针   1:顺时针
            w = 2*pi/T;
            x_fun = @(t) xc+r*cos(w*t);
            y_fun = @(t) yc+r*sin(w*t);
            dx_fun = @(t) -r*w*sin(w*t);
            dy_fun = @(t) r*w*cos(w*t);
            ddx_fun = @(t) -r*w^2*cos(w*t);
            ddy_fun = @(t) -r*w^2*sin(w*t);
            if dir == 1
                x = y_fun(obj.tspan);
                dx = dy_fun(obj.tspan);
                ddx = ddy_fun(obj.tspan);
                y = x_fun(obj.tspan);
                dy = dx_fun(obj.tspan);
                ddy = ddx_fun(obj.tspan);
            else
                x = x_fun(obj.tspan);
                dx = dx_fun(obj.tspan);
                ddx = ddx_fun(obj.tspan);
                y = y_fun(obj.tspan);
                dy = dy_fun(obj.tspan);
                ddy = ddy_fun(obj.tspan);
            end
            z = z0+ v*obj.tspan;
            dz = v.*ones(1,length(obj.tspan));
            ddz = zeros(1,length(obj.tspan));
        end

        function [x,y,z,dx,dy,dz,ddx,ddy,ddz] = EightTraj(obj,xc,yc,zc,rx,ry,T,dir)
            % 8字形轨迹
            % xc    -   中心x坐标
            % yc    -   中心y坐标
            % zc    -   中心z坐标
            % rx    -   方向x上的半径
            % ry    -   方向y上的半径
            % T     -   完成一个轨迹周期的时间
            % dir   -   形状朝向
            w = 2*pi/T;
            if dir == 0     % 朝x方向
                x = xc + rx*sin(w*obj.tspan);
                dx = w*rx*cos(w*obj.tspan);
                ddx = -w^2*rx*sin(w*obj.tspan);

                y = yc + ry*sin(2*w*obj.tspan);
                dy = 2*w*ry*cos(2*w*obj.tspan);
                ddy = -4*w^2*ry*sin(2*w*obj.tspan);
            elseif dir == 1 % 朝y方向
                x = xc + rx*sin(2*w*obj.tspan);
                dx = 2*w*rx*cos(2*w*obj.tspan);
                ddx = -4*w^2*rx*cos(2*w*obj.tspan);

                y = yc + ry*sin(w*obj.tspan);
                dy = w*ry*cos(w*obj.tspan);
                ddy = -w^2*ry*sin(w*obj.tspan);
            end
            z = zc*ones(1,length(obj.tspan));
            dz = zeros(1,length(obj.tspan));
            ddz = zeros(1,length(obj.tspan));
        end

        function [x,y,z,dx,dy,dz,ddx,ddy,ddz] = FixedPoint(obj,xc,yc,zc)
            % 定点轨迹
            % xc    -   定点x坐标
            % yc    -   定点y坐标
            % zc    -   定点z坐标
            x = xc*ones(1,length(obj.tspan));
            dx = zeros(1,length(obj.tspan));
            ddx = zeros(1,length(obj.tspan));
            y = yc*ones(1,length(obj.tspan));
            dy = zeros(1,length(obj.tspan));
            ddy = zeros(1,length(obj.tspan));
            z = zc*ones(1,length(obj.tspan));
            dz = zeros(1,length(obj.tspan));
            ddz = zeros(1,length(obj.tspan));
        end
    end
end