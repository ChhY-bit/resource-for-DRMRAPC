classdef LinearDOB < handle
    %LINEARDOB 线性扰动观测器
    %   LinearDOB 利用MATLAB类提供了线性扰动观测器的实现
    %   初始化：LinearDOB(Ts,A,B,L,Bd,x_ini)
    %   设置增益：LinearDOB.set_gain(L)
    %   更新：LinearDOB.update(u,x)
    %   获取扰动观测值：LinearDOB.d_hat
    %
    %   作者：<C.Yang> Github: <https://github.com/ChhY-bit>
    %   联系邮箱：<yangchenhan.work@gmail.com>
    %            <ych_0872@126.com>
    %   版本：v1.0
    %   更新日期：2026-04-11
    
    properties
        Ts          % 采样时间
        A           % 标称系统矩阵    -   (n x n)
        B           % 标称输入矩阵    -   (n x p)
        Bd          % 扰动输入矩阵    -   (n x q)
        d_hat       % 扰动观测值      -   (q x 1)
        L           % 观测增益        -   (q x n)
        z           % 中间变量        -   (q x 1)
        dimensions  % 维数信息
        x_last      % 上一时刻的状态  -   (n x 1)
        u_last      % 上一时刻的输入  -   (p x 1)
    end
    
    methods
        function obj = LinearDOB(Ts,A,B,L,Bd,x_ini)
            %LINEARDOB 线性扰动观测器的初始化
            %   输入参数：
            %   Ts   - 采样时间
            %   A    - 标称系统矩阵 (n x n)
            %   B    - 标称输入矩阵 (n x p)
            %   L    - 观测增益 (q x n)
            %   Bd   - 扰动输入矩阵 (n x q)，默认为单位矩阵
            %   x_ini - 初始状态 (n x 1)，默认为零向量
            %
            %   输出参数：
            %   obj  - 线性扰动观测器对象
            %
            %   注意：
            %       Bd默认为单位矩阵
            %       x_ini默认为零向量
            %       不具备全局存储功能，需外界记录状态
            if isempty(Bd) || isempty(x_ini)
                Bd = eye(size(A,1));        % 缺省值为单位矩阵
                x_ini = zeros(size(A,1),1);   % 缺省值为零初状态
            end
            % 初始化：
            obj.Ts = Ts;
            obj.dimensions.n = size(A,1);
            obj.dimensions.p = size(B,2);
            obj.dimensions.q = size(Bd,1);
            obj.A = A;
            obj.B = B;
            obj.L = L;
            obj.Bd = Bd;
            obj.x_last = x_ini;
            obj.u_last = zeros(obj.dimensions.p,1);
            obj.z = obj.L*obj.x_last;
            obj.d_hat = zeros(obj.dimensions.q,1);
        end
        
        function set_gain(obj,L)
            %SET_GAIN 设置观测器增益
            %   输入参数：
            %   L    - 观测增益 (q x n)
            %
            %   注意：
            %       若维数不匹配将输出警告
            [m,n] = size(L);
            if m == obj.dimensions.q && n == obj.dimensions.n
                obj.L = L;
            else
                warning("Wrong dimensions of DOB Gain!")
            end
        end
        
        function update(obj,u,x)
            %UPDATE 更新扰动观测
            %   输入参数：
            %   u     - 当前系统输入 (p x 1)
            %   x     - 当前系统状态 (n x 1)
            %
            %   注意：
            %       应当在每一个采样周期内调用此方法
            %       调用后通过 obj.d_hat 获取扰动观测值
            temp = obj.L*(obj.A*obj.x_last + obj.B*obj.u_last + obj.d_hat);
            obj.z = temp*obj.Ts + obj.z;
            obj.d_hat = obj.z - obj.L*x;
            obj.x_last = x;
            obj.u_last = u;
        end
        
    end
end