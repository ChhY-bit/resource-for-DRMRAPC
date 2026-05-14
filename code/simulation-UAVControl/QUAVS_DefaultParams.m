function params = QUAVS_DefaultParams()
    params = zeros(15,1);
    params(1) = 3.13e-5;    %   C_A         -   1.升力系数
    params(2) = 7.50e-7;    %   C_R         -   2.反力矩系数
    params(3) = 0.232;      %   L           -   3.旋翼臂长
    params(4) = 0.52;       %   m           -   4.无人机质量
    params(5) = 6.228e-3;   %   Jxx         -   5.x主轴惯量
    params(6) = 6.225e-3;   %   Jyy         -   6.y主轴惯量
    params(7) = 1.121e-2;   %   Jzz         -   7.z主轴惯量
    % ---------- 扰动因素 ----------
    params(8) = 0.00e-5;    %   Jr          -   8.翼桨惯量
    params(9) = 0.0;    %   mu_x        -   9.x方向平移阻力系数
    params(10) = 0.0;   %   mu_y        -   10.y方向平移阻力系数
    params(11) = 0.0;   %   mu_z        -   11.z方向平移阻力系数
    params(12) = 0.0;    %   lambda_x    -   12.x方向旋转阻力系数
    params(13) = 0.0;    %   lambda_y    -   13.y方向旋转阻力系数
    params(14) = 0.0;    %   lambda_z    -   14.z方向旋转阻力系数
    % ---------- -------- ----------
    params(15) = 9.81;    %   g           -   15.重力加速度
end
