function F_d=QUAVS_Controller_ComputeForce(ax,ay,az,m)
% QUAVS_Controller_ComputeForc
% 设质量为m，输入期望加速度ax,ay,az
% 计算得出期望的升力大小F_d
    F_d = m*norm([ax,ay,az+9.81]);
end