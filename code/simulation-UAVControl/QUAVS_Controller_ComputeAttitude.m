function [phi_d,theta_d]=QUAVS_Controller_ComputeAttitude(ax,ay,az,psi_d)
% QUAVS_Controller_ComputeAttitude
% 输入期望加速度ax,ay,az和期望偏航角psi_d
% 计算得出期望的滚转角phi_d、俯仰角theta_d
    alpha = [ax,ay,az+9.81]./norm([ax,ay,az+9.81]);
    theta_d = atan2(sin(psi_d)*alpha(2)+cos(psi_d)*alpha(1),alpha(3));
    phi_d = atan2(sin(psi_d)*alpha(1)-cos(psi_d)*alpha(2),alpha(3)/cos(theta_d));
end