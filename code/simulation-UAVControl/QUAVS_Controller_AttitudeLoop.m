function w_d=QUAVS_Controller_AttitudeLoop(eta_d,eta_d_dot,eta)
% QUAVS_Controller_AttitudeLoop
% 输入期望姿态eta_d及其导数eta_d_dot、实际姿态eta
% 输出机体坐标系下应达到的角速度w_d
    err = eta_d - eta;
    % 控制参数
    K_eta = diag([1,1,1]); % 应设计为正定的
    
    % 控制律
    [~,W_inv]=W(eta(1),eta(2),eta(3));
    w_d = W_inv*(eta_d_dot+K_eta*err);
end