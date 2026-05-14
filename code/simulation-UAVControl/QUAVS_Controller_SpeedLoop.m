function [Tx,Ty,Tz]=QUAVS_Controller_SpeedLoop(w_d,w,J)
% QUAVS_Controller_SpeedLoop
% 输入期望角速度w_d及实际角速度w、机体惯量矩阵J
% 输出控制力矩Tx,Ty,Tz
    err = w_d - w;
    % 控制参数
    K_w = diag([3,3,3]); % 应设计为正定的
    
    % 控制律
    T = J*K_w*err;
    Tx = T(1);  Ty = T(2);  Tz = T(3);
end