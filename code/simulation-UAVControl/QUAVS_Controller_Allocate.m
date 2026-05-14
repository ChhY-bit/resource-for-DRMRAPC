function [Omg1,Omg2,Omg3,Omg4]=QUAVS_Controller_Allocate(F,Tx,Ty,Tz)
% QUAVS_Controller_Allocate
% 输入控制量F,Tx,Ty,Tz
% 输出分配后的电机转速

    % 参数
    C_A = 3.13e-5;
    C_R = 7.50e-7;
    L = 0.232;
    B = diag([C_A,L*C_A/sqrt(2),L*C_A/sqrt(2),C_R])*...
             [1 1 1 1;1 1 -1 -1;-1 1 1 -1;-1 1 -1 1];
    % 分配
    OMG = sqrt(B\[F;Tx;Ty;Tz]);
    Omg1 = OMG(1); Omg2 = OMG(2); Omg3 = OMG(3); Omg4 = OMG(4);
end