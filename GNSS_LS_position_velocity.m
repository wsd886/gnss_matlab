function [est_r_ea_e,est_v_ea_e,est_clock] = GNSS_LS_position_velocity(...
    GNSS_measurements,no_GNSS_meas,predicted_r_ea_e,predicted_v_ea_e)
%GNSS_LS_position_velocity - 使用非加权迭代最小二乘法计算位置、速度、钟差和钟漂。
% 对位置/钟差和速度/钟漂分别进行独立计算。
%(Calculates position, velocity, clock offset, and clock drift using 
% unweighted iterated least squares. Separate calculations are implemented 
% for position and clock offset and for velocity and clock drift)
%
% 配套书籍 "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," 第二版 软件
% (Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.)
%
% This function created 11/4/2012 by Paul Groves
%
% 输入参数 (Inputs):
%   GNSS_measurements     GNSS 测量数据矩阵:
%     第1列              伪距测量值 (米)
%     第2列              伪距率测量值 (米/秒)
%     第3-5列            卫星 ECEF 位置 (米)
%     第6-8列            卫星 ECEF 速度 (米/秒)
%   no_GNSS_meas          提供测量值的卫星数量
%   predicted_r_ea_e      先验预测的用户 ECEF 位置 (米)，用作迭代初值
%   predicted_v_ea_e      先验预测的用户 ECEF 速度 (米/秒)，用作迭代初值
%
% 输出参数 (Outputs):
%   est_r_ea_e            估计的用户 ECEF 位置 (米)
%   est_v_ea_e            估计的用户 ECEF 速度 (米/秒)
%   est_clock             估计的接收机钟差 (米) 和钟漂 (米/秒) [2x1列向量]
 
% 版权所有 2012, Paul Groves
% 许可证: BSD; 详情请见 license.txt

% 常量定义 (其中一些未来可以改为输入参数)
c = 299792458; % 光速 (米/秒)
omega_ie = 7.292115E-5;  % 地球自转角速度 (弧度/秒)

% 开始 (Begins)

% --- 第一部分：位置和钟差估计 ---
% POSITION AND CLOCK OFFSET

% 初始化迭代的状态预测值
x_pred(1:3,1) = predicted_r_ea_e; % 使用输入的先验位置作为位置初值
x_pred(4,1) = 0; % 钟差初值设为 0 (相对于先验位置的修正)
test_convergence = 1; % 初始化收敛测试变量，确保至少执行一次循环

% 迭代计算，直到收敛 (状态变化小于阈值 0.0001)
while test_convergence > 0.0001
    
    % 遍历所有卫星的测量值
    for j = 1:no_GNSS_meas
        % 1. 预测近似的星地距离 (用于计算信号传播时间)
        delta_r = GNSS_measurements(j,3:5)' - x_pred(1:3); % 卫星位置 - 用户近似位置
        approx_range = sqrt(delta_r' * delta_r); % 计算距离标量
        
        % 2. 计算信号传播期间的地球自转修正矩阵 C_e_I (式 8.36)
        %    用于将信号发射时刻的卫星位置旋转到信号接收时刻
        C_e_I = [1, omega_ie * approx_range / c, 0;...
                 -omega_ie * approx_range / c, 1, 0;...
                 0, 0, 1];
                 
        % 3. 预测伪距测量值 (式 9.143)
        %    首先，将卫星位置旋转到接收时刻的 ECEF 坐标
        delta_r = C_e_I * GNSS_measurements(j,3:5)' - x_pred(1:3);
        %    然后，计算几何距离
        range = sqrt(delta_r' * delta_r);
        %    最后，加上当前估计的钟差 (x_pred(4)) 得到伪距预测值
        pred_meas(j,1) = range + x_pred(4);
        
        % 4. 计算设计矩阵 H 的第 j 行 (式 9.144)
        %    首先计算从用户指向卫星的视线单位向量 (Line of Sight, LOS) 的相反数
        H_matrix (j,1:3) = - delta_r' / range; % 对位置的偏导数
        H_matrix (j,4) = 1; % 对钟差(c*dtu)的偏导数
        
    end % 结束遍历卫星 for j
        
    % 5. 计算非加权最小二乘解 (式 9.35 / 9.141)
    %    计算状态修正量 delta_x = inv(H'*H) * H' * (测量值 - 预测值)
    %    注意：这里 H_matrix 只取了实际使用的 no_GNSS_meas 行
    x_est = x_pred + inv(H_matrix(1:no_GNSS_meas,:)' *...
        H_matrix(1:no_GNSS_meas,:)) * H_matrix(1:no_GNSS_meas,:)' *...
        (GNSS_measurements(1:no_GNSS_meas,1) - pred_meas(1:no_GNSS_meas));
        
    % 6. 测试收敛性
    %    计算本次迭代估计值与上次迭代估计值(预测值)之间的差的范数
    test_convergence = sqrt((x_est - x_pred)' * (x_est - x_pred));
    
    % 7. 更新预测值，用于下一次迭代
    x_pred = x_est;
    
end % 结束迭代 while

% 将最终迭代结果赋值给输出变量
est_r_ea_e(1:3,1) = x_est(1:3); % 估计的用户 ECEF 位置
est_clock(1) = x_est(4); % 估计的接收机钟差 (米)


% --- 第二部分：速度和钟漂估计 ---
% VELOCITY AND CLOCK DRIFT

% 计算地球自转角速度的反对称矩阵 (用于叉乘计算)
Omega_ie = Skew_symmetric([0,0,omega_ie]); % 需要 Skew_symmetric 函数
       
% 初始化迭代的状态预测值
x_pred(1:3,1) = predicted_v_ea_e; % 使用输入的先验速度作为速度初值
x_pred(4,1) = 0; % 钟漂初值设为 0 (相对于先验速度的修正)
test_convergence = 1; % 初始化收敛测试变量

% 迭代计算，直到收敛
while test_convergence > 0.0001
    
    % 遍历所有卫星的测量值
    for j = 1:no_GNSS_meas
        % 1. 计算近似星地距离 (基于第一部分最终估计的位置 est_r_ea_e)
        delta_r = GNSS_measurements(j,3:5)' - est_r_ea_e;
        approx_range = sqrt(delta_r' * delta_r);
        
        % 2. 计算信号传播期间的地球自转修正矩阵 C_e_I
        C_e_I = [1, omega_ie * approx_range / c, 0;...
                 -omega_ie * approx_range / c, 1, 0;...
                 0, 0, 1];
                 
        % 3. 计算几何距离 (基于最终估计位置)
        delta_r = C_e_I * GNSS_measurements(j,3:5)' - est_r_ea_e;
        range = sqrt(delta_r' * delta_r);
        
        % 4. 计算视线单位向量 u_as_e (式 8.41)
        u_as_e = delta_r / range;
        
        % 5. 预测伪距率测量值 (式 9.143 的速率形式)
        %    需要考虑卫星速度、用户速度以及两者因地球自转产生的速度分量
        %    卫星速度项: C_e_I * (卫星ECEF速度 + Omega_ie * 卫星ECEF位置)
        %    用户速度项: (用户ECEF速度预测值 + Omega_ie * 用户ECEF估计位置)
        range_rate = u_as_e' * (C_e_I * (GNSS_measurements(j,6:8)' +...
            Omega_ie * GNSS_measurements(j,3:5)') - (x_pred(1:3) +...
            Omega_ie * est_r_ea_e));        
        %    加上当前估计的钟漂 (x_pred(4)) 得到伪距率预测值
        pred_meas(j,1) = range_rate + x_pred(4);
        
        % 6. 计算设计矩阵 H 的第 j 行 (式 9.144 的速率形式)
        %    伪距率对速度的偏导数是 LOS 向量的相反数
        H_matrix (j,1:3) = - u_as_e'; % 对速度的偏导数
        H_matrix (j,4) = 1; % 对钟漂(c*ddtu)的偏导数
        
    end % 结束遍历卫星 for j
    
    % 7. 计算非加权最小二乘解
    %    计算状态修正量 delta_x = inv(H'*H) * H' * (测量值 - 预测值)
    %    注意：测量值是伪距率 GNSS_measurements(..., 2)
    x_est = x_pred + inv(H_matrix(1:no_GNSS_meas,:)' *...
        H_matrix(1:no_GNSS_meas,:)) * H_matrix(1:no_GNSS_meas,:)' *...
        (GNSS_measurements(1:no_GNSS_meas,2) - pred_meas(1:no_GNSS_meas));
        
    % 8. 测试收敛性
    test_convergence = sqrt((x_est - x_pred)' * (x_est - x_pred));
    
    % 9. 更新预测值，用于下一次迭代
    x_pred = x_est;
    
end % 结束迭代 while

% 将最终迭代结果赋值给输出变量
est_v_ea_e(1:3,1) = x_est(1:3); % 估计的用户 ECEF 速度
est_clock(2) = x_est(4); % 估计的接收机钟漂 (米/秒)

% 结束 (Ends)
