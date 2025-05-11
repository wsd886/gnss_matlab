function [out_profile,out_errors,out_clock] = GNSS_Least_Squares(in_profile,...
    no_epochs,GNSS_config)
%GNSS_Least_Squares - 使用最小二乘定位算法模拟独立的GNSS定位
% (Simulates stand-alone GNSS using a least-squares positioning algorithm)
%
% 配套书籍 "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," 第二版 软件
% (Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.)
%
% This function created 11/4/2012 by Paul Groves
%
% 输入参数 (Inputs):
%   in_profile   真实的运动轨迹数组 (行: 时间历元, 列: 见下方说明)
%   no_epochs    轨迹数据的总历元数
%   GNSS_config  包含GNSS配置参数的结构体
%     .epoch_interval     GNSS历元之间的时间间隔 (秒)
%     .init_est_r_ea_e    初始估计位置 (米; ECEF坐标系)
%     .no_sat             星座中的卫星数量
%     .r_os               卫星轨道半径 (米)
%     .inclination        卫星轨道倾角 (度)
%     .const_delta_lambda 星座的经度偏移 (度)
%     .const_delta_t      星座的时间偏移 (秒)
%     .mask_angle         卫星高度截止角 (度)
%     .SIS_err_SD         信号空间误差标准差 (米)
%     .zenith_iono_err_SD 天顶方向电离层延迟误差标准差 (米)
%     .zenith_trop_err_SD 天顶方向对流层延迟误差标准差 (米)
%     .code_track_err_SD  码跟踪误差标准差 (米)
%     .rate_track_err_SD  测距率(多普勒)跟踪误差标准差 (米/秒)
%     .rx_clock_offset    接收机在t=0时的钟差 (米)
%     .rx_clock_drift     接收机在t=0时的钟漂 (米/秒)
%
% 输出参数 (Outputs):
%   out_profile   导航解算结果，格式同运动轨迹数组
%   out_errors    导航解算误差数组 (格式见下方说明)
%   out_clock     接收机钟差估计数组 (格式见下方说明)
%
% 运动轨迹数组格式 (Format of motion profiles):
%  第1列: 时间 (秒)
%  第2列: 纬度 (弧度)
%  第3列: 经度 (弧度)
%  第4列: 高度 (米)
%  第5列: 北向速度 (米/秒)
%  第6列: 东向速度 (米/秒)
%  第7列: 天向(下)速度 (米/秒) - 注意: 通常指地向速度
%  第8列: 横滚角 (载体坐标系相对于NED坐标系) (弧度)
%  第9列: 俯仰角 (载体坐标系相对于NED坐标系) (弧度)
%  第10列: 偏航角 (载体坐标系相对于NED坐标系) (弧度)
%
% 误差数组格式 (Format of error array):
%  第1列: 时间 (秒)
%  第2列: 北向位置误差 (米)
%  第3列: 东向位置误差 (米)
%  第4列: 天向(下)位置误差 (米) - 注意: 通常指地向误差
%  第5列: 北向速度误差 (米/秒)
%  第6列: 东向速度误差 (米/秒)
%  第7列: 天向(下)速度误差 (米/秒) - 注意: 通常指地向误差
%  第8列: 未使用 (北向姿态误差 (弧度))
%  第9列: 未使用 (东向姿态误差 (弧度))
%  第10列: 未使用 (天向(下)姿态误差 = 航向误差 (弧度))
%
% 接收机钟差数组格式 (Format of receiver clock array):
%  第1列: 时间 (秒)
%  第2列: 估计的钟差 (米)
%  第3列: 估计的钟漂 (米/秒)

% 版权所有 2012, Paul Groves
% 许可证: BSD; 详情请见 license.txt

% 开始 (Begins)

% --- 初始化第一个历元的真实导航状态 ---
old_time = in_profile(1,1); % 获取第一个历元的时间
true_L_b = in_profile(1,2); % 真实纬度 (弧度)
true_lambda_b = in_profile(1,3); % 真实经度 (弧度)
true_h_b = in_profile(1,4); % 真实高度 (米)
true_v_eb_n = in_profile(1,5:7)'; % 真实速度 (北东地, m/s)，转置为列向量
true_eul_nb = in_profile(1,8:10)'; % 真实姿态欧拉角 (rad)，转置为列向量
true_C_b_n = Euler_to_CTM(true_eul_nb)'; % 计算真实的姿态矩阵 C_b^n (从载体到导航)
% 将真实的 NED 坐标系下的位置和速度转换为 ECEF 坐标系
[true_r_eb_e,true_v_eb_e] =...
    pv_NED_to_ECEF(true_L_b,true_lambda_b,true_h_b,true_v_eb_n);

% 初始化 GNSS 处理时间戳和历元计数器
time_last_GNSS = old_time;
GNSS_epoch = 1;

% --- 计算第一个历元的卫星位置和速度 ---
% 调用 Satellite_positions_and_velocities 函数计算该时刻所有卫星的 ECEF 位置和速度
[sat_r_es_e,sat_v_es_e] = Satellite_positions_and_velocities(old_time,...
    GNSS_config);

% --- 初始化 GNSS 偏差 (大气延迟等) ---
% 注意: 此函数假设偏差在仿真期间是恒定的，且基于初始的卫星仰角。
% 因此，这个仿真函数不适用于超过约30分钟的长时间仿真。
% 调用 Initialize_GNSS_biases 函数初始化各种偏差值
GNSS_biases = Initialize_GNSS_biases(sat_r_es_e,true_r_eb_e,true_L_b,...
    true_lambda_b,GNSS_config);

% --- 生成第一个历元的模拟 GNSS 测量值 ---
% 调用 Generate_GNSS_measurements 函数模拟生成伪距、伪距率等测量值
% 输入包括：时间、卫星位置速度、用户真实位置速度、大气偏差、配置参数
[GNSS_measurements,no_GNSS_meas] = Generate_GNSS_measurements(old_time,...
    sat_r_es_e,sat_v_es_e,true_r_eb_e,true_L_b,true_lambda_b,true_v_eb_e,...
    GNSS_biases,GNSS_config);

% --- 执行第一个历元的 GNSS 最小二乘定位解算 ---
% 调用 GNSS_LS_position_velocity 函数进行最小二乘解算
% 输入包括：模拟的测量值、可见卫星数、初始估计位置、初始估计速度(这里设为0)
% 输出：估计的 ECEF 位置、速度和钟差/钟漂
[est_r_eb_e,est_v_eb_e,est_clock] = GNSS_LS_position_velocity(...
    GNSS_measurements,no_GNSS_meas,GNSS_config.init_est_r_ea_e,[0;0;0]);

% 设定估计的姿态等于真实姿态 (因为此函数不解算姿态，姿态误差设为零)
est_C_b_n = true_C_b_n;

% 将估计的 ECEF 位置和速度转换回 NED 坐标系
[est_L_b,est_lambda_b,est_h_b,est_v_eb_n] =...
    pv_ECEF_to_NED(est_r_eb_e,est_v_eb_e);

% --- 记录第一个历元的输出结果 ---
% 记录解算出的轨迹
out_profile(1,1) = old_time;        % 时间
out_profile(1,2) = est_L_b;         % 纬度
out_profile(1,3) = est_lambda_b;    % 经度
out_profile(1,4) = est_h_b;         % 高度
out_profile(1,5:7) = est_v_eb_n';   % 速度 (北东地)
out_profile(1,8:10) = CTM_to_Euler(est_C_b_n')'; % 姿态 (从姿态矩阵转回欧拉角)

% 计算并记录定位和测速误差
% 调用 Calculate_errors_NED 函数计算估计值与真值之间的误差
[delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
    est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n,true_L_b,...
    true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
out_errors(1,1) = old_time;         % 时间
out_errors(1,2:4) = delta_r_eb_n';  % 位置误差 (北东地)
out_errors(1,5:7) = delta_v_eb_n';  % 速度误差 (北东地)
out_errors(1,8:10) = [0;0;0];       % 姿态误差设为零

% 记录估计的钟差和钟漂
out_clock(1,1) = old_time;          % 时间
out_clock(1,2:3) = est_clock(1:2);  % 钟差(米), 钟漂(米/秒)

% --- 初始化进度条 ---
dots = '....................'; % 总共20个点
bars = '||||||||||||||||||||'; % 总共20个竖线
rewind = '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'; % 20个退格符，用于覆盖之前的进度条
fprintf(strcat('Processing: ',dots)); % 打印初始进度条
progress_mark = 0; % 当前进度标记数
progress_epoch = 0; % 上次更新进度条的历元

% --- 主循环：处理从第2个到最后一个历元 ---
for epoch = 2:no_epochs

    % 更新进度条显示 (大约每处理5%的数据更新一次)
    if (epoch - progress_epoch) > (no_epochs/20)
        progress_mark = progress_mark + 1;
        progress_epoch = epoch;
        fprintf(strcat(rewind,bars(1:progress_mark),... % 打印退格符、已完成的竖线
            dots(1:(20 - progress_mark)))); % 打印剩余的点
    end % if epoch

   % 获取当前历元的时间
   time = in_profile(epoch,1);

    % 判断是否到达下一个 GNSS 解算时刻
    if (time - time_last_GNSS) >= GNSS_config.epoch_interval
        GNSS_epoch = GNSS_epoch + 1; % GNSS 输出结果的历元计数器加 1
        time_last_GNSS = time; % 更新上一次 GNSS 解算的时间

        % --- 获取当前历元的真实运动状态 ---
        true_L_b = in_profile(epoch,2);
        true_lambda_b = in_profile(epoch,3);
        true_h_b = in_profile(epoch,4);
        true_v_eb_n = in_profile(epoch,5:7)';
        true_eul_nb = in_profile(epoch,8:10)';
        true_C_b_n = Euler_to_CTM(true_eul_nb)';
        [true_r_eb_e,true_v_eb_e] =...
            pv_NED_to_ECEF(true_L_b,true_lambda_b,true_h_b,true_v_eb_n);

        % --- 计算当前时刻的卫星位置和速度 ---
        [sat_r_es_e,sat_v_es_e] = Satellite_positions_and_velocities(time,...
            GNSS_config);

        % --- 生成当前时刻的模拟 GNSS 测量值 ---
        % 注意: GNSS_biases 在这里被假定为常量，没有随时间更新
        [GNSS_measurements,no_GNSS_meas] = Generate_GNSS_measurements(...
            time,sat_r_es_e,sat_v_es_e,true_r_eb_e,true_L_b,true_lambda_b,...
            true_v_eb_e,GNSS_biases,GNSS_config);

        % --- 执行当前时刻的 GNSS 最小二乘定位解算 ---
        % 使用上一个历元的估计位置和速度作为本次解算的初始值
        [est_r_eb_e,est_v_eb_e,est_clock] = GNSS_LS_position_velocity(...
            GNSS_measurements,no_GNSS_meas,est_r_eb_e,est_v_eb_e);

        % 将估计的 ECEF 位置速度转换回 NED 坐标系
        [est_L_b,est_lambda_b,est_h_b,est_v_eb_n] =...
            pv_ECEF_to_NED(est_r_eb_e,est_v_eb_e);
        % 姿态误差仍设为零
        est_C_b_n = true_C_b_n;

        % --- 记录当前历元的输出结果 ---
        % 记录解算出的轨迹
        out_profile(GNSS_epoch,1) = time;
        out_profile(GNSS_epoch,2) = est_L_b;
        out_profile(GNSS_epoch,3) = est_lambda_b;
        out_profile(GNSS_epoch,4) = est_h_b;
        out_profile(GNSS_epoch,5:7) = est_v_eb_n';
        out_profile(GNSS_epoch,8:10) = CTM_to_Euler(est_C_b_n')';

        % 计算并记录误差
        [delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
            est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n,true_L_b,...
            true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
        out_errors(GNSS_epoch,1) = time;
        out_errors(GNSS_epoch,2:4) = delta_r_eb_n';
        out_errors(GNSS_epoch,5:7) = delta_v_eb_n';
        out_errors(GNSS_epoch,8:10) = [0;0;0]; % 姿态误差设为零

        % 记录估计的钟差和钟漂
        out_clock(GNSS_epoch,1) = time;
        out_clock(GNSS_epoch,2:3) = est_clock(1:2);

        % 更新旧时间戳 (虽然在此循环结构下这行不是必需的)
        old_time = time;
    end % if (time - time_last_GNSS) >= GNSS_config.epoch_interval
end % 结束主循环 for epoch = 2:no_epochs

% 完成进度条显示
fprintf(strcat(rewind,bars,'\n')); % 打印最终的满格进度条并换行

% 结束 (Ends)
