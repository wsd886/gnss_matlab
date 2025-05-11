%GNSS_Demo_1
%脚本: 独立的GNSS最小二乘解算演示
%   使用 Profile_1 (包含两次90度转弯的60秒人工车辆运动轨迹)
%
% 配套书籍 "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," 第二版 软件
% (本书籍 "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems" 第二版配套软件)
%
% 创建于 11/4/12 by Paul Groves
% 版权所有 2012, Paul Groves
% 许可证: BSD; 详情请见 license.txt

%% 初始化与常量定义

% 清理工作环境
clear;      % 清除工作空间变量
close all;  % 关闭所有图形窗口
clc;        % 清除命令行窗口

% 定义常量
deg_to_rad = 0.01745329252; % 角度转弧度的转换系数
rad_to_deg = 1/deg_to_rad; % 弧度转角度的转换系数
micro_g_to_meters_per_second_squared = 9.80665E-6; % 微g转米/秒^2的转换系数 (此脚本中可能未使用)

%% 配置参数 (CONFIGURATION)

% 输入/输出文件名设置
input_profile_name = 'Profile_1.csv'; % 输入的真实运动轨迹文件名 (CSV格式)
output_profile_name = 'GNSS_Demo_1_Profile.csv'; % 输出的计算轨迹文件名 (CSV格式)
output_errors_name = 'GNSS_Demo_1_Errors.csv';   % 输出的定位误差文件名 (CSV格式)

% GNSS模拟与处理相关配置 (存放在 GNSS_config 结构体中)
GNSS_config.epoch_interval = 1; % GNSS解算的历元时间间隔 (秒)

% 初始估计位置 (ECEF坐标系, 单位: 米)
% 这是最小二乘解算开始时使用的用户位置初始猜测值
GNSS_config.init_est_r_ea_e = [0;0;0];

% 模拟的卫星星座参数
GNSS_config.no_sat = 30; % 星座中的卫星总数
GNSS_config.r_os = 2.656175E7; % 卫星轨道半径 (米)
GNSS_config.inclination = 55; % 卫星轨道倾角 (度)
GNSS_config.const_delta_lambda = 0; % 星座的经度偏移 (度)
GNSS_config.const_delta_t = 0; % 星座的时间偏移 (秒)

% 接收机参数
GNSS_config.mask_angle = 10; % 卫星高度截止角 (度)，低于此角度的卫星不用于定位解算

% 模拟的GNSS误差模型标准差 (Standard Deviation, SD)
% 注意: 这些通常指经过模型改正后的 *残余* 误差
GNSS_config.SIS_err_SD = 1; % 卫星信号在空间传播的误差 (米)
GNSS_config.zenith_iono_err_SD = 2; % 天顶方向电离层延迟残差 (米)
GNSS_config.zenith_trop_err_SD = 0.2; % 天顶方向对流层延迟残差 (米)
GNSS_config.code_track_err_SD = 1; % 接收机码环跟踪噪声 (米)，可扩展以模拟多路径效应
GNSS_config.rate_track_err_SD = 0.02; % 接收机载波环(多普勒)跟踪噪声 (米/秒)，可扩展以模拟多路径效应

% 模拟的接收机钟差模型参数
GNSS_config.rx_clock_offset = 10000; % 接收机在 t=0 时刻的初始钟差 (等效为距离，米)
GNSS_config.rx_clock_drift = 100;   % 接收机在 t=0 时刻的初始钟漂 (等效为距离变化率，米/秒)

% 随机数生成器种子设置 (适用于较新版本的MATLAB)
% 设置种子是为了保证每次运行脚本时，模拟产生的随机误差序列相同，便于结果复现。
% 更改种子值会得到不同的随机误差序列。
% 'twister' 指定使用 Mersenne Twister 算法 (mt19937ar)。
rng(1, 'twister');
% 旧版MATLAB语法 (已注释掉):
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));

%% 脚本主体执行 (Begins)

% 从CSV文件读取输入的真实运动轨迹数据
% Read_profile 是一个自定义函数，用于解析CSV文件
% in_profile: 包含轨迹数据的矩阵或结构体
% no_epochs: 轨迹数据的总历元数
% ok: 标志位，表示文件是否成功读取
[in_profile,no_epochs,ok] = Read_profile(input_profile_name);

% 如果文件读取失败，则显示错误信息并终止脚本
if ~ok
    disp(['无法读取输入文件: ', input_profile_name]);
    return;
end %if

% 执行GNSS最小二乘解算
% GNSS_Least_Squares 是核心处理函数，它接收真实轨迹、历元数和配置参数
% 该函数内部会:
% 1. 根据真实轨迹和配置参数，模拟生成带有误差的GNSS测量值(伪距等)
% 2. 使用最小二乘法处理这些模拟的测量值，计算每个历元的接收机位置和可能的钟差
% 3. 计算定位结果与真实轨迹之间的误差
% 输出:
% out_profile: 计算得到的运动轨迹
% out_errors: 定位误差统计
% out_clock: 可能包含计算出的钟差信息 (根据函数实现)
[out_profile,out_errors,out_clock] = GNSS_Least_Squares(in_profile,...
    no_epochs,GNSS_config);

% 绘图显示结果 (注意: 在Octave中可能不完全兼容)
close all; % 关闭所有已打开的图形窗口
Plot_profile(in_profile); % 调用 Plot_profile 函数绘制输入的真实轨迹
Plot_errors(out_errors);  % 调用 Plot_errors 函数绘制计算出的定位误差

% 将输出结果写入CSV文件
Write_profile(output_profile_name,out_profile); % 调用 Write_profile 函数保存计算出的轨迹
Write_errors(output_errors_name,out_errors);   % 调用 Write_errors 函数保存定位误差

% 脚本结束 (Ends)
disp('GNSS_Demo_1 脚本执行完毕。');
disp(['计算得到的轨迹已保存至: ', output_profile_name]);
disp(['定位误差已保存至: ', output_errors_name]);
