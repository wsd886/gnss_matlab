%INS_GNSS_Demo_1
%脚本：紧耦合 INS/GNSS 演示：
%   Profile_1（60 秒的人工车辆运动轨迹，包含两个 90 度转弯）
%   战术级惯性测量单元（IMU）
%
% 本脚本配合《GNSS、惯性及多传感器集成导航原理（第二版）》一书使用
%
% 作者：Paul Groves，创建时间：2012年12月4日

% 版权所有 © 2012, Paul Groves
% 授权协议：BSD；详见 license.txt 文件

% 常量定义
deg_to_rad = 0.01745329252;  % 角度转弧度
rad_to_deg = 1/deg_to_rad;  % 弧度转角度
micro_g_to_meters_per_second_squared = 9.80665E-6;  % 微重力加速度单位换算

% 配置参数
% 输入的真实运动轨迹文件名
input_profile_name = 'Profile_1.csv';
% 输出轨迹与误差文件名
output_profile_name = 'INS_GNSS_Demo_1_Profile.csv';
output_errors_name = 'INS_GNSS_Demo_1_Errors.csv';

% 初始姿态误差（角度转换为弧度；在北、东、地方向）
initialization_errors.delta_eul_nb_n = [-0.05;0.04;1]*deg_to_rad; % 单位：弧度

% 加速度计偏置（微g，转换为 m/s^2；机体系）
IMU_errors.b_a = [900;-1300;800] * micro_g_to_meters_per_second_squared;
% 陀螺仪偏置（度/小时，转换为弧度/秒；机体系）
IMU_errors.b_g = [-9;13;-8] * deg_to_rad / 3600;
% 加速度计比例因子与交叉耦合误差（ppm，转换为无量纲；机体系）
IMU_errors.M_a = [500, -300, 200;...
                 -150, -600, 250;...
                 -250,  100, 450] * 1E-6;
% 陀螺仪比例因子与交叉耦合误差（ppm，转换为无量纲；机体系）
IMU_errors.M_g = [400, -300,  250;...
                    0, -300, -150;...
                    0,    0, -350] * 1E-6;
% 陀螺仪 g 敏感偏置（度/小时/g，转换为 弧度/秒/m/s²；机体系）
IMU_errors.G_g = [0.9, -1.1, -0.6;...
                 -0.5,  1.9, -1.6;...
                  0.3,  1.1, -1.3] * deg_to_rad / (3600 * 9.80665);
% 加速度计噪声根功率谱密度（微g/√Hz，转换为 m/s^1.5）
IMU_errors.accel_noise_root_PSD = 100 * micro_g_to_meters_per_second_squared;
% 陀螺仪噪声根功率谱密度（度/√小时，转换为 rad/s^0.5）
IMU_errors.gyro_noise_root_PSD = 0.01 * deg_to_rad / 60;
% 加速度计量化级别（m/s^2）
IMU_errors.accel_quant_level = 1E-2;
% 陀螺仪量化级别（rad/s）
IMU_errors.gyro_quant_level = 2E-4;

% GNSS历元间隔（单位：秒）
GNSS_config.epoch_interval = 0.5;

% 初始估计位置（ECEF坐标，单位：米）
GNSS_config.init_est_r_ea_e = [0;0;0];

% 卫星数目
GNSS_config.no_sat = 30;
% 卫星轨道半径（单位：米）
GNSS_config.r_os = 2.656175E7;
% 卫星轨道倾角（单位：度）
GNSS_config.inclination = 55;
% 星座经度偏移（单位：度）
GNSS_config.const_delta_lambda = 0;
% 星座时间偏移（单位：秒）
GNSS_config.const_delta_t = 0;

% GNSS观测仰角遮罩（单位：度）
GNSS_config.mask_angle = 10;
% 空间信号误差标准差（单位：米）*当使用修正后的残差
GNSS_config.SIS_err_SD = 1;
% 电离层垂直误差标准差（单位：米）*当使用修正后的残差
GNSS_config.zenith_iono_err_SD = 2;
% 对流层垂直误差标准差（单位：米）*当使用修正后的残差
GNSS_config.zenith_trop_err_SD = 0.2;
% 伪距码跟踪误差标准差（单位：米）*可扩展用于多路径误差
GNSS_config.code_track_err_SD = 1;
% 伪距速率跟踪误差标准差（单位：米/秒）*可扩展用于多路径误差
GNSS_config.rate_track_err_SD = 0.02;
% 初始接收机时钟偏差（单位：米）
GNSS_config.rx_clock_offset = 10000;
% 初始接收机时钟漂移（单位：米/秒）
GNSS_config.rx_clock_drift = 100;

% 初始姿态不确定性（每轴角度，单位转为弧度）
TC_KF_config.init_att_unc = degtorad(1);
% 初始速度不确定性（单位：m/s）
TC_KF_config.init_vel_unc = 0.1;
% 初始位置不确定性（单位：m）
TC_KF_config.init_pos_unc = 10;
% 初始加速度计偏置不确定性（单位：微g 转 m/s^2）
TC_KF_config.init_b_a_unc = 1000 * micro_g_to_meters_per_second_squared;
% 初始陀螺仪偏置不确定性（单位：度/小时 转 rad/s）
TC_KF_config.init_b_g_unc = 10 * deg_to_rad / 3600;
% 初始时钟偏移不确定性（单位：m）
TC_KF_config.init_clock_offset_unc = 10;
% 初始时钟漂移不确定性（单位：m/s）
TC_KF_config.init_clock_drift_unc = 0.1;

% 陀螺仪噪声功率谱密度（单位：deg²/h 转 rad²/s）
TC_KF_config.gyro_noise_PSD = (0.02 * deg_to_rad / 60)^2;
% 加速度计噪声功率谱密度（单位：微g²/Hz 转 m²/s³）
TC_KF_config.accel_noise_PSD = (200 * micro_g_to_meters_per_second_squared)^2;
% 加速度计偏置随机游走PSD（单位：m²/s⁵）
TC_KF_config.accel_bias_PSD = 1.0E-7;
% 陀螺仪偏置随机游走PSD（单位：rad²/s³）
TC_KF_config.gyro_bias_PSD = 2.0E-12;
% 接收机时钟频率漂移PSD（单位：m²/s³）
TC_KF_config.clock_freq_PSD = 1;
% 接收机时钟相位漂移PSD（单位：m²/s）
TC_KF_config.clock_phase_PSD = 1;

% 伪距测量噪声标准差（单位：米）
TC_KF_config.pseudo_range_SD = 2.5;
% 伪距率测量噪声标准差（单位：米/秒）
TC_KF_config.range_rate_SD = 0.1;

% 设置随机数生成器种子，确保每次运行结果一致（在 Octave 中可能无效）
rng(1, 'twister');

% 开始仿真流程

% 从CSV文件中读取真实运动轨迹数据
[in_profile,no_epochs,ok] = Read_profile(input_profile_name);

% 如果读取失败，则退出脚本
if ~ok
    return;
end %if

% 执行紧耦合 ECEF 惯性导航和 GNSS 集成导航仿真
[out_profile,out_errors,out_IMU_bias_est,out_clock,out_KF_SD] =...
    Tightly_coupled_INS_GNSS(in_profile,no_epochs,initialization_errors...
    ,IMU_errors,GNSS_config,TC_KF_config);

% 绘制真实运动轨迹与误差图（Octave 可能不支持）
close all;
Plot_profile(in_profile);
Plot_errors(out_errors);

% 将输出轨迹与误差写入CSV文件
Write_profile(output_profile_name,out_profile);
Write_errors(output_errors_name,out_errors);

% 结束脚本
