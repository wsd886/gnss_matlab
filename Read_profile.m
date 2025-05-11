function [in_profile,no_epochs,ok]= Read_profile(filename)
%Read_profile - 从指定的 .csv 文件格式中读取运动轨迹数据
% (Read_profile - inputs a motion profile in the following .csv format)
% CSV文件列格式说明:
% 第1列: 时间 (秒)
% 第2列: 纬度 (度)
% 第3列: 经度 (度)
% 第4列: 高度 (米)
% 第5列: 北向速度 (米/秒)
% 第6列: 东向速度 (米/秒)
% 第7列: 地向速度 (米/秒)
% 第8列: 横滚角 (载体坐标系 w.r.t NED) (度)
% 第9列: 俯仰角 (载体坐标系 w.r.t NED) (度)
% 第10列: 偏航角 (载体坐标系 w.r.t NED) (度)
%
% 配套书籍 "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," 第二版 软件
% (Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.)
%
% This function created 31/3/2012 by Paul Groves
%
% 输入参数 (Inputs):
%   filename     要读取的文件名 (字符串)
%
% 输出参数 (Outputs):
%   in_profile   从文件中读取的数据数组
%   no_epochs    文件中数据的历元(行)数量
%   ok           逻辑值，指示文件是否具有预期的列数 (true/false)

% 版权所有 2012, Paul Groves
% 许可证: BSD; 详情请见 license.txt

% 开始 (Begins)

% 参数定义
deg_to_rad = 0.01745329252; % 角度转弧度的转换系数

% 从 .csv 文件读取轨迹数据
in_profile = csvread(filename); % 使用 csvread 函数读取文件内容

% 确定文件的大小 (行数和列数)
[no_epochs,no_columns] = size(in_profile); % 获取读取数据的行数和列数

% 检查列数是否正确 (如果不正确则返回)
if no_columns~=10 % 判断列数是否不等于 10
    disp('输入文件的列数错误') % 显示错误提示信息
    ok = false; % 设置输出标志 ok 为 false
else
    ok = true; % 如果列数等于 10，设置输出标志 ok 为 true
    % 将角度单位从度转换为弧度
    % Convert degrees to radians
    in_profile(:,2:3) = deg_to_rad * in_profile(:,2:3); % 转换纬度(第2列)和经度(第3列)
    in_profile(:,8:10) = deg_to_rad * in_profile(:,8:10); % 转换横滚角(第8列)、俯仰角(第9列)和偏航角(第10列)
end % 结束 if 判断

% 结束 (Ends)
