function [cost, results] = geneva_objective(params, return_results)
% 槽轮机构的目标函数，用来评价一组参数好不好
%
% 输入:
%   params - 4个优化参数，分别是:
%       params(1): G点偏转角，范围2~10度
%       params(2): 输入行程，范围25~45度
%       params(3): 速度边界，范围-60~-5
%       params(4): 加速度边界，范围1~30
%   return_results - 要不要返回详细数据，默认不返回
%
% 输出:
%   cost - 目标函数值，越小说明这组参数越好
%   results - 运动曲线数据，画图用的

if nargin < 2
    return_results = false;
end

% 把参数拆出来，方便后面用
Angle_G_rel = params(1);
Input_Stroke = params(2);
v_E = params(3);
a_E = params(4);

% 机构的基本参数，这些是固定的
R = 55;           % 中心距，单位mm
N = 10;           % 槽轮有10个槽
r_pin = 3;        % 圆销半径
Omega_RPM = 60;   % 转速60rpm
Omega = Omega_RPM * 2*pi/60;

% 根据基本参数算一些几何量
d = 2 * R * sin(pi/N);
r_rot = d/2;
O1 = [-R, 0];
Angle_Pitch = 360 / N;
Phi_end = -Angle_G_rel;
Phi_start = 0;

% 解五次多项式的系数，这是运动规律的核心
C0 = Phi_start;
Delta_Phi = Phi_end - C0;
M = [1 1 1; 3 4 5; 6 12 20];
V = [Delta_Phi; v_E; a_E];
X = M \ V;
C3 = X(1); C4 = X(2); C5 = X(3);

% 生成查找表，后面仿真要用
% 点数设成5000，太多了会慢
num_gen = 5000;
Raw_Radius = zeros(1, 2*num_gen);
Raw_Angle = zeros(1, 2*num_gen);
Raw_Type = zeros(1, 2*num_gen);

% 先算停歇段(CD4段)的轨迹
Angle_D4_Abs = -Input_Stroke;
Angle_Extend = Angle_D4_Abs - 200;
theta_dwell = linspace(Angle_Extend, Angle_D4_Abs, num_gen);

for i = 1:num_gen
    p_abs = O1 + r_rot * [cosd(theta_dwell(i)), sind(theta_dwell(i))];
    Raw_Radius(i) = sqrt(p_abs(1)^2 + p_abs(2)^2);
    Raw_Angle(i) = atan2(p_abs(2), p_abs(1));
    Raw_Type(i) = 1;
end

% 再算运动段(D4G段)的轨迹
tau_gen = linspace(0, 1, num_gen);
for i = 1:num_gen
    t = tau_gen(i);
    phi = C0 + C3*t^3 + C4*t^4 + C5*t^5;
    theta_in = Angle_D4_Abs + t * Input_Stroke;
    p_abs = O1 + r_rot * [cosd(theta_in), sind(theta_in)];

    % 坐标变换，从绝对坐标转到相对坐标
    p_rel_x = p_abs(1)*cosd(-phi) - p_abs(2)*sind(-phi);
    p_rel_y = p_abs(1)*sind(-phi) + p_abs(2)*cosd(-phi);

    idx = num_gen + i;
    Raw_Radius(idx) = sqrt(p_rel_x^2 + p_rel_y^2);
    Raw_Angle(idx) = atan2(p_rel_y, p_rel_x);
    Raw_Type(idx) = 2;
end

% 去掉重复的点，建立查找表
[LUT_Radius, unique_idx] = unique(Raw_Radius);
LUT_Angle = unwrap(Raw_Angle(unique_idx));
LUT_Type = Raw_Type(unique_idx);

% 开始仿真，算出整个运动过程
dt = 0.0002;  % 时间步长，不能太大否则不准
Th_Start_Abs = -108;
Th_End_Abs = 0 + 72;
Total_Time = deg2rad(Th_End_Abs - Th_Start_Abs) / Omega;
Time = 0:dt:Total_Time;
Theta_Pin_Seq = Th_Start_Abs + rad2deg(Omega * Time);
Phi_Out = zeros(size(Time));

for i = 1:length(Time)
    theta_pin = Theta_Pin_Seq(i);
    Lx = O1(1) + r_rot * cosd(theta_pin);
    Ly = O1(2) + r_rot * sind(theta_pin);
    Current_R = sqrt(Lx^2 + Ly^2);

    % 限制在查找表范围内
    Current_R = max(min(Current_R, max(LUT_Radius)), min(LUT_Radius));
    Rel_Ang_Rad = interp1(LUT_Radius, LUT_Angle, Current_R, 'pchip');

    Pin_Abs_Rad = atan2(Ly, Lx);
    diff_rad = Pin_Abs_Rad - Rel_Ang_Rad;
    diff_rad = atan2(sin(diff_rad), cos(diff_rad));
    Phi_Out(i) = rad2deg(diff_rad);
end

% 处理一下数据，准备算目标函数
Plot_X = Theta_Pin_Seq - Theta_Pin_Seq(1);
Phi_Plot = Phi_start - Phi_Out;

% 用差分算速度和加速度
Vel = gradient(deg2rad(Phi_Plot)) ./ gradient(Time);
Acc = gradient(Vel) ./ gradient(Time);

% 平滑一下，去掉噪声
Vel_S = smoothdata(Vel, 'gaussian', 15);
Acc_S = smoothdata(Acc, 'gaussian', 40);

% ========== 下面开始算目标函数 ==========
% 思路就是：哪里不好就加惩罚，惩罚越多说明越差

cost = 0;

% 第一条：位移不能倒退
% 槽轮只能往一个方向转，要是出现倒退那肯定不行
disp_diff = diff(Phi_Plot);
negative_disp = disp_diff(disp_diff < 0);
if ~isempty(negative_disp)
    % 倒退了就狠狠惩罚
    cost = cost + 1000 * sum(abs(negative_disp)) + 500 * length(negative_disp);
end

% 第二条：速度要平滑
% 速度突变会产生冲击，机构容易坏
vel_jerk = diff(Vel_S);
vel_jerk_penalty = sum(vel_jerk.^2) * dt * 0.001;
cost = cost + vel_jerk_penalty;

% 108度附近速度要小一点
idx_108 = find(Plot_X >= 108 - 5 & Plot_X <= 108 + 5);
if ~isempty(idx_108)
    vel_at_108 = mean(abs(Vel_S(idx_108)));
    cost = cost + 10 * vel_at_108;
end

% 速度曲线不能有太多尖峰
[~, vel_peaks] = findpeaks(abs(Vel_S), 'MinPeakProminence', 0.5);
if length(vel_peaks) > 2
    cost = cost + 50 * (length(vel_peaks) - 2);
end

% 第三条：加速度也要平滑
% 加速度突变就是所谓的"刚性冲击"，很伤机构
acc_jerk = diff(Acc_S);
acc_jerk_penalty = sum(acc_jerk.^2) * dt * 1e-6;
cost = cost + acc_jerk_penalty;

% 加速度峰值不能太大
max_acc = max(abs(Acc_S));
if max_acc > 500
    cost = cost + 0.1 * (max_acc - 500);
end

% 加速度曲线也不能有太多尖峰
[~, acc_peaks] = findpeaks(abs(Acc_S), 'MinPeakProminence', 50);
if length(acc_peaks) > 3
    cost = cost + 30 * (length(acc_peaks) - 3);
end

% 第四条：整体要光滑
% 位移曲线的曲率不能太大
disp_curv = diff(diff(Phi_Plot));
cost = cost + 0.01 * sum(disp_curv.^2);

% 起点和终点速度要接近0，这样过渡才平滑
vel_start = abs(Vel_S(1));
vel_end = abs(Vel_S(end));
cost = cost + 5 * (vel_start + vel_end);

% 加速度起止点也要小
acc_start = abs(Acc_S(1));
acc_end = abs(Acc_S(end));
cost = cost + 0.05 * (acc_start + acc_end);

% 如果需要返回详细结果，就把数据打包
if return_results
    results.Time = Time;
    results.Plot_X = Plot_X;
    results.Phi_Plot = Phi_Plot;
    results.Vel_S = Vel_S;
    results.Acc_S = Acc_S;
    results.Vel = Vel;
    results.Acc = Acc;
    results.params = params;
    results.cost = cost;
else
    results = [];
end

end
