function [cost, results] = geneva_objective(params, return_results)
% GENEVA_OBJECTIVE 槽轮机构五次多项式运动规律优化目标函数
%
% 输入参数:
%   params(1) = Angle_G_rel  : G点偏转角 (2~10°)
%   params(2) = Input_Stroke : 输入行程 (25~45°)
%   params(3) = v_E          : 速度边界条件 (-60~-5)
%   params(4) = a_E          : 加速度边界条件 (1~30)
%   return_results           : 是否返回详细结果 (可选, 默认false)
%
% 输出:
%   cost    : 目标函数值 (越小越好)
%   results : 包含运动曲线的结构体 (仅当return_results=true时有效)

if nargin < 2
    return_results = false;
end

%% 解析优化变量
Angle_G_rel = params(1);
Input_Stroke = params(2);
v_E = params(3);
a_E = params(4);

%% 固定参数
R = 55;           % 中心距
N = 10;           % 槽数
r_pin = 3;        % 圆销半径
Omega_RPM = 60;   % 转速
Omega = Omega_RPM * 2*pi/60;

%% 几何与边界自动计算
d = 2 * R * sin(pi/N);
r_rot = d/2;
O1 = [-R, 0];
Angle_Pitch = 360 / N;
Phi_end = -Angle_G_rel;
Phi_start = 0;

%% 五次多项式系数求解
C0 = Phi_start;
Delta_Phi = Phi_end - C0;
M = [1 1 1; 3 4 5; 6 12 20];
V = [Delta_Phi; v_E; a_E];
X = M \ V;
C3 = X(1); C4 = X(2); C5 = X(3);

%% 生成几何路书 (Lookup Table)
num_gen = 5000;  % 减少点数以加速计算
Raw_Radius = zeros(1, 2*num_gen);
Raw_Angle = zeros(1, 2*num_gen);
Raw_Type = zeros(1, 2*num_gen);

% A. 生成 CD4 段
Angle_D4_Abs = -Input_Stroke;
Angle_Extend = Angle_D4_Abs - 200;
theta_dwell = linspace(Angle_Extend, Angle_D4_Abs, num_gen);

for i = 1:num_gen
    p_abs = O1 + r_rot * [cosd(theta_dwell(i)), sind(theta_dwell(i))];
    Raw_Radius(i) = sqrt(p_abs(1)^2 + p_abs(2)^2);
    Raw_Angle(i) = atan2(p_abs(2), p_abs(1));
    Raw_Type(i) = 1;
end

% B. 生成 D4G 段
tau_gen = linspace(0, 1, num_gen);
for i = 1:num_gen
    t = tau_gen(i);
    phi = C0 + C3*t^3 + C4*t^4 + C5*t^5;
    theta_in = Angle_D4_Abs + t * Input_Stroke;
    p_abs = O1 + r_rot * [cosd(theta_in), sind(theta_in)];

    p_rel_x = p_abs(1)*cosd(-phi) - p_abs(2)*sind(-phi);
    p_rel_y = p_abs(1)*sind(-phi) + p_abs(2)*cosd(-phi);

    idx = num_gen + i;
    Raw_Radius(idx) = sqrt(p_rel_x^2 + p_rel_y^2);
    Raw_Angle(idx) = atan2(p_rel_y, p_rel_x);
    Raw_Type(idx) = 2;
end

[LUT_Radius, unique_idx] = unique(Raw_Radius);
LUT_Angle = unwrap(Raw_Angle(unique_idx));
LUT_Type = Raw_Type(unique_idx);

%% 物理仿真 (C -> A)
dt = 0.0002;  % 适当增大步长以加速
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

    Current_R = max(min(Current_R, max(LUT_Radius)), min(LUT_Radius));
    Rel_Ang_Rad = interp1(LUT_Radius, LUT_Angle, Current_R, 'pchip');

    Pin_Abs_Rad = atan2(Ly, Lx);
    diff_rad = Pin_Abs_Rad - Rel_Ang_Rad;
    diff_rad = atan2(sin(diff_rad), cos(diff_rad));
    Phi_Out(i) = rad2deg(diff_rad);
end

%% 后处理
Plot_X = Theta_Pin_Seq - Theta_Pin_Seq(1);
Phi_Plot = Phi_start - Phi_Out;

% 计算速度和加速度
Vel = gradient(deg2rad(Phi_Plot)) ./ gradient(Time);
Acc = gradient(Vel) ./ gradient(Time);

% 平滑处理
Vel_S = smoothdata(Vel, 'gaussian', 15);
Acc_S = smoothdata(Acc, 'gaussian', 40);

%% ============ 目标函数计算 ============

cost = 0;

% ---------- 1. 角位移单调性惩罚 ----------
% 计算位移的导数，惩罚任何负值（倒退）
disp_diff = diff(Phi_Plot);
negative_disp = disp_diff(disp_diff < 0);
if ~isempty(negative_disp)
    % 惩罚倒退的幅度和次数
    cost = cost + 1000 * sum(abs(negative_disp)) + 500 * length(negative_disp);
end

% ---------- 2. 角速度平滑性惩罚 ----------
% 2.1 速度连续性 - 惩罚速度的突变（二阶导数）
vel_jerk = diff(Vel_S);
vel_jerk_penalty = sum(vel_jerk.^2) * dt * 0.001;
cost = cost + vel_jerk_penalty;

% 2.2 在108°处速度尽量小
% 找到108°附近的索引
idx_108 = find(Plot_X >= 108 - 5 & Plot_X <= 108 + 5);
if ~isempty(idx_108)
    vel_at_108 = mean(abs(Vel_S(idx_108)));
    cost = cost + 10 * vel_at_108;
end

% 2.3 速度尖点惩罚 - 检测速度曲线的局部极值点过多
[~, vel_peaks] = findpeaks(abs(Vel_S), 'MinPeakProminence', 0.5);
if length(vel_peaks) > 2
    cost = cost + 50 * (length(vel_peaks) - 2);
end

% ---------- 3. 角加速度平滑性惩罚 ----------
% 3.1 加速度的突变（加加速度/jerk）
acc_jerk = diff(Acc_S);
acc_jerk_penalty = sum(acc_jerk.^2) * dt * 1e-6;
cost = cost + acc_jerk_penalty;

% 3.2 加速度峰值惩罚 - 过大的加速度不好
max_acc = max(abs(Acc_S));
if max_acc > 500
    cost = cost + 0.1 * (max_acc - 500);
end

% 3.3 加速度尖点检测
[~, acc_peaks] = findpeaks(abs(Acc_S), 'MinPeakProminence', 50);
if length(acc_peaks) > 3
    cost = cost + 30 * (length(acc_peaks) - 3);
end

% ---------- 4. 整体曲线光滑度 ----------
% 4.1 位移曲线的二阶导数平方和（曲率惩罚）
disp_curv = diff(diff(Phi_Plot));
cost = cost + 0.01 * sum(disp_curv.^2);

% 4.2 起点和终点的平滑过渡
% 检查运动开始和结束时的速度是否接近0
vel_start = abs(Vel_S(1));
vel_end = abs(Vel_S(end));
cost = cost + 5 * (vel_start + vel_end);

% 4.3 加速度起止点平滑
acc_start = abs(Acc_S(1));
acc_end = abs(Acc_S(end));
cost = cost + 0.05 * (acc_start + acc_end);

%% 返回结果（如果需要）
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
