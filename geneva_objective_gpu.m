function [cost, results] = geneva_objective_gpu(params, return_results, use_gpu)
% GENEVA_OBJECTIVE_GPU 槽轮机构优化目标函数 - GPU加速版本
%
% 输入参数:
%   params(1) = Angle_G_rel  : G点偏转角 (2~10°)
%   params(2) = Input_Stroke : 输入行程 (25~45°)
%   params(3) = v_E          : 速度边界条件 (-60~-5)
%   params(4) = a_E          : 加速度边界条件 (1~30)
%   return_results           : 是否返回详细结果 (默认false)
%   use_gpu                  : 是否使用GPU (默认自动检测)

if nargin < 2, return_results = false; end
if nargin < 3
    try
        gpuDevice;
        use_gpu = true;
    catch
        use_gpu = false;
    end
end

%% 解析优化变量
Angle_G_rel = params(1);
Input_Stroke = params(2);
v_E = params(3);
a_E = params(4);

%% 固定参数
R = 55;
N = 10;
Omega_RPM = 60;
Omega = Omega_RPM * 2*pi/60;

%% 几何计算
d = 2 * R * sin(pi/N);
r_rot = d/2;
O1 = [-R, 0];
Phi_end = -Angle_G_rel;
Phi_start = 0;

%% 五次多项式系数
C0 = Phi_start;
Delta_Phi = Phi_end - C0;
M = [1 1 1; 3 4 5; 6 12 20];
V = [Delta_Phi; v_E; a_E];
X = M \ V;
C3 = X(1); C4 = X(2); C5 = X(3);

%% 生成几何路书 (向量化 + GPU)
num_gen = 5000;

% 转换为GPU数组（如果启用）
if use_gpu
    theta_dwell = gpuArray(linspace(-Input_Stroke - 200, -Input_Stroke, num_gen));
    tau_gen = gpuArray(linspace(0, 1, num_gen));
else
    theta_dwell = linspace(-Input_Stroke - 200, -Input_Stroke, num_gen);
    tau_gen = linspace(0, 1, num_gen);
end

% A. CD4段 - 向量化计算
p_abs_x1 = O1(1) + r_rot * cosd(theta_dwell);
p_abs_y1 = O1(2) + r_rot * sind(theta_dwell);
Raw_Radius1 = sqrt(p_abs_x1.^2 + p_abs_y1.^2);
Raw_Angle1 = atan2(p_abs_y1, p_abs_x1);
Raw_Type1 = ones(1, num_gen);

% B. D4G段 - 向量化计算
phi_vec = C0 + C3*tau_gen.^3 + C4*tau_gen.^4 + C5*tau_gen.^5;
theta_in = -Input_Stroke + tau_gen * Input_Stroke;
p_abs_x2 = O1(1) + r_rot * cosd(theta_in);
p_abs_y2 = O1(2) + r_rot * sind(theta_in);

p_rel_x = p_abs_x2.*cosd(-phi_vec) - p_abs_y2.*sind(-phi_vec);
p_rel_y = p_abs_x2.*sind(-phi_vec) + p_abs_y2.*cosd(-phi_vec);

Raw_Radius2 = sqrt(p_rel_x.^2 + p_rel_y.^2);
Raw_Angle2 = atan2(p_rel_y, p_rel_x);
Raw_Type2 = 2 * ones(1, num_gen);

% 合并
Raw_Radius = [Raw_Radius1, Raw_Radius2];
Raw_Angle = [Raw_Angle1, Raw_Angle2];
Raw_Type = [Raw_Type1, Raw_Type2];

% 从GPU取回数据进行unique操作
if use_gpu
    Raw_Radius = gather(Raw_Radius);
    Raw_Angle = gather(Raw_Angle);
end

[LUT_Radius, unique_idx] = unique(Raw_Radius);
LUT_Angle = unwrap(Raw_Angle(unique_idx));
LUT_Type = Raw_Type(unique_idx);

%% 物理仿真 - 向量化
dt = 0.0002;
Th_Start_Abs = -108;
Th_End_Abs = 72;
Total_Time = deg2rad(Th_End_Abs - Th_Start_Abs) / Omega;

if use_gpu
    Time = gpuArray(0:dt:Total_Time);
else
    Time = 0:dt:Total_Time;
end

Theta_Pin_Seq = Th_Start_Abs + rad2deg(Omega * Time);

% 向量化计算
Lx = O1(1) + r_rot * cosd(Theta_Pin_Seq);
Ly = O1(2) + r_rot * sind(Theta_Pin_Seq);
Current_R = sqrt(Lx.^2 + Ly.^2);

% 限制范围
Current_R = max(min(Current_R, max(LUT_Radius)), min(LUT_Radius));

% 取回GPU数据进行插值
if use_gpu
    Current_R = gather(Current_R);
    Lx = gather(Lx);
    Ly = gather(Ly);
    Time = gather(Time);
    Theta_Pin_Seq = gather(Theta_Pin_Seq);
end

Rel_Ang_Rad = interp1(LUT_Radius, LUT_Angle, Current_R, 'pchip');
Pin_Abs_Rad = atan2(Ly, Lx);
diff_rad = Pin_Abs_Rad - Rel_Ang_Rad;
diff_rad = atan2(sin(diff_rad), cos(diff_rad));
Phi_Out = rad2deg(diff_rad);

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

% 1. 角位移单调性惩罚
disp_diff = diff(Phi_Plot);
negative_disp = disp_diff(disp_diff < 0);
if ~isempty(negative_disp)
    cost = cost + 1000 * sum(abs(negative_disp)) + 500 * length(negative_disp);
end

% 2. 角速度平滑性惩罚
vel_jerk = diff(Vel_S);
vel_jerk_penalty = sum(vel_jerk.^2) * dt * 0.001;
cost = cost + vel_jerk_penalty;

% 108°处速度惩罚
idx_108 = find(Plot_X >= 103 & Plot_X <= 113);
if ~isempty(idx_108)
    vel_at_108 = mean(abs(Vel_S(idx_108)));
    cost = cost + 10 * vel_at_108;
end

% 速度尖点惩罚
[~, vel_peaks] = findpeaks(abs(Vel_S), 'MinPeakProminence', 0.5);
if length(vel_peaks) > 2
    cost = cost + 50 * (length(vel_peaks) - 2);
end

% 3. 角加速度平滑性惩罚
acc_jerk = diff(Acc_S);
acc_jerk_penalty = sum(acc_jerk.^2) * dt * 1e-6;
cost = cost + acc_jerk_penalty;

% 加速度峰值惩罚
max_acc = max(abs(Acc_S));
if max_acc > 500
    cost = cost + 0.1 * (max_acc - 500);
end

% 加速度尖点惩罚
[~, acc_peaks] = findpeaks(abs(Acc_S), 'MinPeakProminence', 50);
if length(acc_peaks) > 3
    cost = cost + 30 * (length(acc_peaks) - 3);
end

% 4. 整体曲线光滑度
disp_curv = diff(diff(Phi_Plot));
cost = cost + 0.01 * sum(disp_curv.^2);

% 起止点平滑
vel_start = abs(Vel_S(1));
vel_end = abs(Vel_S(end));
cost = cost + 5 * (vel_start + vel_end);

acc_start = abs(Acc_S(1));
acc_end = abs(Acc_S(end));
cost = cost + 0.05 * (acc_start + acc_end);

%% 返回结果
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
