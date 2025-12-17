clc; clear;

%% 1. --- 用户参数 (User Inputs) ---
R = 55; % 中心距
N = 10; % 槽数
r_pin = 3; % 圆销半径

% [几何变量] G 点偏转角
Angle_G_rel = 7.5; 

% [时间变量] 输入行程 (D4->E)
Input_Stroke = 40; 

% [运动变量] 优化旋钮
v_E = -30; 
a_E = 9;

% 仿真设置
Omega_RPM = 60; 
Omega = Omega_RPM * 2*pi/60;

%% 2. --- 几何与边界自动计算 ---
d = 2 * R * sin(pi/N); r_rot = d/2;
O1 = [-R, 0]; O = [0, 0];
Angle_Pitch = 360 / N; 
Phi_end = -Angle_G_rel;
Phi_start = 0; 

% 五次多项式求解
C0 = Phi_start; Delta_Phi = Phi_end - C0;
M = [1 1 1; 3 4 5; 6 12 20]; V = [Delta_Phi; v_E; a_E];
X = M \ V ;
C3 = X(1); C4 = X(2); C5 = X(3);

%% 3. --- 生成几何路书 (Lookup Table) ---
num_gen = 10000;
Raw_Radius = []; Raw_Angle = []; Raw_Type = []; 

% A. 生成 CD4 段
Angle_D4_Abs = -Input_Stroke;
Angle_Extend = Angle_D4_Abs - 200; 
theta_dwell = linspace(Angle_Extend, Angle_D4_Abs, num_gen);
for i = 1:num_gen
    p_abs = O1 + r_rot * [cosd(theta_dwell(i)), sind(theta_dwell(i))];
    Raw_Radius(end+1) = sqrt(p_abs(1)^2 + p_abs(2)^2);
    Raw_Angle(end+1) = atan2(p_abs(2), p_abs(1));
    Raw_Type(end+1) = 1; 
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
    
    Raw_Radius(end+1) = sqrt(p_rel_x^2 + p_rel_y^2);
    Raw_Angle(end+1) = atan2(p_rel_y, p_rel_x);
    Raw_Type(end+1) = 2; 
end
[LUT_Radius, unique_idx] = unique(Raw_Radius);
LUT_Angle = unwrap(Raw_Angle(unique_idx));
LUT_Type = Raw_Type(unique_idx);

%% 4. --- 物理仿真 (C -> A) ---
dt = 0.0001;
Th_Start_Abs = -108; 
Th_End_Abs = 0 + 72; 
Total_Time = deg2rad(Th_End_Abs - Th_Start_Abs) / Omega;
Time = 0:dt:Total_Time;
Theta_Pin_Seq = Th_Start_Abs + rad2deg(Omega * Time);
Phi_Out = zeros(size(Time));
Region_Flag = zeros(size(Time));

for i = 1:length(Time)
    theta_pin = Theta_Pin_Seq(i);
    Lx = O1(1) + r_rot * cosd(theta_pin);
    Ly = O1(2) + r_rot * sind(theta_pin);
    Current_R = sqrt(Lx^2 + Ly^2);
    
    Current_R = max(min(Current_R, max(LUT_Radius)), min(LUT_Radius));
    Rel_Ang_Rad = interp1(LUT_Radius, LUT_Angle, Current_R, 'pchip');
    Region_Flag(i) = interp1(LUT_Radius, LUT_Type, Current_R, 'nearest');
    
    Pin_Abs_Rad = atan2(Ly, Lx);
    diff_rad = Pin_Abs_Rad - Rel_Ang_Rad;
    diff_rad = atan2(sin(diff_rad), cos(diff_rad));
    Phi_Out(i) = rad2deg(diff_rad);
end

%% 5. --- 后处理 ---
Plot_X = Theta_Pin_Seq - Theta_Pin_Seq(1);
X_A = Plot_X(end); 
Phi_Plot = Phi_start - Phi_Out;
Vel = gradient(deg2rad(Phi_Plot)) ./ gradient(Time);
Acc = gradient(Vel) ./ gradient(Time);
Vel_S = smoothdata(Vel, 'gaussian', 20);
Acc_S = smoothdata(Acc, 'gaussian', 50);
Acc_A = Acc_S(end);

%% 6. --- 绘图 (修改区：黑->蓝->红->绿) ---

% 1. 查找或创建窗口
fig_handle = findobj('Type', 'figure', 'Name', '全周期运动分析对比');
if isempty(fig_handle)
    % 如果窗口不存在，说明是第1次运行
    fig_handle = figure('Color', 'w', 'Position', [100, 50, 800, 900], 'Name', '全周期运动分析对比');
    run_count = 0;
else
    figure(fig_handle); % 激活窗口
    % 2. 核心魔法：数一数图上已经画了几条线
    ax_check = subplot(3,1,1);
    existing_lines = findobj(ax_check, 'Type', 'line');
    run_count = length(existing_lines);
end

% 3. 定义颜色序列 (黑 -> 蓝 -> 红 -> 绿)
switch run_count
    case 0
        This_Color = [0 0 0];       % 1. 黑色 (Black)
        Color_Name = '黑色';
    case 1
        This_Color = [0 0 1];       % 2. 蓝色 (Blue)
        Color_Name = '蓝色';
    case 2
        This_Color = [1 0 0];       % 3. 红色 (Red)
        Color_Name = '红色';
    case 3
        This_Color = [0 0.6 0];     % 4. 绿色 (Green)
        Color_Name = '绿色';
    otherwise
        This_Color = rand(1, 3)*0.7; % 5+. 随机颜色
        Color_Name = '随机色';
end

fprintf('>> 第 %d 次叠加运行，使用颜色: %s\n', run_count+1, Color_Name);

% --- 绘图逻辑 (强制 hold on) ---

% 1. 位移
subplot(3,1,1); hold on;
plot(Plot_X, Phi_Plot, 'Color', This_Color, 'LineWidth', 2);
ylabel('角位移 (deg)'); title('槽轮角位移'); grid on; xlim([0, X_A+5]);

% 2. 速度
subplot(3,1,2); hold on;
plot(Plot_X, Vel_S, 'Color', This_Color, 'LineWidth', 2);
ylabel('角速度 (rad/s)'); title('槽轮角速度'); grid on; xlim([0, X_A+5]);

% 3. 加速度
subplot(3,1,3); hold on;
plot(Plot_X, Acc_S, 'Color', This_Color, 'LineWidth', 2);
ylabel('角加速度 (rad/s^2)'); xlabel('圆销转角(°)'); title('槽轮角加速度'); grid on; xlim([0, X_A+5]);

% 标注 A 点数值 (带背景框)
%text(X_A, Acc_A, sprintf(' Acc: %.2f', Acc_A), ...
  %  'Color', This_Color, 'FontWeight', 'bold', ...
   % 'BackgroundColor', 'w', 'EdgeColor', This_Color, 'Margin', 1);