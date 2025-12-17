%% 槽轮机构五次多项式运动规律 - 遗传算法优化
% 功能: 使用遗传算法优化槽轮机构参数，使运动曲线更加平顺
% 支持: GPU加速 + CPU回退
%
% 优化变量:
%   1. Angle_G_rel  : G点偏转角 (2~10°)
%   2. Input_Stroke : 输入行程 (25~45°)
%   3. v_E          : 速度边界条件 (-60~-5)
%   4. a_E          : 加速度边界条件 (1~30)

clc; clear; close all;

fprintf('============================================\n');
fprintf('   槽轮机构遗传算法优化程序\n');
fprintf('============================================\n\n');

%% ==================== 1. 检测GPU可用性 ====================
gpu_available = false;
try
    gpu_info = gpuDevice;
    gpu_available = true;
    fprintf('[GPU] 检测到GPU: %s\n', gpu_info.Name);
    fprintf('[GPU] 计算能力: %.1f, 显存: %.1f GB\n', ...
        gpu_info.ComputeCapability, gpu_info.TotalMemory/1e9);
catch
    fprintf('[CPU] 未检测到可用GPU，使用CPU模式\n');
end
fprintf('\n');

%% ==================== 2. 优化参数设置 ====================
% 优化变量边界
%        [Angle_G_rel, Input_Stroke, v_E,   a_E]
lb =     [2,           25,           -60,   1];    % 下界
ub =     [10,          45,           -5,    30];   % 上界

% 遗传算法参数 (参考表8)
ga_options = struct();
ga_options.PopulationSize = 50;      % 种群规模
ga_options.MaxGenerations = 100;     % 迭代次数
ga_options.CrossoverFraction = 0.9;  % 交叉概率
ga_options.MutationRate = 0.05;      % 变异概率
ga_options.EliteCount = 2;           % 精英个体数

% 运行次数
num_runs = 4;

fprintf('遗传算法参数设定:\n');
fprintf('  种群规模: %d\n', ga_options.PopulationSize);
fprintf('  迭代次数: %d\n', ga_options.MaxGenerations);
fprintf('  交叉概率: %.2f\n', ga_options.CrossoverFraction);
fprintf('  变异概率: %.2f\n', ga_options.MutationRate);
fprintf('  优化次数: %d\n', num_runs);
fprintf('\n');

%% ==================== 3. 定义目标函数 ====================
% 目标函数包装器
fitness_func = @(x) geneva_objective(x, false);

%% ==================== 4. 配置GA选项 ====================
% 基础选项
options = optimoptions('ga', ...
    'PopulationSize', ga_options.PopulationSize, ...
    'MaxGenerations', ga_options.MaxGenerations, ...
    'CrossoverFraction', ga_options.CrossoverFraction, ...
    'EliteCount', ga_options.EliteCount, ...
    'FunctionTolerance', 1e-6, ...
    'Display', 'iter', ...
    'PlotFcn', [], ...  % 稍后手动绘制
    'UseParallel', true);  % 启用并行计算

% 如果GPU可用，尝试使用GPU加速
if gpu_available
    try
        % 启用并行池
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool('local');
        end
        fprintf('[并行] 已启用并行计算，工作进程数: %d\n\n', pool.NumWorkers);
    catch
        fprintf('[警告] 无法启用并行计算，使用串行模式\n\n');
        options = optimoptions(options, 'UseParallel', false);
    end
else
    % CPU模式也尝试使用并行
    try
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool('local');
        end
        fprintf('[并行] CPU并行计算，工作进程数: %d\n\n', pool.NumWorkers);
    catch
        fprintf('[串行] 使用串行计算模式\n\n');
        options = optimoptions(options, 'UseParallel', false);
    end
end

%% ==================== 5. 多次运行优化 ====================
all_results = cell(num_runs, 1);
all_best_params = zeros(num_runs, 4);
all_best_costs = zeros(num_runs, 1);
all_history = cell(num_runs, 1);

fprintf('开始优化...\n');
fprintf('============================================\n\n');

for run_idx = 1:num_runs
    fprintf('>>> 第 %d/%d 次优化运行 <<<\n', run_idx, num_runs);
    fprintf('--------------------------------------------\n');

    % 记录优化历史
    history = struct('best_cost', [], 'generation', []);

    % 自定义输出函数记录历史
    options = optimoptions(options, 'OutputFcn', @(options, state, flag) ...
        ga_output_func(options, state, flag, run_idx));

    % 运行遗传算法
    tic;
    [best_params, best_cost, exitflag, output] = ga(fitness_func, 4, ...
        [], [], [], [], lb, ub, [], options);
    elapsed_time = toc;

    % 保存结果
    all_best_params(run_idx, :) = best_params;
    all_best_costs(run_idx) = best_cost;

    % 获取详细结果用于绘图
    [~, results] = geneva_objective(best_params, true);
    all_results{run_idx} = results;

    fprintf('\n第 %d 次优化完成!\n', run_idx);
    fprintf('  最优目标值: %.4f\n', best_cost);
    fprintf('  优化参数:\n');
    fprintf('    Angle_G_rel  = %.2f°\n', best_params(1));
    fprintf('    Input_Stroke = %.2f°\n', best_params(2));
    fprintf('    v_E          = %.2f\n', best_params(3));
    fprintf('    a_E          = %.2f\n', best_params(4));
    fprintf('  耗时: %.2f 秒\n', elapsed_time);
    fprintf('--------------------------------------------\n\n');
end

%% ==================== 6. 汇总结果 ====================
fprintf('\n============================================\n');
fprintf('优化结果汇总 (表9格式)\n');
fprintf('============================================\n');
fprintf('%-15s', '优化参数');
for i = 1:num_runs
    fprintf('第%d次\t\t', i);
end
fprintf('\n');
fprintf('------------------------------------------------------------\n');

param_names = {'s/cm', 'h1/cm', 'h2/cm', 'h3/cm'};
display_names = {'Angle_G_rel', 'Input_Stroke', 'v_E', 'a_E'};

for p = 1:4
    fprintf('%-15s', display_names{p});
    for i = 1:num_runs
        fprintf('%.2f\t\t', all_best_params(i, p));
    end
    fprintf('\n');
end
fprintf('------------------------------------------------------------\n');
fprintf('%-15s', '目标函数值');
for i = 1:num_runs
    fprintf('%.2f\t\t', all_best_costs(i));
end
fprintf('\n');

% 找出最优结果
[global_best_cost, best_run_idx] = min(all_best_costs);
global_best_params = all_best_params(best_run_idx, :);

fprintf('\n============================================\n');
fprintf('全局最优解 (第 %d 次运行)\n', best_run_idx);
fprintf('============================================\n');
fprintf('  Angle_G_rel  = %.2f°\n', global_best_params(1));
fprintf('  Input_Stroke = %.2f°\n', global_best_params(2));
fprintf('  v_E          = %.2f\n', global_best_params(3));
fprintf('  a_E          = %.2f\n', global_best_params(4));
fprintf('  目标函数值   = %.4f\n', global_best_cost);
fprintf('============================================\n');

%% ==================== 7. 绘制优化后的运动曲线 ====================
fprintf('\n正在生成运动曲线图...\n');

% 获取最优结果
best_result = all_results{best_run_idx};

% 创建运动曲线图
fig1 = figure('Color', 'w', 'Position', [100, 50, 800, 900], ...
    'Name', '优化后运动曲线');

% 角位移
subplot(3,1,1);
plot(best_result.Plot_X, best_result.Phi_Plot, 'k-', 'LineWidth', 2);
ylabel('角位移 (deg)');
title('槽轮角位移');
grid on;
xlim([0, max(best_result.Plot_X)+5]);

% 角速度
subplot(3,1,2);
plot(best_result.Plot_X, best_result.Vel_S, 'k-', 'LineWidth', 2);
ylabel('角速度 (rad/s)');
title('槽轮角速度');
grid on;
xlim([0, max(best_result.Plot_X)+5]);

% 角加速度
subplot(3,1,3);
plot(best_result.Plot_X, best_result.Acc_S, 'k-', 'LineWidth', 2);
ylabel('角加速度 (rad/s^2)');
xlabel('圆销转角(°)');
title('槽轮角加速度');
grid on;
xlim([0, max(best_result.Plot_X)+5]);

% 保存图像
saveas(fig1, 'optimized_motion_curves.png');
fprintf('运动曲线图已保存: optimized_motion_curves.png\n');

%% ==================== 8. 绘制GA收敛曲线 ====================
fprintf('正在生成GA收敛曲线...\n');

fig2 = figure('Color', 'w', 'Position', [150, 100, 1000, 600], ...
    'Name', 'GA优化收敛曲线');

% 2x2布局显示4次优化的收敛曲线
for i = 1:num_runs
    subplot(2, 2, i);

    % 读取保存的历史数据
    hist_file = sprintf('ga_history_run%d.mat', i);
    if exist(hist_file, 'file')
        load(hist_file, 'ga_history');
        plot(ga_history.generation, ga_history.best_cost, 'k-', 'LineWidth', 1.5);
        xlabel('迭代次数');
        ylabel('最优解');
        title(sprintf('(%c) 第%d次', char('a'+i-1), i));
        grid on;
        xlim([0 100]);
    end
end

% 保存图像
saveas(fig2, 'ga_convergence_curves.png');
fprintf('收敛曲线图已保存: ga_convergence_curves.png\n');

%% ==================== 9. 对比优化前后的结果 ====================
fprintf('\n正在生成优化前后对比图...\n');

% 原始参数 (从my.m中提取)
original_params = [7.5, 40, -30, 9];
[~, original_result] = geneva_objective(original_params, true);

fig3 = figure('Color', 'w', 'Position', [200, 100, 1000, 800], ...
    'Name', '优化前后对比');

% 角位移对比
subplot(3,1,1);
h1 = plot(original_result.Plot_X, original_result.Phi_Plot, 'b--', 'LineWidth', 1.5);
hold on;
h2 = plot(best_result.Plot_X, best_result.Phi_Plot, 'r-', 'LineWidth', 2);
hold off;
ylabel('角位移 (deg)');
title('槽轮角位移对比');
legend([h1, h2], {'优化前', '优化后'}, 'Location', 'best');
grid on;
xlim([0, max(best_result.Plot_X)+5]);

% 角速度对比
subplot(3,1,2);
h1 = plot(original_result.Plot_X, original_result.Vel_S, 'b--', 'LineWidth', 1.5);
hold on;
h2 = plot(best_result.Plot_X, best_result.Vel_S, 'r-', 'LineWidth', 2);
hold off;
ylabel('角速度 (rad/s)');
title('槽轮角速度对比');
legend([h1, h2], {'优化前', '优化后'}, 'Location', 'best');
grid on;
xlim([0, max(best_result.Plot_X)+5]);

% 角加速度对比
subplot(3,1,3);
h1 = plot(original_result.Plot_X, original_result.Acc_S, 'b--', 'LineWidth', 1.5);
hold on;
h2 = plot(best_result.Plot_X, best_result.Acc_S, 'r-', 'LineWidth', 2);
hold off;
ylabel('角加速度 (rad/s^2)');
xlabel('圆销转角(°)');
title('槽轮角加速度对比');
legend([h1, h2], {'优化前', '优化后'}, 'Location', 'best');
grid on;
xlim([0, max(best_result.Plot_X)+5]);

% 保存图像
saveas(fig3, 'optimization_comparison.png');
fprintf('对比图已保存: optimization_comparison.png\n');

%% ==================== 10. 保存优化结果 ====================
fprintf('\n正在保存优化结果...\n');

% 保存到MAT文件
save('optimization_results.mat', ...
    'all_best_params', 'all_best_costs', 'all_results', ...
    'global_best_params', 'global_best_cost', 'best_run_idx', ...
    'ga_options', 'lb', 'ub');

% 导出结果表格到文本文件
fid = fopen('optimization_results.txt', 'w');
fprintf(fid, '槽轮机构遗传算法优化结果\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, '遗传算法参数:\n');
fprintf(fid, '  种群规模: %d\n', ga_options.PopulationSize);
fprintf(fid, '  迭代次数: %d\n', ga_options.MaxGenerations);
fprintf(fid, '  交叉概率: %.2f\n', ga_options.CrossoverFraction);
fprintf(fid, '  变异概率: %.2f\n\n', ga_options.MutationRate);

fprintf(fid, '优化结果 (表9格式):\n');
fprintf(fid, '%-15s', '优化参数');
for i = 1:num_runs
    fprintf(fid, '第%d次\t', i);
end
fprintf(fid, '\n');
fprintf(fid, '----------------------------------------\n');

for p = 1:4
    fprintf(fid, '%-15s', display_names{p});
    for i = 1:num_runs
        fprintf(fid, '%.2f\t', all_best_params(i, p));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\n');

fprintf(fid, '全局最优解 (第 %d 次运行):\n', best_run_idx);
fprintf(fid, '  Angle_G_rel  = %.2f°\n', global_best_params(1));
fprintf(fid, '  Input_Stroke = %.2f°\n', global_best_params(2));
fprintf(fid, '  v_E          = %.2f\n', global_best_params(3));
fprintf(fid, '  a_E          = %.2f\n', global_best_params(4));
fprintf(fid, '  目标函数值   = %.4f\n', global_best_cost);
fclose(fid);

fprintf('结果已保存到:\n');
fprintf('  - optimization_results.mat (MATLAB数据)\n');
fprintf('  - optimization_results.txt (文本报告)\n');

fprintf('\n============================================\n');
fprintf('优化完成!\n');
fprintf('============================================\n');

%% ==================== GA输出函数 ====================
function [state, options, optchanged] = ga_output_func(options, state, flag, run_idx)
    persistent history
    optchanged = false;

    switch flag
        case 'init'
            history = struct('generation', [], 'best_cost', []);
        case 'iter'
            history.generation(end+1) = state.Generation;
            history.best_cost(end+1) = state.Best(end);
        case 'done'
            % 保存历史到文件
            ga_history = history;
            save(sprintf('ga_history_run%d.mat', run_idx), 'ga_history');
    end
end
