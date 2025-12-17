%% 槽轮机构优化 - 快速运行脚本
% 用法: 在MATLAB中直接运行此脚本
%
% 功能:
%   1. 自动检测GPU/CPU环境
%   2. 运行遗传算法优化
%   3. 生成运动曲线和收敛曲线
%   4. 保存优化结果

clc; clear; close all;
global ga_history_temp;

fprintf('========================================\n');
fprintf('  槽轮机构遗传算法优化\n');
fprintf('========================================\n\n');

%% 环境检测
gpu_available = false;
try
    g = gpuDevice;
    gpu_available = true;
    fprintf('[环境] GPU可用: %s (%.1f GB)\n', g.Name, g.TotalMemory/1e9);
catch
    fprintf('[环境] 使用CPU模式\n');
end

%% 参数设置
% 优化变量边界: [Angle_G_rel, Input_Stroke, v_E, a_E]
lb = [2,  25, -60, 1];
ub = [10, 45, -5,  30];

% GA参数 (表8)
PopSize = 50;        % 种群规模
MaxGen = 100;        % 迭代次数
CrossFrac = 0.9;     % 交叉概率
MutRate = 0.05;      % 变异概率
NumRuns = 4;         % 运行次数

fprintf('\nGA参数: 种群=%d, 迭代=%d, 交叉=%.2f, 变异=%.2f\n', ...
    PopSize, MaxGen, CrossFrac, MutRate);
fprintf('优化运行次数: %d\n\n', NumRuns);

%% 目标函数
if gpu_available
    fitness = @(x) geneva_objective_gpu(x, false, true);
else
    fitness = @(x) geneva_objective(x, false);
end

%% GA选项
options = optimoptions('ga', ...
    'PopulationSize', PopSize, ...
    'MaxGenerations', MaxGen, ...
    'CrossoverFraction', CrossFrac, ...
    'EliteCount', 2, ...
    'FunctionTolerance', 1e-6, ...
    'Display', 'iter', ...
    'OutputFcn', @ga_output_callback);

% 尝试启用并行计算
try
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    options = optimoptions(options, 'UseParallel', true);
    fprintf('[并行] 已启用并行计算\n\n');
catch
    fprintf('[串行] 使用串行模式\n\n');
end

%% 多次优化
all_params = zeros(NumRuns, 4);
all_costs = zeros(NumRuns, 1);
all_results = cell(NumRuns, 1);
all_history = cell(NumRuns, 1);

for r = 1:NumRuns
    fprintf('>>> 第 %d/%d 次优化 <<<\n', r, NumRuns);

    % 清空历史记录
    ga_history_temp = struct('gen', [], 'cost', []);

    tic;
    [x_opt, f_opt] = ga(fitness, 4, [], [], [], [], lb, ub, [], options);
    t_elapsed = toc;

    all_params(r,:) = x_opt;
    all_costs(r) = f_opt;
    all_history{r} = ga_history_temp;

    % 获取详细结果
    if gpu_available
        [~, res] = geneva_objective_gpu(x_opt, true, true);
    else
        [~, res] = geneva_objective(x_opt, true);
    end
    all_results{r} = res;

    fprintf('完成! 耗时%.1fs, 目标值=%.2f\n', t_elapsed, f_opt);
    fprintf('参数: Angle_G=%.2f, Stroke=%.2f, v_E=%.2f, a_E=%.2f\n\n', ...
        x_opt(1), x_opt(2), x_opt(3), x_opt(4));
end

%% 找最优解
[best_cost, best_idx] = min(all_costs);
best_params = all_params(best_idx, :);
best_result = all_results{best_idx};

%% 打印结果表格 (表9格式)
fprintf('\n========================================\n');
fprintf('           优化结果 (表9)\n');
fprintf('========================================\n');
fprintf('%-15s', '优化参数');
for i = 1:NumRuns
    fprintf('第%d次\t', i);
end
fprintf('\n----------------------------------------\n');

names = {'Angle_G_rel', 'Input_Stroke', 'v_E', 'a_E'};
for p = 1:4
    fprintf('%-15s', names{p});
    for i = 1:NumRuns
        fprintf('%.2f\t', all_params(i,p));
    end
    fprintf('\n');
end
fprintf('----------------------------------------\n');
fprintf('%-15s', '目标函数值');
for i = 1:NumRuns
    fprintf('%.2f\t', all_costs(i));
end
fprintf('\n');

fprintf('\n全局最优 (第%d次): 目标值=%.4f\n', best_idx, best_cost);
fprintf('  Angle_G_rel  = %.2f\n', best_params(1));
fprintf('  Input_Stroke = %.2f\n', best_params(2));
fprintf('  v_E          = %.2f\n', best_params(3));
fprintf('  a_E          = %.2f\n', best_params(4));

%% 绘图1: 优化后运动曲线
fig1 = figure('Color','w','Position',[100 50 800 900],'Name','优化后运动曲线');

subplot(3,1,1);
plot(best_result.Plot_X, best_result.Phi_Plot, 'k-', 'LineWidth', 2);
ylabel('角位移 (deg)'); title('槽轮角位移'); grid on;
xlim([0 max(best_result.Plot_X)+5]);

subplot(3,1,2);
plot(best_result.Plot_X, best_result.Vel_S, 'k-', 'LineWidth', 2);
ylabel('角速度 (rad/s)'); title('槽轮角速度'); grid on;
xlim([0 max(best_result.Plot_X)+5]);

subplot(3,1,3);
plot(best_result.Plot_X, best_result.Acc_S, 'k-', 'LineWidth', 2);
ylabel('角加速度 (rad/s^2)'); xlabel('圆销转角(°)');
title('槽轮角加速度'); grid on;
xlim([0 max(best_result.Plot_X)+5]);

saveas(fig1, 'optimized_motion.png');
fprintf('\n运动曲线已保存: optimized_motion.png\n');

%% 绘图2: GA收敛曲线
fig2 = figure('Color','w','Position',[150 100 1000 600],'Name','GA收敛曲线');

for i = 1:NumRuns
    subplot(2,2,i);
    h = all_history{i};
    if ~isempty(h.gen)
        plot(h.gen, h.cost, 'k-', 'LineWidth', 1.5);
    end
    xlabel('迭代次数'); ylabel('最优解');
    title(sprintf('(%c) 第%d次', 'a'+i-1, i));
    grid on; xlim([0 MaxGen]);
end

saveas(fig2, 'ga_convergence.png');
fprintf('收敛曲线已保存: ga_convergence.png\n');

%% 绘图3: 优化前后对比
fprintf('\n正在生成对比图...\n');

% 原始参数
orig_params = [7.5, 40, -30, 9];
if gpu_available
    [~, orig_result] = geneva_objective_gpu(orig_params, true, true);
else
    [~, orig_result] = geneva_objective(orig_params, true);
end

fig3 = figure('Color','w','Position',[200 100 1000 800],'Name','优化前后对比');

subplot(3,1,1);
plot(orig_result.Plot_X, orig_result.Phi_Plot, 'b--', 'LineWidth', 1.5); hold on;
plot(best_result.Plot_X, best_result.Phi_Plot, 'r-', 'LineWidth', 2);
ylabel('角位移 (deg)'); title('槽轮角位移对比');
legend('优化前','优化后'); grid on;
xlim([0 max(best_result.Plot_X)+5]);

subplot(3,1,2);
plot(orig_result.Plot_X, orig_result.Vel_S, 'b--', 'LineWidth', 1.5); hold on;
plot(best_result.Plot_X, best_result.Vel_S, 'r-', 'LineWidth', 2);
ylabel('角速度 (rad/s)'); title('槽轮角速度对比');
legend('优化前','优化后'); grid on;
xlim([0 max(best_result.Plot_X)+5]);

subplot(3,1,3);
plot(orig_result.Plot_X, orig_result.Acc_S, 'b--', 'LineWidth', 1.5); hold on;
plot(best_result.Plot_X, best_result.Acc_S, 'r-', 'LineWidth', 2);
ylabel('角加速度 (rad/s^2)'); xlabel('圆销转角(°)');
title('槽轮角加速度对比');
legend('优化前','优化后'); grid on;
xlim([0 max(best_result.Plot_X)+5]);

saveas(fig3, 'optimization_compare.png');
fprintf('对比图已保存: optimization_compare.png\n');

%% 保存结果
save('ga_results.mat', 'all_params', 'all_costs', 'best_params', 'best_cost', ...
    'all_history', 'all_results', 'orig_result', 'best_result');
fprintf('\n结果已保存: ga_results.mat\n');

fprintf('\n========================================\n');
fprintf('优化完成!\n');
fprintf('========================================\n');

%% GA输出回调函数
function [state, options, optchanged] = ga_output_callback(options, state, flag)
    global ga_history_temp;
    optchanged = false;

    switch flag
        case 'init'
            ga_history_temp = struct('gen', [], 'cost', []);
        case 'iter'
            ga_history_temp.gen(end+1) = state.Generation;
            ga_history_temp.cost(end+1) = state.Best(end);
        case 'done'
            % 保存完成
    end
end
