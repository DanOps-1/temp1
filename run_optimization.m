%% 槽轮机构遗传算法优化
% 直接在MATLAB里运行这个脚本就行
% 会自动检测有没有GPU，有就用GPU加速

clc; clear; close all;
global ga_history_temp;

fprintf('========================================\n');
fprintf('  槽轮机构遗传算法优化\n');
fprintf('========================================\n\n');

% 先看看有没有GPU可以用
gpu_available = false;
try
    g = gpuDevice;
    gpu_available = true;
    fprintf('[环境] 检测到GPU: %s (显存%.1fGB)\n', g.Name, g.TotalMemory/1e9);
catch
    fprintf('[环境] 没检测到GPU，用CPU跑\n');
end

% 设置参数范围
% 4个优化变量: G点偏转角、输入行程、速度边界、加速度边界
lb = [2,  25, -60, 1];    % 下限
ub = [10, 45, -5,  30];   % 上限

% 遗传算法的参数，按论文表8设的
PopSize = 50;        % 种群大小，就是每一代有多少个个体
MaxGen = 100;        % 最多跑100代
CrossFrac = 0.9;     % 交叉概率90%
MutRate = 0.05;      % 变异概率5%
NumRuns = 4;         % 跑4次取最好的

fprintf('\nGA参数: 种群=%d, 迭代=%d, 交叉=%.2f, 变异=%.2f\n', ...
    PopSize, MaxGen, CrossFrac, MutRate);
fprintf('总共跑%d次\n\n', NumRuns);

% 选目标函数，有GPU就用GPU版本
if gpu_available
    fitness = @(x) geneva_objective_gpu(x, false, true);
else
    fitness = @(x) geneva_objective(x, false);
end

% 配置GA的选项
options = optimoptions('ga', ...
    'PopulationSize', PopSize, ...
    'MaxGenerations', MaxGen, ...
    'CrossoverFraction', CrossFrac, ...
    'EliteCount', 2, ...           % 每代保留2个最优个体
    'FunctionTolerance', 1e-6, ... % 收敛精度
    'Display', 'iter', ...         % 显示迭代过程
    'OutputFcn', @ga_output_callback);

% 试试能不能开并行，能开就开
try
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    options = optimoptions(options, 'UseParallel', true);
    fprintf('[并行] 并行计算已开启\n\n');
catch
    fprintf('[串行] 并行开不了，单线程跑\n\n');
end

% 准备存结果的变量
all_params = zeros(NumRuns, 4);
all_costs = zeros(NumRuns, 1);
all_results = cell(NumRuns, 1);
all_history = cell(NumRuns, 1);

% 开始优化，跑NumRuns次
for r = 1:NumRuns
    fprintf('>>> 第%d/%d次优化 <<<\n', r, NumRuns);

    % 清空历史记录
    ga_history_temp = struct('gen', [], 'cost', []);

    % 跑GA
    tic;
    [x_opt, f_opt] = ga(fitness, 4, [], [], [], [], lb, ub, [], options);
    t_elapsed = toc;

    % 存结果
    all_params(r,:) = x_opt;
    all_costs(r) = f_opt;
    all_history{r} = ga_history_temp;

    % 拿详细数据，画图要用
    if gpu_available
        [~, res] = geneva_objective_gpu(x_opt, true, true);
    else
        [~, res] = geneva_objective(x_opt, true);
    end
    all_results{r} = res;

    fprintf('完成! 用时%.1f秒, 目标值=%.2f\n', t_elapsed, f_opt);
    fprintf('参数: Angle_G=%.2f, Stroke=%.2f, v_E=%.2f, a_E=%.2f\n\n', ...
        x_opt(1), x_opt(2), x_opt(3), x_opt(4));
end

% 找出最好的那一次
[best_cost, best_idx] = min(all_costs);
best_params = all_params(best_idx, :);
best_result = all_results{best_idx};

% 打印结果表格
fprintf('\n========================================\n');
fprintf('           优化结果汇总\n');
fprintf('========================================\n');
fprintf('%-15s', '参数');
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

fprintf('\n最优解是第%d次: 目标值=%.4f\n', best_idx, best_cost);
fprintf('  Angle_G_rel  = %.2f\n', best_params(1));
fprintf('  Input_Stroke = %.2f\n', best_params(2));
fprintf('  v_E          = %.2f\n', best_params(3));
fprintf('  a_E          = %.2f\n', best_params(4));

% ========== 画图 ==========

% 图1: 优化后的运动曲线
fig1 = figure('Color','w','Position',[100 50 900 1000],'Name','优化后运动曲线');

% 字体设置
titleFontSize = 14;
labelFontSize = 12;
tickFontSize = 10;

ax1 = subplot(3,1,1);
plot(best_result.Plot_X, best_result.Phi_Plot, 'k-', 'LineWidth', 2);
ylabel('角位移 (°)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title('(a) 槽轮角位移', 'FontSize', titleFontSize, 'FontWeight', 'bold');
grid on; xlim([0 max(best_result.Plot_X)+5]);
set(ax1, 'FontSize', tickFontSize, 'LineWidth', 1);
set(ax1, 'Position', [0.15 0.71 0.80 0.24]);

ax2 = subplot(3,1,2);
plot(best_result.Plot_X, best_result.Vel_S, 'k-', 'LineWidth', 2);
ylabel('角速度 (rad/s)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title('(b) 槽轮角速度', 'FontSize', titleFontSize, 'FontWeight', 'bold');
grid on; xlim([0 max(best_result.Plot_X)+5]);
set(ax2, 'FontSize', tickFontSize, 'LineWidth', 1);
set(ax2, 'Position', [0.15 0.40 0.80 0.24]);

ax3 = subplot(3,1,3);
plot(best_result.Plot_X, best_result.Acc_S, 'k-', 'LineWidth', 2);
ylabel('角加速度 (rad/s²)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
xlabel('圆销转角 (°)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title('(c) 槽轮角加速度', 'FontSize', titleFontSize, 'FontWeight', 'bold');
grid on; xlim([0 max(best_result.Plot_X)+5]);
set(ax3, 'FontSize', tickFontSize, 'LineWidth', 1);
set(ax3, 'Position', [0.15 0.08 0.80 0.24]);

saveas(fig1, 'optimized_motion.png');
fprintf('\n运动曲线已保存: optimized_motion.png\n');

% 图2: GA收敛曲线
fig2 = figure('Color','w','Position',[150 100 1100 700],'Name','GA收敛曲线');

positions = [0.08 0.58 0.38 0.35;   % (a)
             0.55 0.58 0.38 0.35;   % (b)
             0.08 0.12 0.38 0.35;   % (c)
             0.55 0.12 0.38 0.35];  % (d)

for i = 1:NumRuns
    ax = subplot(2,2,i);
    h = all_history{i};
    if ~isempty(h.gen)
        plot(h.gen, h.cost/1e5, 'k-', 'LineWidth', 1.5);
    end
    xlabel('迭代次数', 'FontSize', labelFontSize, 'FontWeight', 'bold');
    ylabel('目标函数值 (×10⁵)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
    title(sprintf('(%c) 第%d次优化', 'a'+i-1, i), 'FontSize', titleFontSize, 'FontWeight', 'bold');
    grid on; xlim([0 MaxGen]);
    set(ax, 'FontSize', tickFontSize, 'LineWidth', 1);
    set(ax, 'Position', positions(i,:));
end

saveas(fig2, 'ga_convergence.png');
fprintf('收敛曲线已保存: ga_convergence.png\n');

% 图3: 优化前后对比
fprintf('\n正在生成对比图...\n');

% 原始参数(my.m里的默认值)
orig_params = [7.5, 40, -30, 9];
if gpu_available
    [~, orig_result] = geneva_objective_gpu(orig_params, true, true);
else
    [~, orig_result] = geneva_objective(orig_params, true);
end

fig3 = figure('Color','w','Position',[200 100 900 1000],'Name','优化前后对比');

ax1 = subplot(3,1,1);
plot(orig_result.Plot_X, orig_result.Phi_Plot, 'b--', 'LineWidth', 1.5); hold on;
plot(best_result.Plot_X, best_result.Phi_Plot, 'r-', 'LineWidth', 2); hold off;
ylabel('角位移 (°)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title('(a) 槽轮角位移对比', 'FontSize', titleFontSize, 'FontWeight', 'bold');
grid on; xlim([0 max(best_result.Plot_X)+5]);
text(10, max(best_result.Phi_Plot)*0.85, '— 优化后', 'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');
text(10, max(best_result.Phi_Plot)*0.65, '-- 优化前', 'Color', 'b', 'FontSize', 11, 'FontWeight', 'bold');
set(ax1, 'FontSize', tickFontSize, 'LineWidth', 1);
set(ax1, 'Position', [0.15 0.71 0.80 0.24]);

ax2 = subplot(3,1,2);
plot(orig_result.Plot_X, orig_result.Vel_S, 'b--', 'LineWidth', 1.5); hold on;
plot(best_result.Plot_X, best_result.Vel_S, 'r-', 'LineWidth', 2); hold off;
ylabel('角速度 (rad/s)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title('(b) 槽轮角速度对比', 'FontSize', titleFontSize, 'FontWeight', 'bold');
grid on; xlim([0 max(best_result.Plot_X)+5]);
text(10, max(best_result.Vel_S)*0.85, '— 优化后', 'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');
text(10, max(best_result.Vel_S)*0.65, '-- 优化前', 'Color', 'b', 'FontSize', 11, 'FontWeight', 'bold');
set(ax2, 'FontSize', tickFontSize, 'LineWidth', 1);
set(ax2, 'Position', [0.15 0.40 0.80 0.24]);

ax3 = subplot(3,1,3);
plot(orig_result.Plot_X, orig_result.Acc_S, 'b--', 'LineWidth', 1.5); hold on;
plot(best_result.Plot_X, best_result.Acc_S, 'r-', 'LineWidth', 2); hold off;
ylabel('角加速度 (rad/s²)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
xlabel('圆销转角 (°)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title('(c) 槽轮角加速度对比', 'FontSize', titleFontSize, 'FontWeight', 'bold');
grid on; xlim([0 max(best_result.Plot_X)+5]);
ymax = max(abs([orig_result.Acc_S; best_result.Acc_S]));
text(10, ymax*0.7, '— 优化后', 'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');
text(10, ymax*0.4, '-- 优化前', 'Color', 'b', 'FontSize', 11, 'FontWeight', 'bold');
set(ax3, 'FontSize', tickFontSize, 'LineWidth', 1);
set(ax3, 'Position', [0.15 0.08 0.80 0.24]);

saveas(fig3, 'optimization_compare.png');
fprintf('对比图已保存: optimization_compare.png\n');

% 保存数据
save('ga_results.mat', 'all_params', 'all_costs', 'best_params', 'best_cost', ...
    'all_history', 'all_results', 'orig_result', 'best_result');
fprintf('\n数据已保存: ga_results.mat\n');

fprintf('\n========================================\n');
fprintf('搞定!\n');
fprintf('========================================\n');

% GA的回调函数，用来记录每一代的最优值
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
            % 结束了，不用干啥
    end
end
