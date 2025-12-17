%% 槽轮机构七次多项式优化 - 运行脚本
% 七次多项式比五次多项式多2个自由度（加加速度边界条件）
% 可以更灵活地控制加速度曲线

clc; clear; close all;
global ga_history_temp;

fprintf('================================================\n');
fprintf('  槽轮机构七次多项式遗传算法优化\n');
fprintf('================================================\n\n');

%% 环境检测
gpu_available = false;
try
    g = gpuDevice;
    gpu_available = true;
    fprintf('[环境] GPU可用: %s\n', g.Name);
catch
    fprintf('[环境] 使用CPU模式\n');
end

%% 参数设置 - 6个优化变量
% [Angle_G_rel, Input_Stroke, v_E, a_E, j_0, j_1]
lb = [2,  25, -80, 1,  -100, -100];  % 下界
ub = [15, 60, -5,  50,  100,  100];  % 上界

% GA参数
PopSize = 80;        % 增大种群规模（因为变量增多）
MaxGen = 150;        % 增加迭代次数
CrossFrac = 0.9;
MutRate = 0.05;
NumRuns = 4;

fprintf('\n七次多项式优化变量 (6个):\n');
fprintf('  1. Angle_G_rel  : G点偏转角 [%.0f, %.0f]°\n', lb(1), ub(1));
fprintf('  2. Input_Stroke : 输入行程 [%.0f, %.0f]°\n', lb(2), ub(2));
fprintf('  3. v_E          : 终点速度 [%.0f, %.0f]\n', lb(3), ub(3));
fprintf('  4. a_E          : 终点加速度 [%.0f, %.0f]\n', lb(4), ub(4));
fprintf('  5. j_0          : 起点加加速度 [%.0f, %.0f] (新增)\n', lb(5), ub(5));
fprintf('  6. j_1          : 终点加加速度 [%.0f, %.0f] (新增)\n', lb(6), ub(6));

fprintf('\nGA参数: 种群=%d, 迭代=%d\n', PopSize, MaxGen);
fprintf('优化运行次数: %d\n\n', NumRuns);

%% 目标函数
fitness = @(x) geneva_objective_7poly(x, false);

%% GA选项
options = optimoptions('ga', ...
    'PopulationSize', PopSize, ...
    'MaxGenerations', MaxGen, ...
    'CrossoverFraction', CrossFrac, ...
    'EliteCount', 3, ...
    'FunctionTolerance', 1e-6, ...
    'Display', 'iter', ...
    'OutputFcn', @ga_output_callback);

%% 多次优化
all_params = zeros(NumRuns, 6);
all_costs = zeros(NumRuns, 1);
all_results = cell(NumRuns, 1);
all_history = cell(NumRuns, 1);

for r = 1:NumRuns
    fprintf('>>> 第 %d/%d 次优化 <<<\n', r, NumRuns);

    ga_history_temp = struct('gen', [], 'cost', []);

    tic;
    [x_opt, f_opt] = ga(fitness, 6, [], [], [], [], lb, ub, [], options);
    t_elapsed = toc;

    all_params(r,:) = x_opt;
    all_costs(r) = f_opt;
    all_history{r} = ga_history_temp;

    [~, res] = geneva_objective_7poly(x_opt, true);
    all_results{r} = res;

    fprintf('完成! 耗时%.1fs, 目标值=%.2f\n', t_elapsed, f_opt);
    fprintf('参数: Angle_G=%.2f, Stroke=%.2f, v_E=%.2f, a_E=%.2f\n', ...
        x_opt(1), x_opt(2), x_opt(3), x_opt(4));
    fprintf('      j_0=%.2f, j_1=%.2f\n\n', x_opt(5), x_opt(6));
end

%% 找最优解
[best_cost, best_idx] = min(all_costs);
best_params = all_params(best_idx, :);
best_result = all_results{best_idx};

%% 打印结果
fprintf('\n================================================\n');
fprintf('           七次多项式优化结果\n');
fprintf('================================================\n');
fprintf('%-15s', '优化参数');
for i = 1:NumRuns
    fprintf('第%d次\t\t', i);
end
fprintf('\n------------------------------------------------\n');

names = {'Angle_G_rel', 'Input_Stroke', 'v_E', 'a_E', 'j_0', 'j_1'};
for p = 1:6
    fprintf('%-15s', names{p});
    for i = 1:NumRuns
        fprintf('%.2f\t\t', all_params(i,p));
    end
    fprintf('\n');
end
fprintf('------------------------------------------------\n');
fprintf('%-15s', '目标函数值');
for i = 1:NumRuns
    fprintf('%.2f\t', all_costs(i));
end
fprintf('\n');

fprintf('\n全局最优 (第%d次): 目标值=%.4f\n', best_idx, best_cost);
fprintf('  Angle_G_rel  = %.2f°\n', best_params(1));
fprintf('  Input_Stroke = %.2f°\n', best_params(2));
fprintf('  v_E          = %.2f\n', best_params(3));
fprintf('  a_E          = %.2f\n', best_params(4));
fprintf('  j_0          = %.2f (起点加加速度)\n', best_params(5));
fprintf('  j_1          = %.2f (终点加加速度)\n', best_params(6));

%% 绘图1: 优化后运动曲线
fig1 = figure('Color','w','Position',[100 50 800 900],'Name','七次多项式优化后运动曲线');

subplot(3,1,1);
plot(best_result.Plot_X, best_result.Phi_Plot, 'k-', 'LineWidth', 2);
ylabel('角位移 (deg)'); title('槽轮角位移 (七次多项式)'); grid on;
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

saveas(fig1, 'optimized_motion_7poly.png');
fprintf('\n运动曲线已保存: optimized_motion_7poly.png\n');

%% 绘图2: GA收敛曲线
fig2 = figure('Color','w','Position',[150 100 1000 600],'Name','GA收敛曲线(七次)');

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

saveas(fig2, 'ga_convergence_7poly.png');
fprintf('收敛曲线已保存: ga_convergence_7poly.png\n');

%% 绘图3: 五次 vs 七次多项式对比
fprintf('\n正在生成五次/七次多项式对比图...\n');

% 加载五次多项式结果
if exist('ga_results.mat', 'file')
    load('ga_results.mat', 'best_result');
    result_5poly = best_result;

    fig3 = figure('Color','w','Position',[200 100 1000 800],'Name','五次vs七次多项式对比');

    % 获取七次多项式结果
    result_7poly = all_results{best_idx};

    subplot(3,1,1);
    plot(result_5poly.Plot_X, result_5poly.Phi_Plot, 'b--', 'LineWidth', 1.5); hold on;
    plot(result_7poly.Plot_X, result_7poly.Phi_Plot, 'r-', 'LineWidth', 2); hold off;
    ylabel('角位移 (deg)'); title('槽轮角位移对比 (蓝:五次, 红:七次)');
    grid on;
    xlim([0 max(result_7poly.Plot_X)+5]);

    subplot(3,1,2);
    plot(result_5poly.Plot_X, result_5poly.Vel_S, 'b--', 'LineWidth', 1.5); hold on;
    plot(result_7poly.Plot_X, result_7poly.Vel_S, 'r-', 'LineWidth', 2); hold off;
    ylabel('角速度 (rad/s)'); title('槽轮角速度对比 (蓝:五次, 红:七次)');
    grid on;
    xlim([0 max(result_7poly.Plot_X)+5]);

    subplot(3,1,3);
    plot(result_5poly.Plot_X, result_5poly.Acc_S, 'b--', 'LineWidth', 1.5); hold on;
    plot(result_7poly.Plot_X, result_7poly.Acc_S, 'r-', 'LineWidth', 2); hold off;
    ylabel('角加速度 (rad/s^2)'); xlabel('圆销转角(°)');
    title('槽轮角加速度对比 (蓝:五次, 红:七次)');
    grid on;
    xlim([0 max(result_7poly.Plot_X)+5]);

    saveas(fig3, 'poly5_vs_poly7_compare.png');
    fprintf('对比图已保存: poly5_vs_poly7_compare.png\n');
end

%% 保存结果
save('ga_results_7poly.mat', 'all_params', 'all_costs', 'best_params', 'best_cost', ...
    'all_history', 'all_results', 'best_result');
fprintf('\n结果已保存: ga_results_7poly.mat\n');

fprintf('\n================================================\n');
fprintf('七次多项式优化完成!\n');
fprintf('================================================\n');

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
    end
end
