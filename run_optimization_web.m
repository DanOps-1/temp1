%% 槽轮机构遗传算法优化 - 网页版专用
% 去掉了GPU、并行计算等网页版不支持的功能

clc; clear; close all;
global ga_history_temp;

fprintf('========================================\n');
fprintf('  槽轮机构遗传算法优化 (网页版)\n');
fprintf('========================================\n\n');

% 设置参数范围
lb = [2,  25, -60, 1];    % 下限
ub = [10, 45, -5,  30];   % 上限

% GA参数（和GPU版一致）
PopSize = 50;        % 种群大小
MaxGen = 100;        % 迭代次数
CrossFrac = 0.9;
MutRate = 0.05;
NumRuns = 4;         % 跑4次取最好的

fprintf('GA参数: 种群=%d, 迭代=%d\n', PopSize, MaxGen);
fprintf('总共跑%d次\n\n', NumRuns);

% 目标函数
fitness = @(x) geneva_objective(x, false);

% GA选项（不用并行）
options = optimoptions('ga', ...
    'PopulationSize', PopSize, ...
    'MaxGenerations', MaxGen, ...
    'CrossoverFraction', CrossFrac, ...
    'EliteCount', 2, ...
    'FunctionTolerance', 1e-6, ...
    'Display', 'iter', ...
    'UseParallel', false, ...
    'OutputFcn', @ga_output_callback);

% 存结果
all_params = zeros(NumRuns, 4);
all_costs = zeros(NumRuns, 1);
all_results = cell(NumRuns, 1);
all_history = cell(NumRuns, 1);

% 开始优化
for r = 1:NumRuns
    fprintf('>>> 第%d/%d次优化 <<<\n', r, NumRuns);

    ga_history_temp = struct('gen', [], 'cost', []);

    tic;
    [x_opt, f_opt] = ga(fitness, 4, [], [], [], [], lb, ub, [], options);
    t_elapsed = toc;

    all_params(r,:) = x_opt;
    all_costs(r) = f_opt;
    all_history{r} = ga_history_temp;

    [~, res] = geneva_objective(x_opt, true);
    all_results{r} = res;

    fprintf('完成! 用时%.1f秒, 目标值=%.2f\n', t_elapsed, f_opt);
    fprintf('参数: Angle_G=%.2f, Stroke=%.2f, v_E=%.2f, a_E=%.2f\n\n', ...
        x_opt(1), x_opt(2), x_opt(3), x_opt(4));
end

% 找最优
[best_cost, best_idx] = min(all_costs);
best_params = all_params(best_idx, :);
best_result = all_results{best_idx};

% 打印结果
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

%% 画图
% 图1: 优化后运动曲线
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

saveas(fig1, 'optimized_motion_web.png');
fprintf('\n运动曲线已保存: optimized_motion_web.png\n');

% 图2: GA收敛曲线
fig2 = figure('Color','w','Position',[150 100 1000 600],'Name','GA收敛曲线');

for i = 1:NumRuns
    subplot(2, 2, i);
    h = all_history{i};
    if ~isempty(h.gen)
        plot(h.gen, h.cost, 'k-', 'LineWidth', 1.5);
    end
    xlabel('迭代次数'); ylabel('最优解');
    title(sprintf('(%c) 第%d次', 'a'+i-1, i));
    grid on; xlim([0 MaxGen]);
end

saveas(fig2, 'ga_convergence_web.png');
fprintf('收敛曲线已保存: ga_convergence_web.png\n');

% 图3: 优化前后对比
fprintf('\n正在生成对比图...\n');

orig_params = [7.5, 40, -30, 9];
[~, orig_result] = geneva_objective(orig_params, true);

fig3 = figure('Color','w','Position',[200 100 1000 800],'Name','优化前后对比');

subplot(3,1,1);
plot(orig_result.Plot_X, orig_result.Phi_Plot, 'b--', 'LineWidth', 1.5); hold on;
plot(best_result.Plot_X, best_result.Phi_Plot, 'r-', 'LineWidth', 2); hold off;
ylabel('角位移 (deg)'); title('槽轮角位移对比 (蓝:优化前, 红:优化后)');
grid on; xlim([0 max(best_result.Plot_X)+5]);

subplot(3,1,2);
plot(orig_result.Plot_X, orig_result.Vel_S, 'b--', 'LineWidth', 1.5); hold on;
plot(best_result.Plot_X, best_result.Vel_S, 'r-', 'LineWidth', 2); hold off;
ylabel('角速度 (rad/s)'); title('槽轮角速度对比 (蓝:优化前, 红:优化后)');
grid on; xlim([0 max(best_result.Plot_X)+5]);

subplot(3,1,3);
plot(orig_result.Plot_X, orig_result.Acc_S, 'b--', 'LineWidth', 1.5); hold on;
plot(best_result.Plot_X, best_result.Acc_S, 'r-', 'LineWidth', 2); hold off;
ylabel('角加速度 (rad/s^2)'); xlabel('圆销转角(°)');
title('槽轮角加速度对比 (蓝:优化前, 红:优化后)');
grid on; xlim([0 max(best_result.Plot_X)+5]);

saveas(fig3, 'optimization_compare_web.png');
fprintf('对比图已保存: optimization_compare_web.png\n');

% 保存数据
save('ga_results_web.mat', 'all_params', 'all_costs', 'best_params', 'best_cost', ...
    'all_history', 'all_results', 'orig_result', 'best_result');
fprintf('\n数据已保存: ga_results_web.mat\n');

fprintf('\n========================================\n');
fprintf('搞定!\n');
fprintf('========================================\n');

%% GA回调函数
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
