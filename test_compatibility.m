%% 测试代码兼容性
fprintf('========================================\n');
fprintf('  测试代码兼容性\n');
fprintf('========================================\n\n');

try
    %% 1. 测试目标函数
    fprintf('[1] 测试目标函数... ');
    test_params = [10, 44.89, -20.04, 1.03];
    [cost, result] = geneva_objective(test_params, true);
    fprintf('OK (cost=%.2f)\n', cost);

    %% 2. 测试绘图（无legend）
    fprintf('[2] 测试绘图代码... ');
    fig1 = figure('Visible','off');

    subplot(3,1,1);
    plot(result.Plot_X, result.Phi_Plot, 'b--', 'LineWidth', 1.5); hold on;
    plot(result.Plot_X, result.Phi_Plot*1.05, 'r-', 'LineWidth', 2); hold off;
    ylabel('角位移 (deg)');
    title('槽轮角位移对比 (蓝:优化前, 红:优化后)');
    grid on;
    xlim([0 max(result.Plot_X)+5]);

    subplot(3,1,2);
    plot(result.Plot_X, result.Vel_S, 'b--', 'LineWidth', 1.5); hold on;
    plot(result.Plot_X, result.Vel_S*1.05, 'r-', 'LineWidth', 2); hold off;
    ylabel('角速度 (rad/s)');
    title('槽轮角速度对比 (蓝:优化前, 红:优化后)');
    grid on;
    xlim([0 max(result.Plot_X)+5]);

    subplot(3,1,3);
    plot(result.Plot_X, result.Acc_S, 'b--', 'LineWidth', 1.5); hold on;
    plot(result.Plot_X, result.Acc_S*1.05, 'r-', 'LineWidth', 2); hold off;
    ylabel('角加速度 (rad/s^2)'); xlabel('圆销转角(°)');
    title('槽轮角加速度对比 (蓝:优化前, 红:优化后)');
    grid on;
    xlim([0 max(result.Plot_X)+5]);

    saveas(fig1, 'test_plot.png');
    close(fig1);
    fprintf('OK\n');

    %% 3. 测试GA（简化版，只跑5代）
    fprintf('[3] 测试GA优化（5代）... ');
    fitness = @(x) geneva_objective(x, false);
    lb = [2, 25, -60, 1];
    ub = [10, 45, -5, 30];
    options = optimoptions('ga', ...
        'PopulationSize', 10, ...
        'MaxGenerations', 5, ...
        'Display', 'off');
    [x_opt, f_opt] = ga(fitness, 4, [], [], [], [], lb, ub, [], options);
    fprintf('OK (最优=%.2f)\n', f_opt);

    %% 4. 测试输出回调（简化版）
    fprintf('[4] 测试GA回调函数... ');
    global ga_history_temp;
    ga_history_temp = struct('gen', [], 'cost', []);

    options2 = optimoptions('ga', ...
        'PopulationSize', 10, ...
        'MaxGenerations', 3, ...
        'Display', 'off', ...
        'OutputFcn', @test_callback);
    [x2, f2] = ga(fitness, 4, [], [], [], [], lb, ub, [], options2);
    fprintf('OK (记录了%d代)\n', length(ga_history_temp.gen));

    fprintf('\n========================================\n');
    fprintf('  所有测试通过!\n');
    fprintf('========================================\n');

catch ME
    fprintf('\n\n!!! 错误 !!!\n');
    fprintf('消息: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('位置: %s 第%d行\n', ME.stack(1).name, ME.stack(1).line);
    end
end

%% 回调函数
function [state, options, optchanged] = test_callback(options, state, flag)
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
