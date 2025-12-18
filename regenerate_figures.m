%% 重新生成优化后的图片（使用已保存的数据）
% 用于测试新的绘图格式

clc; close all;

% 加载数据
load('ga_results.mat');

% 字体设置
titleFontSize = 14;
labelFontSize = 12;
tickFontSize = 10;
MaxGen = 100;

fprintf('重新生成图片...\n');

%% 图1: 优化后的运动曲线
fig1 = figure('Color','w','Position',[100 50 900 1000],'Name','优化后运动曲线');

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
fprintf('图1已保存: optimized_motion.png\n');

%% 图2: GA收敛曲线
fig2 = figure('Color','w','Position',[150 100 1100 700],'Name','GA收敛曲线');

positions = [0.08 0.58 0.38 0.35;
             0.55 0.58 0.38 0.35;
             0.08 0.12 0.38 0.35;
             0.55 0.12 0.38 0.35];

NumRuns = length(all_history);
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
fprintf('图2已保存: ga_convergence.png\n');

%% 图3: 优化前后对比
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
fprintf('图3已保存: optimization_compare.png\n');

fprintf('\n全部图片生成完成!\n');
