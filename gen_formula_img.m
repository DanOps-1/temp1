% 生成目标函数公式图片
fig = figure('Color','w','Position',[100 100 700 550]);
axis off;

% 目标函数
text(0.02, 0.95, '目标函数：', 'FontSize', 14, 'FontWeight', 'bold');

text(0.12, 0.85, '$J = J_1 + J_2 + J_3 + J_4$', 'FontSize', 14, 'Interpreter', 'latex');

text(0.12, 0.74, '$J_1 = 1000\sum|\Delta\phi^-| + 500 \cdot N_{neg}$', 'FontSize', 12, 'Interpreter', 'latex');
text(0.75, 0.74, '(1)', 'FontSize', 12);

text(0.12, 0.64, '$J_2 = 0.001 \cdot \Delta t \sum \dot{\omega}^2 + 10 \cdot |\bar{\omega}|_{108^\circ}$', 'FontSize', 12, 'Interpreter', 'latex');
text(0.75, 0.64, '(2)', 'FontSize', 12);

text(0.12, 0.54, '$J_3 = 10^{-6} \cdot \Delta t \sum \dot{\alpha}^2 + 0.1 \cdot (\alpha_{max}-500)$', 'FontSize', 12, 'Interpreter', 'latex');
text(0.75, 0.54, '(3)', 'FontSize', 12);

text(0.12, 0.44, '$J_4 = 0.01\sum \kappa^2 + 5(|\omega_0|+|\omega_T|)$', 'FontSize', 12, 'Interpreter', 'latex');
text(0.75, 0.44, '(4)', 'FontSize', 12);

% 设计变量
text(0.02, 0.32, '设计变量：', 'FontSize', 14, 'FontWeight', 'bold');
text(0.22, 0.32, '$\theta_G$', 'FontSize', 13, 'Interpreter', 'latex');
text(0.30, 0.32, '、', 'FontSize', 13);
text(0.33, 0.32, '$\theta_S$', 'FontSize', 13, 'Interpreter', 'latex');
text(0.41, 0.32, '、', 'FontSize', 13);
text(0.44, 0.32, '$v_E$', 'FontSize', 13, 'Interpreter', 'latex');
text(0.51, 0.32, '、', 'FontSize', 13);
text(0.54, 0.32, '$a_E$', 'FontSize', 13, 'Interpreter', 'latex');

% 约束条件
text(0.02, 0.20, '约束条件：', 'FontSize', 14, 'FontWeight', 'bold');

text(0.15, 0.12, '$\left\{ \begin{array}{l} 2^\circ \leq \theta_G \leq 10^\circ \\ 25^\circ \leq \theta_S \leq 45^\circ \\ -60 \leq v_E \leq -5 \\ 1 \leq a_E \leq 30 \end{array} \right.$', ...
    'FontSize', 12, 'Interpreter', 'latex');
text(0.75, 0.08, '(5)', 'FontSize', 12);

saveas(fig, 'objective_function.png');
fprintf('已保存: objective_function.png\n');
close(fig);
