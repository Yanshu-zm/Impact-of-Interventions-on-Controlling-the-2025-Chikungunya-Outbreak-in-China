% 基于保存的统计数据生成图表
% 首先加载之前保存的统计数据
load('../figure/reproduction_stats_2025-09-06-21-51-59.mat');  % 替换为实际保存的.mat文件路径

% 设置默认样式
gray_color = [0.7 0.7 0.7];
line_width = 1;
ci_alpha = 0.3;  % 置信区间透明度
time_index = time_index';
%% 第一个进口数据图表：无干预 vs 提前1周 vs 提前2周
f1 = figure();
set(f1, "Position", [1000, 1007, 560, 230]);
hold on;

% 1. 无干预场景（对应stats_data.imports(1)）
plot(stats_data.imports(1).cumsum_mean, 'k', 'LineWidth', line_width);
hb = fillyy(time_index, ...
    stats_data.imports(1).cumsum_p5, ...
    stats_data.imports(1).cumsum_p95, ...
    gray_color, ci_alpha);

% 2. 提前1周场景（对应stats_data.imports(2)）
plot(stats_data.imports(2).cumsum_mean, 'blue', 'LineWidth', line_width);
hb1 = fillyy(time_index, ...
    stats_data.imports(2).cumsum_p5, ...
    stats_data.imports(2).cumsum_p95, ...
    'blue', ci_alpha);

% 3. 提前2周场景（对应stats_data.imports(3)）
plot(stats_data.imports(3).cumsum_mean, 'cyan', 'LineWidth', line_width);
hb2 = fillyy(time_index, ...
    stats_data.imports(3).cumsum_p5, ...
    stats_data.imports(3).cumsum_p95, ...
    'cyan', ci_alpha);

% 图表配置
legend([hb, hb1, hb2], {"Without Intervention", "One week early", "Two weeks early"});
xlim([9, 63]);
xticks(1:10:60);
xticklabels(datestr(tickDates, 'mm-dd'));
xtickangle(45);
ylabel("Accumulate exportion");

% 保存图表
% export_file_name1 = fullfile(directory_name, ["import_p1_", currentDateTimeString]);
% saveas(f1, [export_file_name1, ".fig"]);
% exportgraphics(f1, [export_file_name1, ".pdf"]);


%% 第二个进口数据图表：无干预 vs 强度1.5倍 vs 强度1.5倍+提前1周
f2 = figure();
set(f2, "Position", [1000, 1007, 560, 230]);
hold on;

% 1. 无干预场景（对应stats_data.imports(1)）
plot(stats_data.imports(1).cumsum_mean, 'k', 'LineWidth', line_width);
hb = fillyy(time_index, ...
    stats_data.imports(1).cumsum_p5, ...
    stats_data.imports(1).cumsum_p95, ...
    gray_color, ci_alpha);

% 2. 强度1.5倍场景（对应stats_data.imports(4)）
plot(stats_data.imports(4).cumsum_mean, 'green', 'LineWidth', line_width);
hb_s15 = fillyy(time_index, ...
    stats_data.imports(4).cumsum_p5, ...
    stats_data.imports(4).cumsum_p95, ...
    'green', ci_alpha);

% 3. 强度1.5倍+提前1周场景（对应stats_data.imports(5)）
plot(stats_data.imports(5).cumsum_mean, 'magenta', 'LineWidth', line_width);
hb_s15_owh = fillyy(time_index, ...
    stats_data.imports(5).cumsum_p5, ...
    stats_data.imports(5).cumsum_p95, ...
    'magenta', ci_alpha);

% 图表配置
legend([hb, hb_s15, hb_s15_owh], ...
    {"Without Intervention", "Strength x1.5", "Strength x1.5 and one week early"});
xlim([9, 63]);
xticks(1:10:60);
xticklabels(datestr(tickDates, 'mm-dd'));
xtickangle(45);
ylabel("Accumulate exportion");

% 保存图表
% export_file_name2 = fullfile(directory_name, ["import_p2_", currentDateTimeString]);
% saveas(f2, [export_file_name2, ".fig"]);
% exportgraphics(f2, [export_file_name2, ".pdf"]);


%% 第一个新感染数据图表：基准 vs 提前1周干预 vs 提前2周干预
f3 = figure();
set(f3, "Position", [1000, 1007, 560, 230]);
hold on;

% 1. 基准场景（对应stats_data.infections(1)）
plot(stats_data.infections(1).mean, 'k', 'LineWidth', line_width);
hb1 = fillyy(time_index, ...
    stats_data.infections(1).p5, ...
    stats_data.infections(1).p95, ...
    gray_color, ci_alpha);

% 2. 提前1周干预（对应stats_data.infections(2)）
plot(stats_data.infections(2).mean, 'blue', 'LineWidth', line_width);
hb2 = fillyy(time_index, ...
    stats_data.infections(2).p5, ...
    stats_data.infections(2).p95, ...
    'blue', ci_alpha);

% 3. 提前2周干预（对应stats_data.infections(3)）
plot(stats_data.infections(3).mean, 'cyan', 'LineWidth', line_width);
hb3 = fillyy(time_index, ...
    stats_data.infections(3).p5, ...
    stats_data.infections(3).p95, ...
    'cyan', ci_alpha);

% 图表配置
legend([hb1, hb2, hb3], ...
    {"Fitting", "Intervention one week earlier", "Intervention two week earlier"});
xticks(1:10:60);
xticklabels(datestr(tickDates, 'mm-dd'));
xtickangle(45);

% 保存图表
% export_file_name3 = fullfile(directory_name, ["week_ahead-", currentDateTimeString]);
% saveas(f3, [export_file_name3, ".fig"]);
% exportgraphics(f3, [export_file_name3, ".pdf"]);


%% 第二个新感染数据图表：基准 vs 强度1.5倍 vs 强度1.5倍+提前1周
f4 = figure();
set(f4, "Position", [1000, 1007, 560, 230]);
hold on;

% 1. 基准场景（对应stats_data.infections(1)）
plot(stats_data.infections(1).mean, 'k', 'LineWidth', line_width);
hb1 = fillyy(time_index, ...
    stats_data.infections(1).p5, ...
    stats_data.infections(1).p95, ...
    gray_color, ci_alpha);

% 2. 强度1.5倍（对应stats_data.infections(4)）
plot(stats_data.infections(4).mean, 'green', 'LineWidth', line_width);
hb2 = fillyy(time_index, ...
    stats_data.infections(4).p5, ...
    stats_data.infections(4).p95, ...
    'green', 0.1);  % 注意此处透明度与原始一致

% 3. 强度1.5倍+提前1周（对应stats_data.infections(5)）
plot(stats_data.infections(5).mean, 'green', 'LineWidth', line_width);
hb3 = fillyy(time_index, ...
    stats_data.infections(5).p5, ...
    stats_data.infections(5).p95, ...
    'green', ci_alpha);

% 图表配置
legend([hb1, hb2, hb3], ...
    {"Fitting", "1.5 stronger", "1.5 stronger + one week ahead"});
xticks(1:10:60);
xticklabels(datestr(tickDates, 'mm-dd'));
xtickangle(45);

% % 保存图表
% export_file_name4 = fullfile(directory_name, ["strength-", currentDateTimeString]);
% saveas(f4, [export_file_name4, ".fig"]);
% exportgraphics(f4, [export_file_name4, ".pdf"]);

hold off;
