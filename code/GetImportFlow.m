function weekdayResult = GetImportFlow(target_city)

%%
file_name = "data/迁入+迁出人口统计表（2025.5.13-7.28）.xlsx";
opts = detectImportOptions(file_name);
opts.VariableNamingRule = 'preserve';  % 关键设置：保留原始列标题
data = readtable(file_name, opts);
    

% disp('原始列名：');
% disp(data.Properties.VariableNames);
filter_logic = strcmp(data.('地点'), target_city) & strcmp(data.('迁移类型'), '迁入');

% 应用筛选条件，得到符合条件的子表格
filtered_data = data(filter_logic, :);

% % 显示筛选结果
% disp('筛选后的结果：');
% disp(filtered_data);

%%
result = filtered_data(:, {'迁入日期', '人次'});

% 3. 将迁入日期转换为datetime类型（如果原先是字符串）
% 检查当前日期列类型
if iscell(result.('迁入日期'))
    % 若为字符串单元格数组，转换为datetime
    result.('迁入日期') = datetime(result.('迁入日期'), 'Format', 'yyyy-MM-dd');
elseif ischar(result.('迁入日期'))
    % 若为字符数组，转换为datetime
    result.('迁入日期') = datetime(result.('迁入日期'), 'Format', 'yyyy-MM-dd');
end
data = result;
%%
% Read the CSV file
% data = readtable('from_foshan_to_shenzhen.csv', 'HeaderLines', 0); % 0表示没有标题行

% 为列指定名称（根据实际数据结构）
data.Properties.VariableNames = {'Date', 'Flow'};
% Assuming the columns are named 'Date' and 'Flow' - adjust if needed
% If columns don't have headers, use:
% data = readtable('flow_data.csv', 'VariableNames', {'Date', 'Flow'});

% Convert date strings to datetime format
data.Date = datetime(data.Date, 'Format', 'yyyy/MM/dd');

% Extract weekday (1=Sunday, 2=Monday, ..., 7=Saturday in MATLAB)
% Convert to 1=Monday, 2=Tuesday, ..., 7=Sunday for better readability
% data.WeekdayNum = mod(weekday(data.Date ) - 2, 7) + 1;

% Create weekday names
weekdayNames = { 'Sunday','Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday',};
% data.Weekday = categorical(data.WeekdayNum, 1:7, weekdayNames);
data.Weekday = weekday(data.Date ) ;
% Extract week number (of the year)
% data.WeekNum = week(data.Date);

% Calculate mean and std for each weekday across all weeks
[weekdayGroups, weekdayNames] = findgroups(data.Weekday);
weekdayStats = splitapply(@(x) [mean(x) std(x)], data.Flow, weekdayGroups);

% Calculate mean and std for each week
% [weekGroups, weekNums] = findgroups(data.WeekNum);
% weekStats = splitapply(@(x) [mean(x) std(x)], data.Flow, weekGroups);

% Display weekday results
% disp('Weekday Statistics (Mean and Std):');
weekdayResult = table(weekdayNames, weekdayStats(:,1), weekdayStats(:,2), ...
     'VariableNames', {'Weekday', 'MeanFlow', 'StdFlow'});
% disp(weekdayResult);

% Display weekly results
% disp('\nWeekly Statistics (Mean and Std):');
% weekResult = table(weekNums, weekStats(:,1), weekStats(:,2), ...
%     'VariableNames', {'WeekNumber', 'MeanFlow', 'StdFlow'});
% disp(weekResult);
% 
% % Optional: Visualize the results
% figure;
% bar(weekdayStats(:,1));
% hold on;
% errorbar(weekdayStats(:,1), weekdayStats(:,2), 'r.');
% set(gca, 'XTickLabel', weekdayNames, 'XTick', 1:7);
% title('Mean Flow by Weekday with Std Deviation');
% ylabel('Flow');
% xtickangle(45);
% grid on;
end