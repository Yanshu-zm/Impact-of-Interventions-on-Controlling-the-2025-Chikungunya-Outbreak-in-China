function flow_data = GetImportFlowNation(target_city,target_col)

%%
file_name = "data/flow_20241008_20241014.xlsx";
opts = detectImportOptions(file_name);
opts.VariableNamingRule = 'preserve';  % 关键设置：保留原始列标题
data = readtable(file_name, opts);
    
if ~strcmp(target_col,'输出城市') &&  ~strcmp(target_col,'输入城市')
    error('wrong target_col')
end
% disp('原始列名：');
% disp(data.Properties.VariableNames);
filter_logic = strcmp(data.(target_col), target_city);

% 应用筛选条件，得到符合条件的子表格
flow_data = data(filter_logic, :);
end