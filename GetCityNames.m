function cities = GetCityNames()

    %%
    file_name = "数据/flow_20241008_20241014.xlsx";
    opts = detectImportOptions(file_name);
    opts.VariableNamingRule = 'preserve';  % 关键设置：保留原始列标题
    data = readtable(file_name, opts);
    
    cities = unique(data.("输出城市"));
end