function [dateindex, temp_num,rainfall_num] = read_tem_jiangmen()
    
    filenames = {
        'data/Jiangmen_202507.xlsx', ...
        'data/Jiangmen_202508.xlsx', ...
        'data/Jiangmen_202509.xlsx', ...
        'data/Jiangmen_202510.xlsx', ...
        'data/Jiangmen_202410.xlsx', ...
        'data/Jiangmen_202411.xlsx', ...
        'data/Jiangmen_202412.xlsx', ...
        'data/Jiangmen_202501.xlsx'};

    T7 = readtable(filenames{1},VariableNamingRule='preserve');
    T8 = readtable(filenames{2},VariableNamingRule='preserve');
    T9 = readtable(filenames{3},VariableNamingRule='preserve');
    T10 = readtable(filenames{4},VariableNamingRule='preserve');
    T10p = readtable(filenames{5},VariableNamingRule='preserve');
    T11p = readtable(filenames{6},VariableNamingRule='preserve');
    T12p = readtable(filenames{7},VariableNamingRule='preserve');
    T1p = readtable(filenames{8},VariableNamingRule='preserve');

    T10p = T10p(T10p.("日期") >= datetime(2024, 10, 24),:);
    
    %最后会翻转
    T = [T1p;T12p;T11p;T10p;T10;T9];

    temp_str = T.("平均温度");
    rainfall = T.("24h降水量");
    dateindex = T.("日期");
    temp_num = zeros(height(T), 1);  
    rainfall_num =  zeros(height(T), 1);
    
    % 4. 遍历每个字符串，提取数字并转换
    for i = 1:height(T)
        str = temp_str{i};  % 取出单元格的字符串
    
        % 正则表达式匹配：
        % [-+]?  ：可选的正负号
        % \d+    ：1个或多个整数部分
        % \.?    ：可选的小数点
        % \d*    ：0个或多个小数部分
        num_match = regexp(str, '[-+]?\d+\.?\d*', 'match');
        
        if ~isempty(num_match)  % 若匹配到数字
            temp_num(i) = str2double(num_match{1});  % 转换为数值
        else
            temp_num(i) = NaN;  % 无数字时设为NaN（可根据需求调整）
        end
    end
    
    for i = 1:height(T)
        str = rainfall{i};  % 取出单元格的字符串
        num_match = regexp(str, '[-+]?\d+\.?\d*', 'match');
        
        if ~isempty(num_match)  % 若匹配到数字
            rainfall_num(i) = str2double(num_match{1});  % 转换为数值
        else
            rainfall_num(i) = NaN;  % 无数字时设为NaN（可根据需求调整）
        end
    end
    rainfall_num = flip(rainfall_num); % 数据是从日期大到小的，转过来
    temp_num = flip(temp_num);         % 转过来
end