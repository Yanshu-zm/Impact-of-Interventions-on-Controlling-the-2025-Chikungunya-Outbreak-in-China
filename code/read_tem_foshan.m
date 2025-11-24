function [dateindex, temp_num,rainfall_num] = read_tem_foshan(filename)
    T = readtable(filename,VariableNamingRule='preserve');
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
    rainfall_num = flip(rainfall_num);
    temp_num = flip(temp_num);
end