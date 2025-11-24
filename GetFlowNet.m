
function [cities, flow_matrix,pss] = GetFlowNet()
    cities = GetCityNames();
    ncites = length(cities);
    
    % make contagiout matrix
    flow_matrix = zeros(ncites,ncites); %% 源头，目的地
    for n = (1:ncites)
        flow_data = GetImportFlowNation(cities{n},"输入城市");
        for j = 1:size(flow_data,1)
            source_citytmp = flow_data.("输出城市")(j);
            flow_volumn = flow_data.("七天累计人流数")(j)/7;
            index = strcmp(cities, source_citytmp{1});
            flow_matrix(index,n) = flow_volumn;
        end
    end
    pss = readtable('数据/常住人口.csv');
    pss = table2array(pss(:,2)) * 1e4;
end