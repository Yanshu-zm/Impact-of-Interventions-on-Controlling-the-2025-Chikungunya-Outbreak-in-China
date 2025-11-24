function [cities, FlowNet, pss] = GetFlowNetV2()

%%
file_name = "data/广东与外省地市间人口流动（20250728_0803）.csv";
opts = detectImportOptions(file_name);
opts.VariableNamingRule = 'preserve';  % 关键设置：保留原始列标题
data = readtable(file_name, opts);
    
% if ~strcmp(target_col,'输出城市') &&  ~strcmp(target_col,'输入城市')
%     error('wrong target_col')
% end

data.('输入城市') = strcat( data.('到达省份'), data.('到达城市'));
data.('输出城市') = strcat( data.('出发省份'), data.('出发城市'));

cities = unique(data.("输入城市"));
is_guangdong = startsWith(cities, '广东省');

guangdong_cities = cities(is_guangdong);
other_cities = cities(~is_guangdong);

cities = [guangdong_cities; other_cities];

ncities = length(cities);
FlowNet = zeros(ncities,ncities,7);  %% 源头，目的地

data.weekday = weekday(datetime(num2str(data.("日期")), 'Format', 'yyyyMMdd'));

for wd = 1: 7
    index = data.weekday == wd;
    flow_data_wd = data(index,:);
    for n = (1:ncities)
        flow_data = flow_data_wd(strcmp(flow_data_wd.("输入城市"),cities{n}),:);
        for j = 1:size(flow_data,1)
            source_citytmp = flow_data.("输出城市")(j);
            flow_volumn = flow_data.("人次")(j);
            index = strcmp(cities, source_citytmp{1});
            FlowNet(index,n,wd) = flow_volumn;
        end
    end
end

%% fill inner Guangdong flow
[citiesV1, flow_matrix, pss] = GetFlowNet();
GD_V1 = {'东莞市', '佛山市' , '广州市' ,'深圳市'};
GD_V2 = {'广东省东莞市', '广东省佛山市' , '广东省广州市' ,'广东省深圳市'};
is_in_list_v1 = ismember(citiesV1,GD_V1);
is_in_list_v2 = ismember(cities,GD_V2);
FlowNet(is_in_list_v2,is_in_list_v2,:) =  repmat(flow_matrix(is_in_list_v1,is_in_list_v1), [1, 1, 7]);

%%
pss = readtable('数据/常住人口_v2.csv');
pss = table2array(pss(:,2)) * 1e4;
pss = cat(1,pss,ones(length(ncities)-length(pss),1));
end