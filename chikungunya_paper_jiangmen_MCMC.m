clear 
addpath(genpath("../utils"));
addpath(genpath("../dengue_parameters"));
addpath(genpath("../mcmcstat/mcmcstat-master"))

%% parameters
param = get_config();

param.ie = 0;%0.0535 / 365; %-310220 / 1400000000;
param.mu_h = 0;%7/1000; % Human death rate
param.mu_n = 0;%8/1000; % Human birth rate

ps_jiangmen = 4822600; %17500000;

param.sym_ratio = 0.9368; % 2012 CDC https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0012254
param.delta_h = 1/6; %i.e., rate of recovery.
param.gamma_h = 1/3; %Intrinsic incubation period of human
param.gamma_v = 1/8; %Intrinsic incubation period of vector
param.beta_v = 0.67;
param.beta_h = makedist('Uniform','Lower', 0.33,'Upper',0.37);

%% data 
% 
[dateindex,daily_temperature, daily_rainfall] = read_tem_jiangmen( );
a14days_rainfall = zeros(size(daily_rainfall));
num_days = length(daily_rainfall);

for i = 1:num_days
    start_idx = max(1, i - 14);
    end_idx = i;
    a14days_rainfall(i) = sum(daily_rainfall(start_idx:end_idx))*14 /(end_idx - start_idx + 1);
end
a14days_rainfall = a14days_rainfall(1:num_days);

z = 0.28; %%%
Rmin = 1;
Rmax = 123;
carrying_capacity = carrying_capacity_Tpart(daily_temperature, param) .* ...
    carrying_capacity_Rpart_Briere(a14days_rainfall, Rmin, Rmax, z);
%%carrying_capacity = smooth(carrying_capacity,3);
carrying_capacity = smoothdata(carrying_capacity, 'movmean', 3);
disp("-----size of carrying_capacity-----");
disp(size(carrying_capacity));
assert(size(carrying_capacity,2)==1)
%% running disease model in foshan/ sources
source_city = '江门';
if strcmp(source_city, '江门')
    daily_case = readtable('数据/jiangmen_Daily_cases.xlsx',VariableNamingRule='preserve');
    if ~isa(daily_case.Date, 'datetime')
        daily_case.Date = datetime(daily_case.Date);  % 自动解析日期格式
    end
    
    data_begin_date = datetime(2025, 9, 20);
    local_infection = daily_case.Cases(daily_case.Date >= data_begin_date, :)';
    dayofLocalObs = datenum(2025, 9, 20) - datenum(2025, 9, 1) + 1;

    start_day = 4;
    startDate = datetime(2025, 9, start_day)
    local_infection_dateindex = (dayofLocalObs:dayofLocalObs+length(local_infection)-1);
    dayOfStartDate = datenum(2025, 9, start_day) - datenum(2025, 9, 1) + 1;
    TIMELENGHT = length(dateindex);
    param.infection_rate_decline_begin1 = days(datetime(2025, 9, 19) - datetime(2025, 9, 1)) + 1;
    param.infection_rate_decline_begin2 = days(datetime(2025, 10, 6) - datetime(2025, 9, 1)) + 1;
end

param.mu_v_ratio = 1;
EPS = 1e-6;
use_mc = false;
if (use_mc)
    % MC
    num_iteration = 1e3;
    param_chain = zeros(num_iteration,3); % init_infection, infection_rate_decline1, infection_rate_decline2
    TIMELENGHT = length(dateindex);

    mu_v_increase_prior1 = 4;
    mu_v_increase_sigma1 = 10;
    mu_v_increase_prior2 = 8;
    mu_v_increase_sigma2 = 10;
    smooth_rate = 6000;
    prior_mean_ii = 4000;
    prior_sigma_ii = 1000;
    
    TIMELENGHT = 153;
    pre_log_lik = -inf;
    acc_rate = 0;
    ModelingOutputs = zeros(TIMELENGHT,7,num_iteration);
    NewInfections = zeros(TIMELENGHT,num_iteration);
    pre_MInfections = zeros(TIMELENGHT,7);
    pre_NewINfection = zeros(TIMELENGHT,1);
    f = waitbar(0, 'Starting');

    for k = 1 : num_iteration
        waitbar(k/num_iteration, f, sprintf('Progress: %d %%', floor(k/num_iteration*100)));
        import_infection_c = max(1,normrnd(prior_mean_ii/param.sym_ratio, prior_sigma_ii));
        mu_v_increase1 = max(0.00001,normrnd(mu_v_increase_prior1, mu_v_increase_sigma1));
        param.mu_v_increase1 = mu_v_increase1;
        
        mu_v_increase2 = max(0.00001,normrnd(mu_v_increase_prior2, mu_v_increase_sigma2));
        param.mu_v_increase2 = mu_v_increase2;

        [NewInfection, R0_array, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c,daily_temperature,a14days_rainfall,carrying_capacity,param);
        log_lik = - sum((NewInfection(dayofLocalObs:dayofLocalObs+length(local_infection)-1) - local_infection').^2)/length(local_infection) / smooth_rate;
        assert(isscalar(log_lik),"error_size_loglik")
        sample_rate = rand();
        % disp(strcat(num2str(k),':','pre_log_lik',num2str(pre_log_lik),';log_lik:',num2str(log_lik),'----',num2str(exp(log_lik - pre_log_lik))))
        if sample_rate < exp(log_lik - pre_log_lik) % NewInfection(dayOfYear)>0 &&
            pre_log_lik = log_lik;
            acc_rate = acc_rate+1;
            ModelingOutputs(:,:,k) = ModelingOutput;
            NewInfections(:,k) = NewInfection;
            pre_MInfections = ModelingOutput;
            pre_NewINfection = NewInfection;
            param_chain(k,:) = [import_infection_c, param.mu_v_increase1, param.mu_v_increase2];
        else
            ModelingOutputs(:,:,k) = pre_MInfections;
            NewInfections(:,k) = pre_NewINfection;
            param_chain(k,:) =  param_chain(k-1,:);
        end
    end
    close(f);
   
    disp(mean(param_chain(:,1),1))
    disp(mean(param_chain(:,2:3),1))

    acc_rate
else
    % MCMC
    mu_v_increase_prior1 = 4;
    mu_v_increase_sigma1 = 5;
    mu_v_increase_prior2 = 8;
    mu_v_increase_sigma2 = 5;

    prior_mean_ii = 4000;
    prior_sigma_ii = 1000;
    smooth_rate = 3000;

    model.ssfun      = @ssfun;
    options.nsimu    = 10000; 
    options.adaptint = 2000;
    %  {'par2',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
    mcmc_params = {
        {'init_infection', prior_mean_ii, 0.01*prior_mean_ii, 10*prior_mean_ii,prior_mean_ii,prior_sigma_ii};
        {'mu_v_increase1', mu_v_increase_prior1, 1, 10, mu_v_increase_prior1, mu_v_increase_sigma1};
        {'mu_v_increase2', mu_v_increase_prior2, 1, 10, mu_v_increase_prior2, mu_v_increase_sigma2};
    };
    % mcmc_params = {
    %     {'init_infection', prior_mean_ii, 0.01*prior_mean_ii, 10*prior_mean_ii};
    %     {'mu_v_increase1', mu_v_increase_prior1, 1, 10, mu_v_increase_prior1};
    %     {'mu_v_increase2', mu_v_increase_prior2, 1, 10, mu_v_increase_prior2};
    % };

    data = struct();
    data.ps_jiangmen = ps_jiangmen;              % 模型配置
    data.dayOfStartDate = dayOfStartDate;        % 模拟起始日期索引
    data.daily_temperature = daily_temperature;    % 每日温度
    data.a14days_rainfall = a14days_rainfall;      % 14天降雨量
    data.carrying_capacity = carrying_capacity;    % 承载能力
    data.param = param;                           % 固定参数集
    data.dayofLocalObs = dayofLocalObs;           % 本地观测起始索引
    data.local_infection = local_infection;       % 本地感染观测值
    data.smooth_rate = smooth_rate;               % 平滑系数
    [res,chain] = mcmcrun(model,data,mcmc_params,options);
    param_chain = chain;
    figure(1); clf; mcmcplot(chain);
    %%
    num_iteration = options.nsimu;
    out = mcmcpred(res,chain(options.adaptint:end,:),[],data,@f_model,num_iteration); %
    run_sims = out.ysaveout{1,1}{1}; % niter, time, nstates
    run_sims = permute(run_sims,[2,3,1]); % to (TIMELENGHT,7,num_iteration);
    NewInfections = squeeze(run_sims(:,1,:));
    ModelingOutputs = run_sims(:,2:end,:);

    % ModelingOutputs = zeros(TIMELENGHT,7,num_iteration);
    % NewInfections = zeros(TIMELENGHT,num_iteration);
end
%%
currentDateTime = datetime('now');
currentDateTimeString = datestr(currentDateTime, 'yyyy-mm-dd-HH-MM-SS');
%%
directory_name = "result/result_jiangmen";
export_file_name = strcat(directory_name, "/fitting_chain-", currentDateTimeString,".mat");
save(export_file_name,"param_chain");

%% fitting performance
if (true)
    f = figure();
    set(f,"Position",[1000,1007,560,230]);
    hold on;
    CI_plot(mean(NewInfections'), prctile(NewInfections',5)  , prctile(NewInfections',95))
    scatter(local_infection_dateindex, local_infection,Marker=".",color='red');
    %legend('Obs', 'Modeling');
    %xlabel('Days');
    xlim([1,153]);
    ylabel('Number of infections');
    
    tickDates =  datetime(2025, 9, 0) + caldays((1:10:153));
    xticks((1:10:153))
    xticklabels(datestr(tickDates, 'mm-dd'));
    xtickangle(45);

    directory_name = "result/result_jiangmen";
    export_file_name = strcat(directory_name, "/simulation_fitting-", currentDateTimeString);
    saveas(f, strcat(export_file_name, ".fig"));
    exportgraphics(f,strcat(export_file_name, ".pdf"));
end

%% Simulations
if (true)
    % Intervention one week early
    ModelingOutputs_owh = zeros(TIMELENGHT,7,num_iteration);
    NewInfections_owh = zeros(TIMELENGHT,num_iteration);

    param_new = param;
    param_new.infection_rate_decline_begin1 = param_new.infection_rate_decline_begin1 - 7; %% C
    param_new.infection_rate_decline_begin2 = param_new.infection_rate_decline_begin2 - 7; %% C
    for k = 1 : num_iteration
        import_infection_c =  param_chain(k,1);
        param_new.mu_v_increase1 = param_chain(k,2);
        param_new.mu_v_increase2 = param_chain(k,3);

        [NewInfection, R0_array, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c,daily_temperature,a14days_rainfall,carrying_capacity, ...
            param_new);
        ModelingOutputs_owh(:,:,k) = ModelingOutput;
        NewInfections_owh(:,k) = NewInfection;
    end
%% Intervention TWO week early
    ModelingOutputs_twh = zeros(TIMELENGHT,7,num_iteration);
    NewInfections_twh = zeros(TIMELENGHT,num_iteration);

    param_new = param;
    param_new.infection_rate_decline_begin1 = param_new.infection_rate_decline_begin1 - 14; %% C
    param_new.infection_rate_decline_begin2 = param_new.infection_rate_decline_begin2 - 14; %% C
    for k = 1 : num_iteration
        import_infection_c =  param_chain(k,1);
        param_new.mu_v_increase1 =  param_chain(k,2);
        param_new.mu_v_increase2 = param_chain(k,3);

        [NewInfection, R0_array, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c,daily_temperature,a14days_rainfall,carrying_capacity, ...
            param_new);
        ModelingOutputs_twh(:,:,k) = ModelingOutput;
        NewInfections_twh(:,k) = NewInfection;
    end

  %% Intervention One week lately
    ModelingOutputs_l1w = zeros(TIMELENGHT,7,num_iteration);
    NewInfections_l1w   = zeros(TIMELENGHT,num_iteration);
    
    param_new = param;
   
    param_new.infection_rate_decline_begin1 = param_new.infection_rate_decline_begin1 + 7;  % +7
    param_new.infection_rate_decline_begin2 = param_new.infection_rate_decline_begin2 + 7;
    
    for k = 1:num_iteration
        import_infection_c        = param_chain(k,1);
        param_new.mu_v_increase1  = param_chain(k,2);
        param_new.mu_v_increase2  = param_chain(k,3);
    
        [NewInfection, R0_array, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c,daily_temperature,a14days_rainfall,carrying_capacity, ...
            param_new);
        ModelingOutputs_l1w(:,:,k) = ModelingOutput;
        NewInfections_l1w(:,k)     = NewInfection;
    end
    
  %% Intervention Two week lately
    ModelingOutputs_l2w = zeros(TIMELENGHT,7,num_iteration);
    NewInfections_l2w   = zeros(TIMELENGHT,num_iteration);
    
    param_new = param;
    param_new.infection_rate_decline_begin1 = param_new.infection_rate_decline_begin1 + 14;  % +14
    param_new.infection_rate_decline_begin2 = param_new.infection_rate_decline_begin2 + 14;
    
    for k = 1:num_iteration
        import_infection_c        = param_chain(k,1);
        param_new.mu_v_increase1  = param_chain(k,2);
        param_new.mu_v_increase2  = param_chain(k,3);
    
        [NewInfection, R0_array, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c,daily_temperature,a14days_rainfall,carrying_capacity, ...
            param_new);
        ModelingOutputs_l2w(:,:,k) = ModelingOutput;
        NewInfections_l2w(:,k)     = NewInfection;

    end
%% Intervention strength get 1.5x
    ModelingOutputs_s15 = zeros(TIMELENGHT,7,num_iteration);
    NewInfections_s15 = zeros(TIMELENGHT,num_iteration);

    param_new = param;
    for k = 1 : num_iteration
        import_infection_c =  param_chain(k,1);
        param_new.mu_v_increase1 = param_chain(k,2) * 1.5; %% C
        param_new.mu_v_increase2 = param_chain(k,3) * 1.5; %% C

        [NewInfection, R0_array, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c,daily_temperature,a14days_rainfall,carrying_capacity, ...
            param_new);
        ModelingOutputs_s15(:,:,k) = ModelingOutput;
        NewInfections_s15(:,k) = NewInfection;
    end
%% Intervention strength get 1.5x & Intervention one week early
    ModelingOutputs_s15_owh = zeros(TIMELENGHT,7,num_iteration);
    NewInfections_s15_owh = zeros(TIMELENGHT,num_iteration);

    param_new = param;
    param_new.infection_rate_decline_begin1 = param_new.infection_rate_decline_begin1 - 7; %% C
    param_new.infection_rate_decline_begin2 = param_new.infection_rate_decline_begin2 - 7; %% C

    for k = 1 : num_iteration
        import_infection_c =  param_chain(k,1);
        param_new.mu_v_increase1 = param_chain(k,2) * 1.5; %% C
        param_new.mu_v_increase2 = param_chain(k,3) * 1.5; %% C

        [NewInfection, R0_array, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c,daily_temperature,a14days_rainfall,carrying_capacity, ...
            param_new);
        ModelingOutputs_s15_owh(:,:,k) = ModelingOutput;
        NewInfections_s15_owh(:,k) = NewInfection;
    end
 %% 新增：0.5 倍强度
    ModelingOutputs_s05 = zeros(TIMELENGHT,7,num_iteration);
    NewInfections_s05   = zeros(TIMELENGHT,num_iteration);
    
    param_new = param;
    for k = 1 : num_iteration
        import_infection_c       = param_chain(k,1);
        param_new.mu_v_increase1 = param_chain(k,2) * 0.5;   % 0.5×
        param_new.mu_v_increase2 = param_chain(k,3) * 0.5;
        [NewInfection, ~, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c,...
                                   daily_temperature,a14days_rainfall,carrying_capacity,param_new);
        ModelingOutputs_s05(:,:,k) = ModelingOutput;
        NewInfections_s05(:,k)     = NewInfection;
    end
    
 %% 新增：0.5 倍强度 + 提前一周
    ModelingOutputs_s05_owh = zeros(TIMELENGHT,7,num_iteration);
    NewInfections_s05_owh   = zeros(TIMELENGHT,num_iteration);
    
    param_new = param;
    param_new.infection_rate_decline_begin1 = param_new.infection_rate_decline_begin1 + 7;
    param_new.infection_rate_decline_begin2 = param_new.infection_rate_decline_begin2 + 7;
    
    for k = 1 : num_iteration
        import_infection_c       = param_chain(k,1);
        param_new.mu_v_increase1 = param_chain(k,2) * 0.5;
        param_new.mu_v_increase2 = param_chain(k,3) * 0.5;
        [NewInfection, ~, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c,...
                                   daily_temperature,a14days_rainfall,carrying_capacity,param_new);
        ModelingOutputs_s05_owh(:,:,k) = ModelingOutput;
        NewInfections_s05_owh(:,k)     = NewInfection;
    end

    %% 10. 新增：完全无干预场景
    % —— 蚊媒额外死亡率 = 1（不额外增加）
    % —— 干预开始时间 > 模拟长度，永远不会触发
    param_noIntv = param;
    param_noIntv.mu_v_increase1 = 0;
    param_noIntv.mu_v_increase2 = 0;
    param_noIntv.infection_rate_decline_begin1 = 999;
    param_noIntv.infection_rate_decline_begin2 = 999;
    
    % 预分配
    ModelingOutputs_noIntv = zeros(TIMELENGHT,7,num_iteration);
    NewInfections_noIntv   = zeros(TIMELENGHT,num_iteration);
    
    for k = 1:num_iteration
        import_infection_c       = param_chain(k,1);   % 用同一条链
        [NewInfection, ~, ModelingOutput] = simulate( ...
            ps_jiangmen,dayOfStartDate,import_infection_c, ...
            daily_temperature,a14days_rainfall,carrying_capacity, ...
            param_noIntv);
        ModelingOutputs_noIntv(:,:,k) = ModelingOutput;
        NewInfections_noIntv(:,k)   = NewInfection;
    end
end
%%
currentDateTime = datetime('now');
currentDateTimeString = datestr(currentDateTime, 'yyyy-mm-dd-HH-MM-SS');
scenarioCell = {'Base','owh','twh','s15','s15_owh','l1w','l2w','s05','s05_owh','noIntv'};
dataCell     = {NewInfections, NewInfections_owh, NewInfections_twh, ...
                NewInfections_s15, NewInfections_s15_owh, ...
                NewInfections_l1w, NewInfections_l2w, ...
                NewInfections_s05, NewInfections_s05_owh,NewInfections_noIntv};

cumAll = cellfun(@(x)sum(x,1), dataCell, 'UniformOutput', false); % 9×1 cell，每 cell 1×1000
ref    = cumAll{1};                 % 以 Base 为参照
redTbl = table();             % 准备输出表
for i = 2:numel(scenarioCell)
    scen = scenarioCell{i};
    y    = cumAll{i};         % 当前场景 1×1000
    
    pct = y./ ref * 100;  %SI中计算基于base的case倍数

    % 汇总 median + 95% CI
    redTbl.(scen) = [median(pct); prctile(pct,[2.5 97.5])'];
end
writetable(redTbl, fullfile(directory_name, ...
              sprintf('Local_cases_compared_with_baseline_jm.xlsx', currentDateTimeString)));

for i = 2:numel(scenarioCell)
    scen = scenarioCell{i};
    y    = cumAll{i};         % 当前场景 1×1000
    
    if any(strcmp(scen, {'owh','twh','s15','s15_owh'}))
        pct2 = (ref - y) ./ ref * 100;                 % 减少百分比

    elseif any(strcmp(scen, {'l1w','l2w','s05','s05_owh'}))
        pct2 = (y - ref) ./ ref * 100;              % 增加百分比（相对 Base）

    elseif strcmp(scen, 'noIntv')
        pct2 = (y - ref) ./ y * 100;  


    end

    % 汇总 median + 95% CI
    redTbl.(scen) = [median(pct2); prctile(pct2,[2.5 97.5])'];
end
writetable(redTbl, fullfile(directory_name, ...
              sprintf('cumReductionPct_95CI_%s.xlsx', currentDateTimeString)));

%% ======  统一刻度  ======
nDay      = TIMELENGHT;                    % 
tickStep  = 10;                            % 每 10 天一个刻度
tickIdx   = 1:tickStep:nDay;               % [1 11 21 ... 71]
tickDates = datetime(2025,9,1) + caldays(tickIdx-1);  % 对应日历
%% ----- 读入江门真实每日病例 -----
realTbl   = readtable('数据/jiangmen_real_daily_case.xlsx','Sheet',1);   % 文件就在当前目录
excelDate = datenum(realTbl.date);                                      % Excel 序号日期
realCase  = realTbl.jiangmen;                                          % 病例数
% 转成"2025-07-01 为第 1 天"的索引
refDay    = datenum(2025,9,1);                                     % 模型第 1 天
realDay   = excelDate - refDay + 1;                                % 与模型 time_index 完全对齐
% 只保留落在模拟区间内的点
valid     = realDay >= 1 & realDay <= TIMELENGHT;
realDay   = realDay(valid);
realCase  = realCase(valid);
%% figure2a
line_width = 1;
ci_alpha = 0.3;  % 置信区间透明度
%time_index = time_index';
time_index = 1:TIMELENGHT;
f = figure();
set(f,"Position",[1000,1007,560,230]);
hold on;grid minor;

f1= [0 0 0];
f2 = [56,152,211]/255;
f3 = [0,73,146]/255;
f4 = [237,90,96]/255;
f5 = [169,14,16]/255;
f6 = [0 0 0];  
f7 = [88,0,1]/255;
% --- 原曲线 ---
hb0 = fillyy(time_index, ...
             prctile(NewInfections',5),  ...
             prctile(NewInfections',95), ...
             f1, 0.5);
plot(mean(NewInfections'), 'Color', f1, 'LineWidth', 1.5);

hb5 = fillyy(time_index, ...
             prctile(NewInfections_noIntv',5),  ...
             prctile(NewInfections_noIntv',95), ...
             f6, 0.15);
plot(mean(NewInfections_noIntv'), 'Color', f6, 'LineWidth', 1.5,'LineStyle','--');

hb1 = fillyy(time_index, ...
             prctile(NewInfections_owh',5),  ...
             prctile(NewInfections_owh',95), ...
             f2, 0.25);
plot(mean(NewInfections_owh'), 'Color', f2, 'LineWidth', 1.5);

hb2 = fillyy(time_index, ...
             prctile(NewInfections_twh',5),  ...
             prctile(NewInfections_twh',95), ...
             f3, 0.25);
plot(mean(NewInfections_twh'), 'Color', f3, 'LineWidth', 1.5);

% --- 新增延迟曲线 ---
hb3 = fillyy(time_index, ...
             prctile(NewInfections_l1w',5),  ...
             prctile(NewInfections_l1w',95), ...
             f4, 0.25);
plot(mean(NewInfections_l1w'), 'Color', f4, 'LineWidth', 1.5);

hb4 = fillyy(time_index, ...
             prctile(NewInfections_l2w',5),  ...
             prctile(NewInfections_l2w',95), ...
             f5, 0.25);
plot(mean(NewInfections_l2w'), 'Color', f5, 'LineWidth', 1.5);

plot(realDay, realCase, 'Color',f7, 'Marker','o','MarkerSize',5,'LineWidth', 1.5, 'LineStyle', 'none');


% % ---------- legend & axes ----------
% legend([hb0,hb1,hb2,hb3,hb4,hb5], ...
%        {' Actual Intervention','One week early','Two weeks early', ...
%         'One week late','Two weeks late','Without Intervention'}, ...
%        'Location','northeastoutside');
xlim([1 nDay]);
xticks(tickIdx);
xticklabels(datestr(tickDates,'dd-mmm'));
xtickangle(45);
ylabel('Number of daily infections');
ylim([0 6000]);  
% 导出
export_file_name = strcat(directory_name, "/week_early_late-", currentDateTimeString);
saveas(f, strcat(export_file_name, ".fig"));
exportgraphics(f, strcat(export_file_name, ".pdf"), 'ContentType', 'vector');

%% figure2b
f = figure();
set(f,"Position",[1000,1007,560,230]);
hold on;grid minor;

% 颜色定义
% gray = [0.7 0.7 0.7];
% green= [0 0.6 0];
% magen= [0.8 0.2 0.8];
% blue = [0 0.4 0.8];

f1= [0 0 0];
f2 = [194,116,183]/255;
f3 = [132,25,119]/255;
f4 = [255,145,36]/255;
f5 = [194,64,0]/255;
f6 = [0 0 0]; 

% 1. 原拟合
hb0 = fillyy(time_index, ...
             prctile(NewInfections',5),  ...
             prctile(NewInfections',95), ...
             f1, 0.5);
plot(mean(NewInfections'),  'Color', f1, 'LineWidth', 1.5);

% ---------- 6. NEW: Without intervention ----------
hb5 = fillyy(time_index, ...
             prctile(NewInfections_noIntv',5),  ...
             prctile(NewInfections_noIntv',95), ...
             f6, 0.15);
plot(mean(NewInfections_noIntv'), 'Color', f6, 'LineWidth', 1.5,'LineStyle','--');

% 2. 1.5×
hb1 = fillyy(time_index, ...
             prctile(NewInfections_s15',5),  ...
             prctile(NewInfections_s15',95), ...
             f2, 0.25);
plot(mean(NewInfections_s15'), 'Color', f2, 'LineWidth', 1.5);

% 3. 1.5× + 提前一周
hb2 = fillyy(time_index, ...
             prctile(NewInfections_s15_owh',5),  ...
             prctile(NewInfections_s15_owh',95), ...
             f3, 0.25);
plot(mean(NewInfections_s15_owh'), 'Color', f3, 'LineWidth', 1.5);

% 4. 0.5×
hb3 = fillyy(time_index, ...
             prctile(NewInfections_s05',5),  ...
             prctile(NewInfections_s05',95), ...
             f4, 0.25);
plot(mean(NewInfections_s05'), 'Color', f4, 'LineWidth', 1.5);

% 5. 0.5× + 提前一周
hb4 = fillyy(time_index, ...
             prctile(NewInfections_s05_owh',5),  ...
             prctile(NewInfections_s05_owh',95), ...
             f5, 0.25);   % 浅青色
plot(mean(NewInfections_s05_owh'), 'Color', f5, 'LineWidth', 1.5);

plot(realDay, realCase, 'Color',f7, 'Marker','o','MarkerSize',5,'LineWidth', 1.5, 'LineStyle', 'none');

% legend([hb0,hb1,hb2,hb3,hb4,hb5], ...
%        {'Actual Intervention','Strength ×1.5','Strength ×1.5 + 1 week early', ...
%         'Strength ×0.5','Strength ×0.5 + 1 week late','Without Intervention'})
xlim([1 nDay]);
xticks(tickIdx);
xticklabels(datestr(tickDates,'dd-mmm'));
xtickangle(45);
ylabel('Number of daily infections');
ylim([0 6000]);          % 强制 y 轴 0–1000
export_file_name = strcat(directory_name, "/strength_05_15-", currentDateTimeString);
saveas(f, strcat(export_file_name, ".fig"));
exportgraphics(f, strcat(export_file_name, ".pdf"), 'ContentType', 'vector');
%% 
if (true)
     %%
    [cities, FlowNet, pss] = GetFlowNetV2();
    ncites = length(cities);
    JIANGMEN_INDEX = 11;
    disp(cities{JIANGMEN_INDEX})

    immigrant_risk_ratio_simulation = 0.1;
    %%
    num_other_cites = ncites-21; % 广东21省市在前面
    E_import_nums = zeros( num_other_cites, TIMELENGHT, num_iteration); % 江门源头，目的地，时间
    AI_import_nums = zeros( num_other_cites, TIMELENGHT, num_iteration); % 江门源头，目的地，时间
    
    for t = 9 : TIMELENGHT
       wd = weekday(datetime(2025,7,0)+caldays(t));
       flow = FlowNet(JIANGMEN_INDEX,22:end,wd);
    
       % cal the number of outputs
       prevalance = ModelingOutputs(t-1,5,:)./ps_jiangmen; % 改用E
       prevalance_asymI = (1-param.sym_ratio) * ModelingOutputs(t-1,6,:)./ps_jiangmen;
       E_import_nums(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
       AI_import_nums(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';

    end

    % f = figure();
    % set(f,"Position",[600,600,560,230]);
    % hold on;
    % time_index = (1:TIMELENGHT);
    % gray_color = [0.7 0.7 0.7];
    % 
    % import_nums = squeeze(sum(E_import_nums+AI_import_nums,1));
    % plot(squeeze(cumsum(mean(import_nums,2))),'k',LineWidth=1);
    % hb = mfillyy(time_index, ...
    %     squeeze(cumsum(prctile(import_nums',5))), ...
    %     squeeze(cumsum(prctile(import_nums',95))), ...
    %     gray_color, 0.3);

    %%
E_import_nums_owh = zeros( num_other_cites, TIMELENGHT, num_iteration);
AI_import_nums_owh = zeros( num_other_cites, TIMELENGHT, num_iteration); 

E_import_nums_twh = zeros( num_other_cites, TIMELENGHT, num_iteration); 
AI_import_nums_twh = zeros( num_other_cites, TIMELENGHT, num_iteration); 

E_import_nums_s15 = zeros( num_other_cites, TIMELENGHT, num_iteration); 
AI_import_nums_s15 = zeros( num_other_cites, TIMELENGHT, num_iteration); 

E_import_nums_s15_owh = zeros( num_other_cites, TIMELENGHT, num_iteration); 
AI_import_nums_s15_owh = zeros( num_other_cites, TIMELENGHT, num_iteration); 
% 1. 延迟 1 周
E_import_nums_l1w  = zeros(num_other_cites, TIMELENGHT, num_iteration);
AI_import_nums_l1w = zeros(num_other_cites, TIMELENGHT, num_iteration);
% 2. 延迟 2 周
E_import_nums_l2w  = zeros(num_other_cites, TIMELENGHT, num_iteration);
AI_import_nums_l2w = zeros(num_other_cites, TIMELENGHT, num_iteration);
% 3. 强度 0.5
E_import_nums_s05  = zeros(num_other_cites, TIMELENGHT, num_iteration);
AI_import_nums_s05 = zeros(num_other_cites, TIMELENGHT, num_iteration);
% 4. 强度 0.5 + 提前 1 周
E_import_nums_s05_owh  = zeros(num_other_cites, TIMELENGHT, num_iteration);
AI_import_nums_s05_owh = zeros(num_other_cites, TIMELENGHT, num_iteration);

%% 无干预输入风险
E_import_nums_noIntv  = zeros(num_other_cites, TIMELENGHT, num_iteration);
AI_import_nums_noIntv = zeros(num_other_cites, TIMELENGHT, num_iteration);
for t = 9:TIMELENGHT
    wd = weekday(datetime(2025,9,0)+caldays(t));
    flow = FlowNet(JIANGMEN_INDEX,22:end,wd);
    prevalance     = ModelingOutputs_noIntv(t-1,5,:)./ps_jiangmen;
    prevalance_asy = (1-param.sym_ratio)*ModelingOutputs_noIntv(t-1,6,:)./ps_jiangmen;
    E_import_nums_noIntv(:,t,:)  = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)';
    AI_import_nums_noIntv(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asy)';
end

for t = 9 : TIMELENGHT
   wd = weekday(datetime(2025,9,0)+caldays(t));
   flow = FlowNet(JIANGMEN_INDEX,22:end,wd);

   % cal the number of outputs
   prevalance = ModelingOutputs(t-1,5,:)./ps_jiangmen; % 改用E
   prevalance_asymI = (1-param.sym_ratio) * ModelingOutputs(t-1,6,:)./ps_jiangmen;
   E_import_nums(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
   AI_import_nums(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';

   % one week +
   prevalance = ModelingOutputs_owh(t-1,5,:)./ps_jiangmen; % 改用E
   prevalance_asymI = (1-param.sym_ratio) * ModelingOutputs_owh(t-1,6,:)./ps_jiangmen;
   E_import_nums_owh(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
   AI_import_nums_owh(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';

   % two week +
   prevalance = ModelingOutputs_twh(t-1,5,:)./ps_jiangmen; % 改用E
   prevalance_asymI = (1-param.sym_ratio) * ModelingOutputs_twh(t-1,6,:)./ps_jiangmen;
   E_import_nums_twh(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
   AI_import_nums_twh(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';
   

   % s1.5
   prevalance = ModelingOutputs_s15(t-1,5,:)./ps_jiangmen; % 改用E
   prevalance_asymI = (1-param.sym_ratio) * ModelingOutputs_s15(t-1,6,:)./ps_jiangmen;
   E_import_nums_s15(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
   AI_import_nums_s15(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';

   % s1.5 + one week
   prevalance = ModelingOutputs_s15_owh(t-1,5,:)./ps_jiangmen; % 改用E
   prevalance_asymI = (1-param.sym_ratio) * ModelingOutputs_s15_owh(t-1,6,:)./ps_jiangmen;
   E_import_nums_s15_owh(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
   AI_import_nums_s15_owh(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';

    % ---- 延迟 1 周 ----
   prevalance = ModelingOutputs_l1w(t-1,5,:)./ps_jiangmen;
   prevalance_asymI = (1-param.sym_ratio)*ModelingOutputs_l1w(t-1,6,:)./ps_jiangmen;
   E_import_nums_l1w(:,t,:)  = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
   AI_import_nums_l1w(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';

    % ---- 延迟 2 周 ----
   prevalance = ModelingOutputs_l2w(t-1,5,:)./ps_jiangmen;
   prevalance_asymI = (1-param.sym_ratio)* ModelingOutputs_l2w(t-1,6,:)./ps_jiangmen;
   E_import_nums_l2w(:,t,:)  = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
   I_import_nums_l2w(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';
    
    % ---- 强度 0.5 ----
   prevalance = ModelingOutputs_s05(t-1,5,:)./ps_jiangmen;
   prevalance_asymI = (1-param.sym_ratio)*ModelingOutputs_s05(t-1,6,:)./ps_jiangmen;
   E_import_nums_s05(:,t,:)  = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
   AI_import_nums_s05(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';

    % ---- 强度 0.5 + 提前 1 周 ----
   prevalance = ModelingOutputs_s05_owh(t-1,5,:)./ps_jiangmen;
   prevalance_asymI = (1-param.sym_ratio)*ModelingOutputs_s05_owh(t-1,6,:)./ps_jiangmen;
   E_import_nums_s05_owh(:,t,:)  = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance)' ;
   AI_import_nums_s05_owh(:,t,:) = immigrant_risk_ratio_simulation * squeeze(flow)' * squeeze(prevalance_asymI)';

end
% %%
% importScenCell = {'Base','owh','twh','s15','s15_owh','l1w','l2w','s05','s05_owh','noIntv'};
% E_AI1 = E_import_nums   + AI_import_nums;
% E_AI2 = E_import_nums_owh + AI_import_nums_owh;
% E_AI3 = E_import_nums_twh + AI_import_nums_twh;
% E_AI4 = E_import_nums_s15 + AI_import_nums_s15;
% E_AI5 = E_import_nums_s15_owh + AI_import_nums_s15_owh;
% E_AI6 = E_import_nums_l1w + AI_import_nums_l1w;
% E_AI7 = E_import_nums_l2w + AI_import_nums_l2w;
% E_AI8 = E_import_nums_s05 + AI_import_nums_s05;
% E_AI9 = E_import_nums_s05_owh + AI_import_nums_s05_owh;
% E_AI10= E_import_nums_noIntv + AI_import_nums_noIntv;   % 新增
% importDataCell  = {E_AI1, E_AI2,E_AI3,E_AI4,E_AI5,E_AI6,E_AI7,E_AI8,E_AI9,E_AI10};
% 
% totImport = @(x) squeeze(sum(x,[1 2]));   % 返回 1×1000
% cumImportAll = cellfun(totImport, importDataCell, 'UniformOutput', false);
% % cumImportAll = cellfun(@(x) squeeze(sum(x,[1 2])), importDataCell, 'UniformOutput', false);
% refImp = cumImportAll{1};          % Base 作为参照
% redImpTbl = table();
% for jj = 2:numel(importScenCell)
%     pctImp = (refImp - cumImportAll{jj}) ./ refImp * 100;
%     redImpTbl.(importScenCell{jj}) = [median(pctImp); prctile(pctImp,[2.5 97.5])'];
% end
% writetable(redImpTbl, fullfile(directory_name, ...
%               sprintf('importCumReductionPct_allScen_%s.xlsx', currentDateTimeString)));

%%
%%figure3a

f = figure();
set(f,"Position",[1000,1007,560,230]);
hold on;grid minor;
time_index = (1:TIMELENGHT);


% 颜色顺序
f1= [0 0 0];
f2 = [56,152,211]/255;
f3 = [0,73,146]/255;
f4 = [237,90,96]/255;
f5 = [169,14,16]/255;
f6 = [0 0 0];  

% —— 1. 实际干预 ——
import_nums = squeeze(sum(E_import_nums+AI_import_nums,1));
yMean  = cumsum(mean(import_nums,2))+1;      % ← 取 log10
yLow   = cumsum(prctile(import_nums',5))+1;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb0 = fillyy(time_index, yLow, yHigh, f1, 0.5);
plot(time_index, yMean, 'Color', f1, 'LineWidth',1.5);

% —— 0. 无干预 ——
import_nums = squeeze(sum(E_import_nums_noIntv+AI_import_nums_noIntv,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hbNo = fillyy(time_index, yLow, yHigh, f1, 0.15);
plot(time_index, yMean, 'Color', f1, 'LineWidth',1.5, 'LineStyle','--');

% —— 2. 提前 1 周 ——
import_nums = squeeze(sum(E_import_nums_owh+AI_import_nums_owh,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb1 = fillyy(time_index, yLow, yHigh, f2, 0.25);
plot(time_index, yMean, 'Color', f2, 'LineWidth',1.5);

% —— 3. 提前 2 周 ——
import_nums = squeeze(sum(E_import_nums_twh+AI_import_nums_twh,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb2 = fillyy(time_index, yLow, yHigh, f3, 0.25);
plot(time_index, yMean, 'Color', f3, 'LineWidth',1.5);

% —— 4. 延迟 1 周 ——
import_nums = squeeze(sum(E_import_nums_l1w+AI_import_nums_l1w,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb3 = fillyy(time_index, yLow, yHigh, f4, 0.25);
plot(time_index, yMean, 'Color', f4, 'LineWidth',1.5);

% —— 5. 延迟 2 周 ——
import_nums = squeeze(sum(E_import_nums_l2w+AI_import_nums_l2w,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb4 = fillyy(time_index, yLow, yHigh, f5, 0.25);
plot(time_index, yMean, 'Color', f5, 'LineWidth',1.5);

% 修改 legend 顺序，把无干预放最前
legend([hbNo,hb0,hb1,hb2,hb3,hb4], ...
       {'Without intervention','Actual Intervention','One week early', ...
        'Two weeks early','One week late','Two weeks late'}, ...
       'Location','northeastoutside');
xlim([1 nDay]);
xticks(tickIdx);
xticklabels(datestr(tickDates,'dd-mmm')); xtickangle(45);
ylabel('Cumulative exported infections');
ylim([0 600]); 
export_file_name = strcat(directory_name, "/import_timing_all-", currentDateTimeString);
saveas(f, strcat(export_file_name, ".fig"));
exportgraphics(f, strcat(export_file_name, ".pdf"), 'ContentType', 'vector');

%%figure3b

f = figure();
set(f,"Position",[1000,1007,560,230]);
hold on;grid minor;
time_index = (1:TIMELENGHT);

% gray  = [0.7 0.7 0.7];
% green = [0 0.6 0];
% magen= [0.8 0.2 0.8];
% blue  = [0 0.4 0.8];
% lime  = [0.6 0.9 0.2];        % 0.5×
% sky   = [0.2 0.8 0.8];        % 0.5×+提前1周

f1= [0 0 0];
f2 = [194,116,183]/255;
f3 = [132,25,119]/255;
f4 = [255,145,36]/255;
f5 = [194,64,0]/255;
f6 = [0 0 0]; 

% —— 1. 实际干预 ——
import_nums = squeeze(sum(E_import_nums+AI_import_nums,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb0 = fillyy(time_index, yLow, yHigh, f1, 0.5);
plot(time_index, yMean, 'Color', f1, 'LineWidth',1.5);

% —— 0. 无干预 ——
import_nums = squeeze(sum(E_import_nums_noIntv+AI_import_nums_noIntv,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hbNo = fillyy(time_index, yLow, yHigh, f1, 0.15);
plot(time_index, yMean, 'Color', f1, 'LineWidth',1.5, 'LineStyle','--');

% —— 2. 1.5×强度 ——
import_nums = squeeze(sum(E_import_nums_s15+AI_import_nums_s15,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb1 = fillyy(time_index, yLow, yHigh, f2, 0.25);
plot(time_index, yMean, 'Color', f2, 'LineWidth',1.5);

% —— 3. 1.5× + 提前1周 ——
import_nums = squeeze(sum(E_import_nums_s15_owh+AI_import_nums_s15_owh,1)); 
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb2 = fillyy(time_index, yLow, yHigh, f3, 0.25);
plot(time_index, yMean, 'Color', f3, 'LineWidth',1.5);

% —— 4. 0.5×强度 ——
import_nums = squeeze(sum(E_import_nums_s05+AI_import_nums_s05,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb3 = fillyy(time_index, yLow, yHigh, f4, 0.25);
plot(time_index, yMean, 'Color', f4, 'LineWidth',1.5);

% —— 5. 0.5× + 提前1周 ——
import_nums = squeeze(sum(E_import_nums_s05_owh+AI_import_nums_s05_owh,1));
yMean  = cumsum(mean(import_nums,2))+1 ;
yLow   = cumsum(prctile(import_nums',5))+1 ;
yHigh  = cumsum(prctile(import_nums',95))+1 ;
hb4 = fillyy(time_index, yLow, yHigh, f5, 0.25);
plot(time_index, yMean, 'Color', f5, 'LineWidth',1.5);



legend([hbNo,hb0,hb1,hb2,hb3,hb4], ...
       {'Without intervention','Actual Intervention','Strength ×1.5', ...
        'Strength ×1.5 + 1 week early','Strength ×0.5','Strength ×0.5 + 1 week late'}, ...
       'Location','northeastoutside');
xlim([1 nDay]);
xticks(tickIdx);
xticklabels(datestr(tickDates,'dd-mmm')); xtickangle(45);
ylabel('Cumulative exported infections');
ylim([0 600]);
export_file_name = strcat(directory_name, "/import_strength_all-", currentDateTimeString);
saveas(f, strcat(export_file_name, ".fig"));
exportgraphics(f, strcat(export_file_name, ".pdf"), 'ContentType', 'vector');
end

%% ========== 新增：挑出 341 个广东外城市 ==========
guangdong21  = contains(cities,'广东省');      % 逻辑真值数组
outProvIdx   = find(~guangdong21);            % 341 个外省城市下标
outProvCities= cities(outProvIdx);            % 对应城市名

%% ---------- ----------
% 把 10 个场景的三维数组打包成 cell，方便循环调用
E_AI_Cell = { ...
    E_import_nums   + AI_import_nums, ...                 % 1  Actual
    E_import_nums_owh + AI_import_nums_owh, ...            % 2  owh
    E_import_nums_twh + AI_import_nums_twh, ...            % 3  twh
    E_import_nums_s15 + AI_import_nums_s15, ...            % 4  s15
    E_import_nums_s15_owh + AI_import_nums_s15_owh, ...    % 5  s15_owh
    E_import_nums_l1w + AI_import_nums_l1w, ...            % 6  l1w
    E_import_nums_l2w + AI_import_nums_l2w, ...            % 7  l2w
    E_import_nums_s05 + AI_import_nums_s05, ...            % 8  s05
    E_import_nums_s05_owh + AI_import_nums_s05_owh, ...    % 9  s05_owh
    E_import_nums_noIntv + AI_import_nums_noIntv};         % 10 noIntv

nCase = 10;
resultTable = table((1:nCase)', nan(nCase,1), ...
            'VariableNames',{'Case','NumCities_ge1'});

for k = 1:nCase
    E_AI    = E_AI_Cell{k};                  % 取出当前场景 341×62×1000
    cumLast = squeeze(sum(E_AI,2));          % 341×1000  每人累计
    cumMean = mean(cumLast,2);               % 341×1     均值
    idx1    = find(cumMean >= 1);            % 过线城市下标
    resultTable.NumCities_ge1(k) = length(idx1);
end

fn = sprintf('result/result_jiangmen/CrossCity_ge1_%s.xlsx', datestr(now,'yyyymmdd_HHMMSS'));
writetable(resultTable, fn);
fprintf('10 个场景的过线城市数已写入：%s\n', fn);

%% ---------- 过线城市 ----------
E_AI = E_import_nums   + AI_import_nums;
cumLast  = squeeze(sum(E_AI,2));                            % 341×1000，每人62天累计
cumMean = mean(cumLast,2);  
idx1 = find(cumMean >= 1);
name = outProvCities(idx1);

%% ---------- 画图：≥1 例 ----------
% 2. 预分配颜色（11 条，用原来灰-蓝-青-橙-红循环）
color11 = [
     
    % 126,0,2 ;  
          %   红
    169,14,16 ;      %   
    % 214,34,54;  
    237,90,96 ;      %    
    246,156,164; 
    %254,204,201;%  
   
    % 189,230,254 ;     %   
    142,204,236 ;        %   
    56,152,211 ;        %   
    %0,112,179;        %  
    0,73,146;  
      
    %0,36,93;
    ]/255;            % 蓝

finalCases = zeros(1,numel(idx1));
for k = 1:numel(idx1)
    iCity = idx1(k);
    cum3D = squeeze(sum(E_AI(iCity,:,:),1));  % 62×1000
    mu    = mean(cum3D,2);                  % 均值
    tmp   = cumsum(mu);                     % 先存变量
    finalCases(k) = tmp(end);               % 最后一天累积
end
[~, sortIdx] = sort(finalCases,'descend');           % 由小到大

% 重排一切
idx1      = idx1(sortIdx);
% color11   = color11(sortIdx,:);


f1 = figure('Color','w');
set(f1,"Position",[1000,1007,560,230]);
hold on;  grid minor;
time_index = 1:TIMELENGHT;

% 预分配句柄数组
h = gobjects(1,11);

for k = 1:length(idx1)
    iCity   = idx1(k);
    cum3D   = squeeze(sum(E_AI(iCity,:,:),1));  % 62×1000
    mu      = mean(cum3D,2);
    low95   = prctile(cum3D',5);
    up95    = prctile(cum3D',95);

 % 1) 均值线——要句柄
    h(k) = plot(time_index, cumsum(mu), ...
                'Color', color11(k,:), ...
                'LineWidth', 2);

    % 2) 置信带——不要句柄
    fillyy(time_index, ...
           cumsum(low95), ...
           cumsum(up95), ...
           color11(k,:), 0.2);
end


% 5. 坐标轴装饰
xlim([1 nDay]);
xticks(tickIdx);
xticklabels(datestr(tickDates,'dd-mmm'));
xtickangle(45);
ylabel('Cumulative imported infections');


% 只把均值线放进图例（中文 name11 不再使用）
% legend(h, name(sortIdx), 'Location','northeastoutside','FontSize',8);

% 7. 导出
export_file_name = fullfile(directory_name, ...
    sprintf('cumImport_ge1_11cities_CI_%s', currentDateTimeString));
saveas(f1, export_file_name+".fig");
exportgraphics(f1, export_file_name+".pdf", 'ContentType','vector');

%% ---------- 把过线城市名单写 Excel ----------
tbl1 = table(name(sortIdx), cumMean(idx1), ...
    'VariableNames', {'City','CumMean_day62'});

excelFile = fullfile(directory_name, ...
    sprintf('Cities_over_threshold_%s2.xlsx', currentDateTimeString));
writetable(tbl1, excelFile, 'Sheet', 'GE1')
%% Save data for plots
% save(strcat(currentDateTimeString,"_all_working_data.mat"))
% Save only essential statistics for large variables to reduce file size
% Initialize a structure to store statistics
stats_data = struct();

% Process import number variables (with cumulative sums)
% Process import number variables (with cumulative sums)
import_variable_pairs = {
    {'E_import_nums', 'AI_import_nums'};
    {'E_import_nums_owh', 'AI_import_nums_owh'};
    {'E_import_nums_twh', 'AI_import_nums_twh'};
    {'E_import_nums_s15', 'AI_import_nums_s15'};
    {'E_import_nums_s15_owh', 'AI_import_nums_s15_owh'};
    {'E_import_nums_l1w', 'AI_import_nums_l1w'};      % 新增
    {'E_import_nums_l2w', 'AI_import_nums_l2w'};      % 新增
    {'E_import_nums_s05', 'AI_import_nums_s05'};      % 新增
    {'E_import_nums_s05_owh', 'AI_import_nums_s05_owh'};
    {'E_import_nums_noIntv', 'AI_import_nums_noIntv'};% 新增
};
%  9×1 cell 数组，与行对应
import_scenario_names = {
    'Actual Intervention';...
    'One week early';...
    'Two weeks early';...
    'Strength x1.5';...
    'Strength x1.5 and one week early';...
    'One week late';...
    'Two weeks late';...
    'Strength x0.5';...
    'Strength x0.5 and one week early';...
    'Without Intervention'
};

% Calculate and store statistics for each import scenario

nDay = TIMELENGHT;
nSamp = num_iteration;

% 预分配
dayMean_I  = zeros(nDay, numel(import_variable_pairs));
dayLower_I = zeros(nDay, numel(import_variable_pairs));
dayUpper_I = zeros(nDay, numel(import_variable_pairs));
cumMean_I  = zeros(nDay, numel(import_variable_pairs));

for i = 1:numel(import_variable_pairs)
    % 先拼成一行，再 eval
    combined = eval([import_variable_pairs{i}{1} ' + ' import_variable_pairs{i}{2}]);         % 62×1000
    squeezed = squeeze(sum(combined,1));   % 把城市维压掉
    dayMean_I(:,i)  = mean(squeezed, 2);
    dayLower_I(:,i) = prctile(squeezed, 5, 2);
    dayUpper_I(:,i) = prctile(squeezed, 95, 2);
    cumMean_I(:,i)  = cumsum(dayMean_I(:,i));
end
% 拼表
tblI = array2table((1:nDay)', 'VariableNames', {'DayIndex'});
for i = 1:numel(import_variable_pairs)
    nm = strrep(import_scenario_names{i}, ' ', '_');
    tblI.([nm '_mean'])   = dayMean_I(:,i);
    tblI.([nm '_lower'])  = dayLower_I(:,i);
    tblI.([nm '_upper'])  = dayUpper_I(:,i);
    tblI.([nm '_cumsum']) = cumMean_I(:,i);
end

%% ---------- 2. 生成本地 NewInfections 统计 ---------- %%
infection_variables = { ...
    'NewInfections','NewInfections_noIntv', 'NewInfections_owh', 'NewInfections_twh', ...
    'NewInfections_s15', 'NewInfections_s15_owh', ...
    'NewInfections_l1w', 'NewInfections_l2w', ...
    'NewInfections_s05', 'NewInfections_s05_owh' ...
};
scenario_name = { ...
    'Fitting'; 'Without intervention';'One week early'; 'Two weeks early'; ...
    'Strength x1.5'; 'Strength x1.5 + one week early'; ...
    'One week late'; 'Two weeks late'; ...
    'Strength x0.5'; 'Strength x0.5 + one week early' ...
};

dayMean_N  = zeros(nDay, numel(infection_variables));
dayLower_N = zeros(nDay, numel(infection_variables));
dayUpper_N = zeros(nDay, numel(infection_variables));
cumMean_N  = zeros(nDay, numel(infection_variables));

for i = 1:numel(infection_variables)
    data = eval(infection_variables{i});   % 62×1000
    dayMean_N(:,i)  = mean(data, 2);
    dayLower_N(:,i) = prctile(data, 5, 2);
    dayUpper_N(:,i) = prctile(data, 95, 2);
    cumMean_N(:,i)  = cumsum(dayMean_N(:,i));
end

tblN = array2table((1:nDay)', 'VariableNames', {'DayIndex'});
for i = 1:numel(infection_variables)
    nm = strrep(scenario_name{i}, ' ', '_');
    tblN.([nm '_mean'])   = dayMean_N(:,i);
    tblN.([nm '_lower'])  = dayLower_N(:,i);
    tblN.([nm '_upper'])  = dayUpper_N(:,i);
    tblN.([nm '_cumsum']) = cumMean_N(:,i);
end

%% ---------- 3. 写 Excel 双 sheet ---------- %%
excelFile = fullfile(directory_name, ...
            sprintf('Import_and_Local_%s.xlsx', currentDateTimeString));
writetable(tblI, excelFile, 'Sheet', 'ImportRisk');
writetable(tblN, excelFile, 'Sheet', 'LocalInfection');
%%
function [dSv, dEv, dIv, dSh, dEh, dIh, dRh, Ihn] = SEI_SEIR_dev(t, Sv, Ev, Iv, Sh, Eh, Ih, Rh, TP, RF, CC, Nv, Nh, param)
    ie = gpv(param.ie);
    mu_h = gpv(param.mu_h);
    mu_n = gpv(param.mu_n);
    delta_h = gpv(param.delta_h);
    gamma_h = gpv(param.gamma_h);
    mu_v = 1./imu_v(TP, param);
    beta_v = gpv(param.beta_v);
    beta_h = gpv(param.beta_v);
    infection_rate_v = b(TP,param)*beta_v.*Ih./Nh;
    infection_rate_h = b(TP,param)*beta_h.*Iv./Nh;
    Nv = Sv+ Ev+ Iv;
    newVector = EFD(TP, param).*pEA(TP, param).*MDR(TP, param).*imu_v(TP, param).*(1-Nv./(CC*Nh))*Nv;
    newVector = max(newVector, 0);
    mu_v_increase =  param.mu_v_ratio;
    if t >= param.infection_rate_decline_begin1 &&  t <= param.infection_rate_decline_begin2
         mu_v_increase = param.mu_v_increase1;
    end
    if t >= param.infection_rate_decline_begin2
        mu_v_increase = param.mu_v_increase2;
    end
    mu_v = min(1-1e-8,mu_v_increase * mu_v);
    dSv = newVector - ...
        infection_rate_v.*Sv - mu_v.*Sv;
    dEv = infection_rate_v.*Sv - (gamma_v(TP, param)+mu_v) * Ev;
    gamma_v_tmp = param.gamma_v; %gamma_v(TP, param)
    % Ihn = gamma_v(TP, param).*Ev;
    dIv = gamma_v_tmp.*Ev-mu_v*Iv;
    dSh = (mu_n + ie)*Nh - infection_rate_h*Sh - (mu_h + ie)*Sh;
    dEh = infection_rate_h.*Sh - (gamma_h + mu_h + ie)*Eh;
    dIh = gamma_h*Eh - (delta_h + mu_h + ie)*Ih;
    Ihn = param.sym_ratio *gamma_h*Eh;
    dRh = delta_h*Ih - (mu_h + ie)*Rh;
end

function  [Sv, Ev, Iv, Sh, Eh, Ih, Rh, Ihn] = SEI_SEIR(t, Sv, Ev, Iv, Sh, Eh, Ih, Rh, TP, RF, CC, Nv, Nh, param)
    
    [dSv, dEv, dIv, dSh, dEh, dIh, dRh, Ihn] = SEI_SEIR_dev(t, Sv, Ev, Iv, Sh, Eh, Ih, Rh, TP, RF, CC, Nv, Nh, param);
    Sv = Sv + dSv; Sv = clip(Sv, 0, 2*Nh);
    Ev = Ev + dEv; Ev = clip(Ev, 0, 2*Nh);
    Iv = Iv + dIv; Iv = clip(Iv, 0, 2*Nh);
    Sh = Sh + dSh; Sh = clip(Sh, 0, Nh);
    Eh = Eh + dEh; Eh = clip(Eh, 0, Nh);
    Ih = Ih + dIh; Ih = clip(Ih, 0, Nh);
    Rh = Rh + dRh; Rh = clip(Rh, 0, Nh); 
end

function [MInfections, R0_array,ModelingOutput] = simulate(ps,dayOfYear,import_infection,daily_temperature,a14days_rainfall,carrying_capacity,param)
    timelengdth = 153;
    [Sv, Ev, Iv, Sh, Eh, Ih, Rh, Nv, Nh] = deal(ps*2, 0, 0, ps, 0, 0, 0, ps*2, ps);
    MInfections = zeros(timelengdth,1);
    Vpopulations = zeros(timelengdth,1);
    ModelingOutput = zeros(timelengdth,7);
    flag_begin = false;
    R0_array = zeros(timelengdth,1);
    EPS = 1e-4;
    for t = 1:timelengdth
        if t == dayOfYear
            Ih = import_infection;
            flag_begin = true;
        end

        if flag_begin
            TP = daily_temperature(t);
            RF = a14days_rainfall(t);
            CC = carrying_capacity(t)+EPS;
            
            [Sv, Ev, Iv, Sh, Eh, Ih, Rh, Ihn] = SEI_SEIR(t, Sv, Ev, Iv, Sh, Eh, Ih, Rh, TP, RF, CC, Nv, Nh, param);
            MInfections(t) = Ihn;
            Vpopulations(t) = Sv + Ev + Iv;
            ModelingOutput(t,:) =  [Sv, Ev, Iv, Sh, Eh, Ih, Rh];
            % R0_array(t) = cal_r0_v2(t, TP,  RF, CC, Nv, Nh, param);
        end
    
    end
end

function loglik_all = ssfun(local_param, data)
    param = data.param;
    dayofLocalObs = data.dayofLocalObs;
    ps_jiangmen = data.ps_jiangmen;        
    daily_temperature = data.daily_temperature;    
    a14days_rainfall = data.a14days_rainfall;     
    carrying_capacity = data.carrying_capacity;   
    local_infection = data.local_infection;
    smooth_rate = data.smooth_rate; 
    dayOfStartDate = data.dayOfStartDate;
    import_infection_c =   local_param(1);
    param.mu_v_increase1 = local_param(2);
    param.mu_v_increase2 = local_param(3);

    [NewInfection, R0_array, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c, ...
      daily_temperature,a14days_rainfall,carrying_capacity,param);
    
    loglik_all = sum((NewInfection(dayofLocalObs:dayofLocalObs+length(local_infection)-1) - local_infection').^2)/length(local_infection) / smooth_rate;
    
    % 14-19: sum 1664
    assert(dayofLocalObs==20,'dayofLocalObs changed, pls review the sum likelihood')
    sum_loglik = (sum(NewInfection(dayofLocalObs-6:dayofLocalObs-1)) - 1664).^2/5/ smooth_rate;
    loglik_all = loglik_all + 0.1 * sum_loglik;

    % important: is -loglik 
    assert(isscalar(loglik_all))
end

function ObsInfections = f_model(data, local_param)
    param = data.param;
    dayofLocalObs = data.dayofLocalObs;
    ps_jiangmen = data.ps_jiangmen;        
    daily_temperature = data.daily_temperature;    
    a14days_rainfall = data.a14days_rainfall;     
    carrying_capacity = data.carrying_capacity;   
    local_infection = data.local_infection;
    smooth_rate = data.smooth_rate; 
    dayOfStartDate = data.dayOfStartDate;

    import_infection_c =   local_param(1);
    param.mu_v_increase1 = local_param(2);
    param.mu_v_increase2 = local_param(3);

    [NewInfection, R0_array, ModelingOutput] = simulate(ps_jiangmen,dayOfStartDate,import_infection_c, ...
      daily_temperature,a14days_rainfall,carrying_capacity,param);
    ObsInfections = cat(2,NewInfection,ModelingOutput);
end