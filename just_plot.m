load("result/0803-佛山.mat");

if (true)
    % Stochastic simulation with MC
    if strcmp(source_city,'佛山')
        shenzhen_real_import_risk = [0	1	0	0	0	0	0	0	0	2	1	1	0	0];
        shenzhen_real_import_risk = [0	1	0	0	0	0	0	0 0,0,0,0,4,0,0,0,0];
    elseif strcmp(source_city,'广州')
        shenzhen_real_import_risk = [0	0	0	0	0	0	0	0	0	1	0	0	0	0];
    else
        shenzhen_real_import_risk = [0	0	0	0	0	0	0	0	0	0	0	0	0	0];
    end
    weekdayflow = GetImportFlow(strcat(source_city,'市'));
    import_num = zeros(52*7,num_iteration/10);
    import_num_tongbao = zeros(52*7,num_iteration/10);
    import_num_jiuzhen = zeros(52*7,num_iteration/10);
    if use_stochastic
        for k = 1 : num_iteration/10
            for t = dayOfYear:52*7
                % ModelingOutput(t,:) =  [Sv, Ev, Iv, Sh, Eh, Ih, Rh];
                % prevalance = ModelingOutputs(t,6,k*10)/ps;
                prevalance = ModelingOutputs(t,5,k*10)/ps; % 改用E
                cur_date = datetime(2025,1,0) + caldays(t);
                wd = weekday(cur_date);
                import_num(t,k) = immigrant_risk_ratio * normrnd(weekdayflow{wd,2},weekdayflow{wd,3}) * prevalance;
                import_num_tongbao(t,k) = immigrant_risk_ratio_tongbao * normrnd(weekdayflow{wd,2},weekdayflow{wd,3}) * prevalance;
                import_num_jiuzhen(t,k) = immigrant_risk_ratio_jiuzhen * normrnd(weekdayflow{wd,2},weekdayflow{wd,3}) * prevalance;

                % import_num(t,k) = weekdayflow{wd,2} * prevalance;
            end
        end
    else
        for t = 1:52*7
            % ModelingOutput(t,:) =  [Sv, Ev, Iv, Sh, Eh, Ih, Rh];
            prevalance = ModelingOutput(t,6)/ps;
            cur_date = datetime(2025,1,0) + caldays(t);
            wd = weekday(cur_date);
            import_num(t) = weekdayflow{wd,2} * prevalance;
        end
    end

    % sum(mean(import_num(dayOfYear:dayOfYear+14,:)'))
    % sum(mean(import_num(dayOfYear:dayOfYear+30,:)'))
    %%
    Names = {'Statistic', 'Mean', 'Median', '5%', '95%'};
    Items = {};
    Means = [];
    Medians = [];
    c5 = [];
    c95 = [];

    disp("23号至31号的预估数量")
    est_import= import_num(dayOfYear+days(datetime(2025,7,23)-startDate): ...
                        dayOfYear+days(datetime(2025,7,31)-startDate),:)';
    Items = [Items,{'省统计深圳本地病例'}];
    Means = [Means, sum(mean(est_import)) ];
    Medians = [Medians,sum(median(est_import)) ];
    c5 = [c5,sum(prctile(est_import,5)) ];
    c95 = [c95, sum(prctile(est_import,95)) ];
    
    est_import= import_num_tongbao(dayOfYear+days(datetime(2025,7,23)-startDate): ...
                        dayOfYear+days(datetime(2025,7,31)-startDate),:)';
    Items = [Items,{'通报协查病例'}];
    Means = [Means, sum(mean(est_import)) ];
    Medians = [Medians,sum(median(est_import)) ];
    c5 = [c5,sum(prctile(est_import,5)) ];
    c95 = [c95, sum(prctile(est_import,95)) ];

    est_import= import_num_jiuzhen(dayOfYear+days(datetime(2025,7,23)-startDate): ...
                        dayOfYear+days(datetime(2025,7,31)-startDate),:)';
    Items = [Items,{'深圳就诊'}];
    Means = [Means, sum(mean(est_import)) ];
    Medians = [Medians,sum(median(est_import)) ];
    c5 = [c5,sum(prctile(est_import,5)) ];
    c95 = [c95, sum(prctile(est_import,95)) ];

    % 创建表格
    T = table( ...
    Items', ...          % 转置为列向量
    Means', ...          % 转置为列向量
    Medians', ...        % 转置为列向量
    c5', ...             % 转置为列向量
    c95', ...            % 转置为列向量
    'VariableNames', Names ...  % 设置列名
    );
    
    % 显示表格
    disp(T);

    %%
    figure();
    subplot(2,1,1);

    hold on;
    title(source_city)
    % plot(MInfections(dayOfYear:end));
    CI_plot(mean(NewInfections(dayOfYear:end,:)'), ...
        prctile(NewInfections(dayOfYear:end,:)',5), ...
        prctile(NewInfections(dayOfYear:end,:)',95))
    plot(local_infection,'Color','r');
    % legend('Obs', 'Modeling');
    xlabel('Days');
    xlim([1,30]);
    ylabel('Incidence');

    subplot(2,1,2);
    hold on;
    boxObj = boxplot(import_num(dayOfYear:end,:)','symbol', ''); %'PlotStyle','compact',
    scatter((1:length(shenzhen_real_import_risk)),shenzhen_real_import_risk,'r*');
    xlabel('Days');
    ylabel('Number of imports');
    xlim([1,30]);
    ylim([0,5])

    tickDates = startDate + caldays((1:10:60));
    xticks((1:10:60))
    xticklabels(datestr(tickDates, 'mm-dd'));
    xtickangle(45);
    % save(strcat("result/0803-",source_city,'.mat'))
end