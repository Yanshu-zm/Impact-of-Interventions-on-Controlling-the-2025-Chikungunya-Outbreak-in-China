if (true)
    [cities, FlowNet, pss] = GetFlowNetV2();
    ncites = length(cities);
    susceptible_cites = {'广东省佛山市','广东省广州市','广东省东莞市','广东省深圳市'};
    nsusc = length(susceptible_cites);
    indextmp = ismember(cities,susceptible_cites);
    inner_flow = FlowNet(indextmp,indextmp,:);
    outer_flow = FlowNet(indextmp,22:end,:);
    pss = pss(indextmp);

    % simulaton infections
    nation_import_nums = zeros(nsusc, ncites, 52*7); % 源头，目的地，时间
    asyminfection_nation_import_nums = zeros(nsusc, ncites, 52*7); % 源头，目的地，时间

    nation_models_outputs = zeros(nsusc, 52*7, size(ModelingOutput,2),num_iteration);
    nation_models_incidence = zeros(nsusc, 52*7,num_iteration);

    foshan_index = strcmp(susceptible_cites, "广东省佛山市");
    nation_models_outputs(foshan_index,:,:,:) = ModelingOutputs;
    
    f = waitbar(0, 'Starting');
    for iter = 1: num_iteration
        param.infection_rate_decline = param_chain(iter,2);
        waitbar(iter/num_iteration, f, sprintf('Progress: %d %%', floor(iter/num_iteration*100)));
        for t = dayOfYear-20:364
           prevalance = nation_models_outputs(:,t-1,5,iter)./pss; % 改用E
           prevalance_asymI = (1-param.sym_ratio) * nation_models_outputs(:,t-1,6,iter)./pss;

           wd = weekday(datetime(2025,1,0)+caldays(t));
           inner_flow_matrix = inner_flow(:,:,wd);
           outer_flow_matrix = outer_flow(:,:,wd);
           
           % trans in GD
           for n = (1:nsusc)
                if strcmp(susceptible_cites{n},'广东省佛山市')
                    continue
                end
                TP = daily_temperature(t);
                RF = a14days_rainfall(t);
                CC = carrying_capacity(t)+EPS;

                output_vector = squeeze(nation_models_outputs(n, t-1, :, iter))';
                Sv = output_vector(1);
                Ev = output_vector(2);
                Iv = output_vector(3);
                Sh = output_vector(4);
                Eh = output_vector(5);
                Ih = output_vector(6);
                Rh = output_vector(7);
    
                % Eh = Eh + sum(nation_import_nums(:,n,t,iter));
                % Ih = Ih + sum(asyminfection_nation_import_nums(:,n,t,iter));

                Eh = Eh + sum(nation_import_nums(:,n,t,iter));
                Ih = Ih + sum(asyminfection_nation_import_nums(:,n,t,iter));
                
                pstmp = pss(n);
                Nv = pstmp*2;Nh = pstmp;
                [Sv, Ev, Iv, Sh, Eh, Ih, Rh, Ihn] = SEI_SEIR(t, Sv, Ev, Iv, Sh, Eh, Ih, Rh, TP, RF, CC, Nv, Nh, param);
                
                nation_models_outputs(n,t,:,iter) = [Sv, Ev, Iv, Sh, Eh, Ih, Rh];
                nation_models_incidence(n,t,iter) = Ihn;
            
           end

                           nation_import_nums(:,n, t,iter) = immigrant_risk_ratio * flow_matrix(:,n) .* prevalance;
                asyminfection_nation_import_nums(:,n, t,iter) = immigrant_risk_ratio * flow_matrix(:,n) .* prevalance_asymI;


        end
    end
    close(f);
    result_import_table_E = table(cities,squeeze(sum(nation_import_nums,[1,3,4])/num_iteration)','VariableNames',{'输入城市','输入潜伏病例数'});
    result_import_table_asymI = table(cities,squeeze(sum(asyminfection_nation_import_nums,[1,3,4])/num_iteration)','VariableNames',{'输入城市','输入无症状病例数'});

    result_outbreak_sum = table(cities,squeeze(sum(nation_models_incidence,[2,3])/num_iteration),'VariableNames',{'城市','总计患病数'});
    
    begintmp = datenum(2025,7,28) - datenum(2025,1,0);
    endtmp = datenum(2025,8,4) - datenum(2025,1,0);
    result_outbreak_sum = table(cities,squeeze(sum(nation_models_incidence(:,begintmp:endtmp,:),[2,3])/num_iteration),'VariableNames',{'城市','总计患病数'});

    save(strcat('result/','startDay',num2str(start_day),'national-simulations-v2.mat'))
end