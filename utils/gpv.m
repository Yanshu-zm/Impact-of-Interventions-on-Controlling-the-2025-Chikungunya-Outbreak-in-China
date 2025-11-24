% 获取参数值的统一函数get_param_value
function val = gpv(param_var)
    if isa(param_var, 'prob.NormalDistribution') || ...
       isa(param_var, 'prob.UniformDistribution')
        % 如果是概率分布对象，则进行采样
        val = random(param_var);
    elseif isa(param_var, 'function_handle')
        % 如果是函数句柄，调用函数获取样本
        val = param_var();
    else
        % 否则直接返回数值
        val = param_var;
    end
end