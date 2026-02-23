function delta = compute_delta(BI, A, tau)
    % COMPUTE_DELTA 从BI序列计算delta(t)
    %   delta = COMPUTE_DELTA(BI, A) 使用默认tau=3计算delta
    %   delta = COMPUTE_DELTA(BI, A, tau) 使用指定tau计算delta
    %
    % 输入:
    %   BI  - 一维数组，按时间顺序排列的BI值，允许NaN（缺失值）
    %   A   - 缩放参数，标量
    %   tau - 时间窗口参数，默认3
    %
    % 输出:
    %   delta - 一维数组，对应每个有效t的delta(t)值
    
    if nargin < 3
        tau = 3; % 按题目默认tau=3
    end
    if nargin < 2
        error('必须提供缩放参数A');
    end
    
    N = length(BI);
    t_start = tau + 1;          % 最小t，确保t-tau >=1
    t_end   = N - tau;          % 最大t，确保t+tau <= N
    
    if t_start > t_end
        error('BI序列长度不足！需要至少 %d 个数据点，当前只有 %d 个。', 2*tau + 1, N);
    end
    
    % 初始化输出
    num_t = t_end - t_start + 1;
    delta = zeros(num_t, 1);
    
    for t_idx = 1:num_t
        t = t_start + t_idx - 1; % 当前t在BI数组中的索引
        
        % 计算分子: sum_{i=1}^tau [sum_{k=t}^{t+i} BI(k)]
        numerator = 0;
        for i = 1:tau
            numerator = numerator + nansum(BI(t : t+i)); % 用nansum忽略NaN
        end
        
        % 计算分母: sum_{i=1}^tau [sum_{k=t-i}^{t-1} BI(k)]
        denominator = 0;
        for i = 1:tau
            denominator = denominator + nansum(BI(t - i : t - 1));
        end
        
        % 计算delta，避免分母为0
        if denominator == 0
            delta(t_idx) = NaN;
            warning('分母为0，t=%d处的delta设为NaN', t);
        else
            delta(t_idx) = -A * log(numerator / denominator);
        end
    end
end
