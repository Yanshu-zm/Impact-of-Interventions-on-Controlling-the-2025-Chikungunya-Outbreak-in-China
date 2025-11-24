function v = indicator_after(current_t, threshold_t, before_value,after_value)
    if current_t < threshold_t
        v = before_value;
    else
        v = after_value;
    end
end

