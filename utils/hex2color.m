function color = hex2color(hex_color)
    color = sscanf(hex_color(2:end),'%2x%2x%2x',[1 3])/255;
end

