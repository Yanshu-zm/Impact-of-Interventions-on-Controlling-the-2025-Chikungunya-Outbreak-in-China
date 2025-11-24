function tabel = f_readcsv(filename)
    opts = detectImportOptions(filename); 
    opts.VariableNamingRule = 'preserve';
    tabel = table2array(readtable(filename,opts));
end
