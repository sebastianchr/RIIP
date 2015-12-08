function intensity=readtif(filename)
            RawDataIP = imread(filename);
            intensity = (single(RawDataIP)+0.5)*(2^16-1)/42948;
end