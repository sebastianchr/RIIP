function intensity=readgel(filename)
            RawDataIP = imread(filename);
            intensity = single(RawDataIP).^2/42948;
end