function writegel(intensity,filename)
            ip_gel=(intensity*42948).^0.5;
            ip_gel=uint16(ip_gel); %unsigned 16bit integer
            
            
            numrows=size(ip_gel,1);
            numcols=size(ip_gel,2);
            
            fobj=Tiff(filename,'w');
            fobj.setTag('Photometric',Tiff.Photometric.MinIsWhite);
            fobj.setTag('Compression',Tiff.Compression.None);
            fobj.setTag('BitsPerSample',16);
            fobj.setTag('ResolutionUnit',Tiff.ResolutionUnit.Centimeter);
            fobj.setTag('SampleFormat',Tiff.SampleFormat.UInt);
            fobj.setTag('ImageLength',numrows);
            fobj.setTag('ImageWidth',numcols);
            fobj.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
            fobj.write(ip_gel);
            
            fobj.close;
            
end