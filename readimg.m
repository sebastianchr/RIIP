function [intensity,info]=readimg(varargin)

if nargin ==1
    filename=varargin{1};
    disp('Looking for INF-file')
    filename_inf=regexprep(filename,'.img','.inf');
    info=readinf(filename_inf);
elseif nargin==2
    filename=varargin{1};
    info=varargin{2};
else
    error('Incorrect number of input arguments')
end



fid_IMG=fopen(filename,'r','b');  %Last parameter defines bit-reading order. 'b' is big endian. and 'l' is little endian i.e. Intel style.
if fid_IMG==-1; disp('inf file not found'); return; end
switch info.gradation
    case 8
        dat=fread(fid_IMG, 'uint8');
    case 16
        dat=fread(fid_IMG, 'uint16');
end
IPunits=reshape(dat, info.pixelnum,info.rasternum);
clear dat
fclose(fid_IMG);


%convert to intensity
%quantify data. formula from: 
%http://beamline.harima.riken.jp/bl45xu/web_old/Info/BAS2500imgSpec.pdf
%This page also contains som information on the file formats.

intensity=(info.res1/100)^2*4000/info.sensitivity*10.^(info.latitude*(IPunits/(2^info.gradation-1)-1/2));

end