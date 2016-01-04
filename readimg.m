function [PSL,info]=readimg(varargin)

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
QL=reshape(dat, info.pixelnum,info.rasternum);
clear dat
fclose(fid_IMG);


%convert to intensity
%quantify data. formula from: 
%http://beamline.harima.riken.jp/bl45xu/web_old/Info/BAS2500imgSpec.pdf
%This page also contains som information on the file formats.

%PSL: Photo Stimulated luminescence
PSL=(info.res1/100)*(info.res2/100)*4000/info.sensitivity*10.^(info.latitude*(QL/(2^info.gradation-1)-1/2));

%to agree with data integrate from spring8 the date must be multiplied by
%100. 0.5 must also be subtracted. The 0.5 maybe comes because the PSL value
%represents the upper limit of a bin. By subtracting a half one obtains 
%something that is closer to the average value
PSL=PSL*100-0.5;


%EXCEPTION if QL=0 then PSL defines as PSL=0
PSL(QL==0)=0;

%The integration program at spring8 aparently does not remove these QL=0
%values. That can be seen at high angle where the intensity is 51 (= the number of pixels)

end