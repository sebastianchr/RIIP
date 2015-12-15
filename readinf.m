function info=readinf(filename)
fid_INF=fopen(filename);

if fid_INF==-1; error('inf-file (%s) not found',filename); end

info.header=fgetl(fid_INF);                  %1
info.original_name=fgetl(fid_INF);           %2
info.IP_type=fgetl(fid_INF);                 %3
info.res1=str2double(fgetl(fid_INF));        %4
info.res2=str2double(fgetl(fid_INF));        %5
info.gradation=str2double(fgetl(fid_INF));   %6
info.pixelnum=str2double(fgetl(fid_INF));    %7
info.rasternum=str2double(fgetl(fid_INF));   %8
info.sensitivity=str2double(fgetl(fid_INF)); %9
info.latitude=str2double(fgetl(fid_INF));    %10
info.date_and_time1=fgetl(fid_INF);                  %11
info.date_and_time2=fgetl(fid_INF);                  %12
info.over_flow_pixels=str2double(fgetl(fid_INF));    %13
fclose(fid_INF);

end