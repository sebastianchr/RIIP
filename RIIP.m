% function RIIP(varargin)
%Description: 
%Author: S. Christensen, Aarhus University (2015)
%Last edit: 10/11/2015 - 13:00
%How to use:
%
%
%
%Report all bugs and requests to sebastian@chem.au.dk

% Recent changes
% 06/11/15
% - Changes the default values for some settings to appropiate values.
% 05/11/15
% - checks that keywords only appear once in the input file
% 04/11/15
% - removed the separate "flip-data-step". Know there are buttons for rotating
% the data in the "define-image-plate-step". The data must rotated BEFORE
% selecting image plates.
% - In configuration file the keyword "flip" has been replaced with
% "rotate". "rotate" is an integere denoting how many time 90 degrees the
% data should be rotated counterclockwise, e.g 0 and 4 correspondes to no
% change.
%- In "define-image-plate-step" added keyboard short cuts. Find them in the
%menu-point "Hot Keys". When the zoom tool is selected. To zoom out, hold
%down shift and mouseclick or right click to select the desired action
%2/11/15
% - added button on plot of integrated data to save settings and integrated
% data
%1/11/15
%- Fit of beam center to curved path uses either circle or debye cone
%depending on the choice of 'path shape'
%- To fit the beam center from curve click somewhere in the reflection.
%10/11/15
%-mask certain number of pixels long strip along the edges
%08/12/15
%-implement different pixelsizes in X and Z-direction. NEEDS CHECKING
%especially for debye cone 
% - implemented difference in pixelsize in integration weighting 
% (by length and area)


%todo:
%I think there might be a problem with the relative beam center if the
%long side of the image plate is turning downwards.
%Error in the integration step for certain lengths of imageplates. This bug
%is however difficult to reproduce
%Automatic mask edges
%

%Fixed: 
% width of imageplate is now update when users changes value of the
% belonging field
% Fixed proplems reading 'line' and 'rectangle' from configuration file 
% In IP-select: the selected borders (blue polygon) was incorrectly read
% from the positions of the points instead of the red polygon.
% re-added the 2theta axis to the post-integration 2d-analysis


function RIIP(varargin)
if nargin ==0
    % No input is given. The program enters interactive mode where most
    % settings can be set interactively.
    interactive_mode=2 ;
    % set 'interactive_mode' to 2 for all options and to 1 for limited 
    % options
elseif nargin==1
    %Settings are read from a configuration file
    settings0=ascii2struct(varargin{1});
    interactive_mode=0;
    % If an essential setting is not contained in the configuration file
    % and no default value exist then the user is prompted for input
else
    % The settings can also be input directly into the function 
    settings0=struct(varargin{:});
    interactive_mode=0;
end


% set(0,'DefaultAxesFontSize', 16);
scrsz = get(0,'ScreenSize');
%surpress a warning that is issued when the tiff file is read. Someone
%should check if this has any consequences
warning('off', 'MATLAB:imagesci:rtifc:notPhotoTransformed')

%% Default settings
default=[];
% Instrument settings
default.pixel_size=[24.941e-3/0.9926 24.941e-3];%mm
default.camera_radius=1200.5; %mm
default.conversion_factor=4.5;


%Data file settings
default.filename='';
default.decay_constants='new'; %options: none, new, old, (a0 t1 t2) In the last option the user provides the numeric paramters for a double exponential decay t1 and t2 are in minutes.
default.exposure_time=60;
default.waiting_time=3.5;
default.reading_time=3.9;
default.scanner_correction='nanna'; %options: none, nanna, kasper
default.rotate=[];


default.beam_center_guess=[];

%Image plate specific settings
default.image_plate_borders=[];  
default.image_plate_beam_center=[];  
default.unwarp_image_plate='oneline'; %options: none, parallelogram, tetragon, twolines, onelines
default.mask_polygon=[];
default.mask_inout=[];  

% Integration settings
default.path_shape='debye_cone';  %options: line, circle, debye_cone
default.weighting_method='area';  %options: unity, length, area


%Save-settings
default.save_settings= [];
default.save_result='$';  %0: do not save,
                          %$: create output file name from data filename,
                              %any other string: save as that
%Visual/info settings
default.resolution_reduction=5;
default.verbose=2;      %lower value means less information while running
%0: no plots
%1: mainly 1D plots
%2: also 2D plots of masked data and possibility for analysis
%3: everything

if interactive_mode>0;
    %launch user interface
    settings=default;
    settings.filename='';
    settings.save_result='';
    settings.save_settings='';
    
    if interactive_mode>1
        %function to input instrument parameter
        
        % ######### This part you can comment and uncomment ###############
        settings=input_instrument_parameters(settings);
        %settings = input_integration_parameters(settings);
        %settings = input_position_correction_method(settings);
        % #################################################################
        
        %settings = input_unwarp_method(settings);
    end
    
else
    % As standard the default values are used
    settings=default;
    % If a parameter is defined in the settings-file (as loaded into 
    % settings0), then this value is used instead
    if (isfield(settings0,'pixel_size'));            settings.pixel_size=settings0.pixel_size;  end  %in mm
    if (isfield(settings0,'camera_radius'));         settings.camera_radius=settings0.camera_radius;  end  %in mm
    if (isfield(settings0,'conversion_factor'));     settings.conversion_factor=settings0.conversion_factor;  end
    if (isfield(settings0,'image_plate_borders'));   settings.image_plate_borders=settings0.image_plate_borders;   end
    if (isfield(settings0,'image_plate_beam_center'));settings.image_plate_beam_center=settings0.image_plate_beam_center;   end
    if (isfield(settings0,'rotate'));                settings.rotate=settings0.rotate;   end
    if (isfield(settings0,'resolution_reduction'));  settings.resolution_reduction=settings0.resolution_reduction;   end
    if (isfield(settings0,'path_shape'));            settings.path_shape=settings0.path_shape;   end
    if (isfield(settings0,'weighting_method'));      settings.weighting_method=settings0.weighting_method;   end
    if (isfield(settings0,'filename'));              settings.filename=settings0.filename;  end
    if (isfield(settings0,'beam_center_guess'));     settings.beam_center_guess=settings0.beam_center_guess;   end
    if (isfield(settings0,'mask_polygon'));          settings.mask_polygon=settings0.mask_polygon;  end
    if (isfield(settings0,'mask_inout'));            settings.mask_inout=settings0.mask_inout;  end
    if (isfield(settings0,'decay_constants'));       settings.decay_constants=settings0.decay_constants;  end
    if (isfield(settings0,'exposure_time'));         settings.exposure_time=settings0.exposure_time;  end
    if (isfield(settings0,'waiting_time'));          settings.waiting_time=settings0.waiting_time;  end
    if (isfield(settings0,'reading_time'));          settings.reading_time=settings0.reading_time;  end
    if (isfield(settings0,'scanner_correction'));    settings.scanner_correction=settings0.scanner_correction;  end
    if (isfield(settings0,'unwarp_image_plate'));    settings.unwarp_image_plate=settings0.unwarp_image_plate;  end
    if (isfield(settings0,'save_settings'));         settings.save_settings=settings0.save_settings;  end
    if (isfield(settings0,'save_result'));           settings.save_result=settings0.save_result;  end
    if (isfield(settings0,'verbose'));               settings.verbose=settings0.verbose;  end

end

standard_colors=[0 0 1; 0 0.5 0; 1 0 0 ; 0 0.75 0.75; ...
                 0.75 0 0.75; 0.75 0.75 0 ; 0.25 0.25 0.25];

%i: refers to data file
%j: refers to imageplate of datafile i

global precision verbose dfs pixsize R resred
precision = 'single';
verbose =settings.verbose;
dfs=scrsz(3:4)*0.70;
if numel(settings.pixel_size)==1
    pixsize=repmat(settings.pixel_size,1,2);
elseif numel(settings.pixel_size)==2
    pixsize=settings.pixel_size;
else
    error('Incorrect number of elements in pixel_size')
end
R=settings.camera_radius;
resred=settings.resolution_reduction;

%ensure that certain settings are in a cell-array if there is only a single
%value;
if ischar(settings.filename)
    temp={}; temp{1}=settings.filename;
    settings.filename=temp;
end



%% Read raw data. 
% Data is a cell-array. Each entry is a saparate data file as
% scanned by the scanner. Settings is updated to contain the filenames of
% data files if these are manually selected by the user during the call of
% "read_raw_data"
%The data files must be selected from low to high angle.

if verbose >0
    fprintf('Starting to load data\n')
end
[data, settings]=read_raw_data(settings);
number_of_data_files=length(settings.filename);

if verbose >0
    fprintf('Finished loading data\n')
end

%% Apply corrections to the scanned data
if verbose >0
    fprintf(['Starting to correct imageplate data for intensity decay and for '...
        'inaccurate scanner position\n'])
end
%Decay correction. Corrects for the decay of the imageplate while scanning
data=correct_for_decay(data,settings);

%Correct the scanner position. The position recorded by the scanner is 
% slightly incorrect in the Z-direction. THe position error-appear to be 
% stable between scans. Therefore the Z-axis is corrected. Data array 
% returned below is interpolated to correspond to the true Z-position 
if ~strcmpi(settings.scanner_correction,'none')
    for i=1:number_of_data_files
        data{i}=correct_scanner_position(data{i},...
                                  settings.scanner_correction);
    end
end

%% Rotate the data
% Scanned IP must be oriented so that lowest 2theta corresponds to low
% Z-pixel-value. If not the case the data can be rotated in the 
% for each data file, an entry in settings.rotate indicate if data should be
% left as-is (=0) or rotated counter-clockwise by 90*settings.rotate(i)
% degrees

if ~isempty(settings.rotate)
     %change the orientation of the raw data
     if length(settings.rotate)>=number_of_data_files
         for i=1:number_of_data_files
             data{i}=rot90(data{i},settings.rotate(i));
         end
     else
         error('Keyword ''rotate'' not defined for all scans')
     end
end

if verbose >0
    fprintf('Finished making scanner corrections\n')
end

%% Determine the position of image plates
% The position of each IP is defined by 4 points.
% if the positions are not defined in  settings.image_plate_borders then a
% userinterface is launched. 
ip_positions=cell(number_of_data_files,1);
data_ip=cell(number_of_data_files,1);
number_of_image_plates=zeros(number_of_data_files,1);
ip_borders_cut=cell(number_of_data_files,1);
ip_length=cell(number_of_data_files,1);

if isempty(settings.image_plate_borders)
    %The user selects the position of imageplates graphically. This is
    %important to do accurately since this is used to define the coordinate
    %system for the integration.
    for i=1:number_of_data_files
        [data{i}, ip_positions{i},settings.rotate(i)]=define_image_plate_graphic(data{i},settings);
        %image plates must be selected in the order from low angle to high
        %angle. This cannot be changed later
    end
    settings.image_plate_borders=ip_positions;
end

%count number of image plates pr data file
for i=1:number_of_data_files
    number_of_image_plates(i)=length(settings.image_plate_borders{i});
end
% Check that image plates have been selected before continuing or throw an
% error. Only the first datafile is checked, since it must contain at least
% one image plate which contains the beam center.
if number_of_image_plates(1)==0
    error('No image plates selected for the first data file')
end

%Determine the length of each image plate from the coordinates of the IP
%borders
for i=1:number_of_data_files
    ip_length{i}=zeros(number_of_image_plates(i),1);
    for j=1:number_of_image_plates(i)
        [~, ip_length{i}(j,1)]=get_ip_origin_and_length(...
            settings.image_plate_borders{i}{j});
    end
end


%% Optimize IP-borders
% At this point, I could make a function to optimice the image plate 
% borders

%% Cut out individual image plates
% The full data arrays are divided into smaller sections. Each ection
% containing only a single image_plate
for i=1:number_of_data_files
    [data_ip{i}, ip_borders_cut{i}]=cut_raw_data(data{i}, ...
        settings.image_plate_borders{i});
    % data_ip is nested cell-array. 1st level refers to each raw-data-file.
    % The second layer contains all the separate image plates for each data
    % file.
end
clear data %clean up to save memory
%% Unwarp the recorded image plate
% If the image plate was misaligned in the scanner, this function can rotate
% the image plate so it aligned with the integration coordinate system. 
if ~strcmpi(settings.unwarp_image_plate,'none')
    for i=1:number_of_data_files
        for j=1:number_of_image_plates(i)
%           also output unwarped ip_borders
            [data_ip{i}{j}, ip_borders_cut{i}{j}]= ...
                transform_to_rectangular_grid(data_ip{i}{j}, ...
                ip_borders_cut{i}{j},settings.unwarp_image_plate);
        end
    end
end

%% Find the center of the beam
%If data files and image plates have been selcted in the correct order,
%the beam center should always be in the first data file and on the first 
%image plate.
i_bc=1; j_bc=1;

%checks if a guess (or refined value) for the beam center has been defined.
%if not the user will be asked to select the position of the beam center

bc_ip=settings.image_plate_beam_center;
if isempty(bc_ip)
    try
        bc_guess=settings.beam_center_guess;
    catch err   %#ok
        bc_guess=[];
    end
    if isempty(bc_guess)
        [xguess, zguess]=...
            get_beam_center_from_user_graphical(data_ip{i_bc}{j_bc});
        bc_guess=[xguess, zguess];
    end
    
    %Refine the position of the beam center
    bc_refined=refine_beam_center(data_ip{i_bc}{j_bc},bc_guess);
    
    bc_ip=cell(number_of_data_files,1); % the beam center relative to the cut-out IP
    off_set_temp=0;
    for i=1:number_of_data_files
        for j=1:number_of_image_plates(i)
            bc_ip{i}(j,:)=bc_refined-[0 off_set_temp];
            off_set_temp=off_set_temp+ip_length{i}(j)+1;
        end
    end
else
    if verbose >0
        disp('Will use beam center as defined in input')
    end
end

settings.image_plate_beam_center=bc_ip;

%% Mask the data
% First the masks are a number of polygons (mask_polygon) and a value 
% (mask_inout). The later dictates whether the data inside or outside the
% polygon is kept. These polygons are converted into an actual mask of same
% dimensions as data_ip{i}{j}. It is an array of 1 and 0 where 1 denotes
% that a particular pixel should be kept.
if isempty(settings.mask_polygon)
    % In this step the user can define a mask for each image plate. The mask
    % will determine which part of the data will be integrated. Each mask is
    % defined by a polygon ("mask_poly{i}{j}{k}") and a scalar
    % "inout{i}{j}(k)". This scalar can be 0 or 1 and determines the area
    % inside or outside the polygon is to be excluded from the integration.
    mask_poly=cell(number_of_data_files,1);
    inout    =cell(number_of_data_files,1);
    mask     =cell(number_of_data_files,1);
    if interactive_mode==1
       for i=1:number_of_data_files
            for j=1:number_of_image_plates(i)
                %user graphically selects the desired masks.
                mask{i}{j}=ones(size(data_ip{i}{j}));
            end
        end
    else
        for i=1:number_of_data_files
            for j=1:number_of_image_plates(i)
                %user graphically selects the desired masks.
                [mask_poly{i}{j}, ...
                    inout{i}{j}, ...
                    mask{i}{j}]=define_mask(data_ip{i}{j},ip_borders_cut{i}{j},settings);
            end
        end
    end
    settings.mask_polygon=mask_poly;
    settings.mask_inout=inout;
else
    %     If the coordinates for mask-polygons are given in the input-files,
    %     will these be converted into the actual mask with same
    %     dimensions as data_ip{i}{j}
    
    mask=cell(number_of_data_files,1);
    for i=1:number_of_data_files
        mask{i}=cell(number_of_image_plates(i),1);
        for j=1:number_of_image_plates(i)
            mask{i}{j}=define_mask_from_polygons(data_ip{i}{j},settings.mask_polygon{i}{j},settings.mask_inout{i}{j});
        end
    end
end

% Save the settings to an ascii file
save_settings_to_file(settings,true)

%% Do data integration
% The 2D data is integrated to produce a 1D diffraction pattern for each
% image plate. The resulting 1D-patterns are normalized to the integration
% path, and the scanner correction is applied to the pixel-positions.
integrated_data=integrate_data(data_ip,mask,settings);
% "integrated_data" is a nested cell-array containing a file and
% image-plate level

if verbose > 0
    %Plot all diffractograms
    h=cell(number_of_data_files,1);
    h_esd=cell(number_of_data_files,1);
    fig_int=figure('Outerposition',[ 0  50   dfs],...
        'Name','Integrated data vs 2theta',...
        'NumberTitle','off');
    
    
    ax_int=axes('Parent',fig_int,...
        'units','normalized',...
        'position',[0.1 0.15 0.85 0.80 ]);
    hold(ax_int,'on')
    
    for i=1:number_of_data_files
        h{i}=zeros(number_of_image_plates(i),1);
        h_esd{i}=zeros(number_of_image_plates(i),2);
        for j=1:number_of_image_plates(i)
            x=integrated_data{i}{j}(:,1);
            y=integrated_data{i}{j}(:,2);
            e=integrated_data{i}{j}(:,3);
            h{i}(j)=plot(ax_int,x,y,'-k');
            %plot upper and lower limit
            h_esd{i}(j,1:2)=plot(ax_int,[x, x], [y-e, y+e],'-',...
                'color',standard_colors(sum(number_of_image_plates(1:i-1))+j,:));
            % can also plot zero-width error ticks with
            % plot([x, x]', [y-e, y+e]', '-k');
            %however this is very slow. Probably because of large
            %overhead when each errorbar gets its own handle
            clear x y e
        end
    end
    xlabel('2\theta [degree]')
    ylabel('Integrated intensity [photons]')
    set(gca,'fontsize',16)

    
    h1 = uibuttongroup('visible','on','units','pixels',...
    'Position',[5 5 230 40]);
h11 = uicontrol('Position',[10 10 100 20],...
    'parent',h1,'String','Save data',...
    'Callback',@(src,evt) save_integrated_data(integrated_data,settings,false)); %#ok
h12 = uicontrol('Position',[120 10 100 20],...
    'parent',h1,'String','Save settings',...
    'Callback',@(src,evt) save_settings_to_file(settings,false)); %#ok
    
    
end

%% Save
% Save the integrated data. Each diffractogram is saved to a separate file
save_integrated_data(integrated_data,settings,true)


end  %This concludes the main program


%% Sub functions
%Below are the sub-functions defined.
%% Settings userselection 
function settings=input_instrument_parameters(settings)
%Get instrument settings. If more than one data-file has been loaded, these
% parameters are assumed to the equal for all of them.
prompt = {'Camera radius (mm)' ...
    'Pixel size X (mm)' ...
    'Pixel size Z (mm)' ...
    'Conversion factor'};

name = 'Instrument corrections';
numlines = 1;
defaultanswer = {num2str(settings.camera_radius), ...
    num2str(settings.pixel_size(1)), ...
    num2str(settings.pixel_size(2)), ...
    num2str(settings.conversion_factor)};
options.Resize='on';
options.WindowStyle='normal';
answer=inputdlg(prompt,name,numlines,...
    defaultanswer,options);

if isempty(answer)
    %check if "Cancel" was pressed. If so, issue this error and end program
    error('User did not define mandatory settings. Program ended.')
end

% The inputs are given ads strings. COnvert here to numbers.
for ii = 1:length(answer)
    answer{ii} = str2double(answer{ii});
end

%Save parameters fot the settings.structure
settings.camera_radius=answer{1};
settings.pixel_size(1)=answer{2};
settings.pixel_size(2)=answer{3};
settings.conversion_factor=answer{4};
end


function [settings] = input_integration_parameters(settings)
%Choose methods for integration settings. If more than one data-file has
% been loaded, these parameters are assumed to the equal for all of them.

str = {'Debye cone','Circle', 'Line'};
str_code = {'debye_cone','circle','line'};
[s,v] = listdlg('PromptString',...
    'Select shape of integration path across imageplate',...
    'SelectionMode','single',...
    'ListString',str,...
    'Listsize',[160 80]);
if v
    settings.path_shape=str_code{s};
else
    error('Did not choose a shape for integration path')
end

str = {'area','length', 'unity'};
str_code = str;
[s,v] = listdlg('PromptString','Select weighting method for integration',...
    'SelectionMode','single',...
    'ListString',str,...
    'Listsize',[160 80]);
if v
    settings.weighting_method=str_code{s};
else
    error('Did not choose a shape for border')
end

settings.save_settings=[];
end

function [settings] = input_unwarp_method(settings)
%Choose methods for integration settings. If more than one data-file has
% been loaded, these parameters are assumed to the equal for all of them.

str = {'Tetragon','Parallelogram', 'None'};
str_code = {'tetragon','parallelogram','none'};
[s,v] = listdlg('PromptString',...
    'Select a shape for unwarping the scanned image plate',...
    'SelectionMode','single',...
    'ListString',str,...
    'Listsize',[160 80]);
if v
    settings.unwarp_image_plate=str_code{s};
else
    error('Did not select a shape for unwarping')
end
end

function [settings] = input_position_correction_method(settings)
%Choose methods for integration settings. If more than one data-file has
% been loaded, these parameters are assumed to the equal for all of them.

str = {'Nanna - Poly.','Kasper -Spline','None',};
str_code = {'nanna','kasper','none'};
[s,v] = listdlg('PromptString',...
    'Select function for correcting the scanner position',...
    'SelectionMode','single',...
    'ListString',str,...
    'Listsize',[160 80]);
if v
    settings.scanner_correction=str_code{s};
else
    error('Did not choose a function for scanner position correction')
end
end

function [exposure_time, ...
    waiting_time, ...
    reading_time]=input_scanner_parameters(default)
% Provide integration parameters for decay-correction. 
% A set of these parameters must be inputed for each data-file.

prompt = {'Enter exposure time (min)',...
    'Enter waiting time (min)',...
    'Enter reading time (min)'};

defaultanswer = {num2str(default.exposure_time),...
                 num2str(default.waiting_time),...
                 num2str(default.reading_time)};

options=struct('Resize','on','WindowStyle','normal');

answer=inputdlg(prompt,'Instrument corrections',1,defaultanswer,options);

if isempty(answer)
    error('User did not define mandatory settings. Program ended.')
end

for ii = 1:length(answer)
    answer{ii} = str2double(answer{ii});
end

exposure_time=answer{1};
waiting_time=answer{2};
reading_time=answer{3};

end

%% Read data from file 
function [DataIP, settings]=read_raw_data(settings)
% Load the raw 2D data. If no filenames are contained in settings.filename
% the user will be prompted to select the desired files. Multiple files can
% be selected. The files must be selected in order of low to high angle.
% For each loaded file the user is prompted to input scanning time
% parameters.
% If the filenames are defined in settings.filename, the user is not prompted
% it is assumed that scanner parameters have been defined correctly. If not
% it will give an error somewhere else. FIXLATER.
RawFileName=settings.filename;
global verbose
if isempty(RawFileName{1})
    RawFileName={};
    questdlg(['You must select the data files in order of low angle '...
        'to high angle'],'-','OK','OK');
    for i=1:10
        %Open GUI to select the tif file (or dat)
        [tempFileName,tempPathName] = uigetfile({'*.gel';'*.tif';'*.dat'},...
            'Select Image Plate File (tif/gel/dat)');
        
        % End script execution cleanly if the user press cancel
        if tempFileName == 0
            if i==0
                error('You pressed cancel. No file selected.')
            else
                %Use the filenames that have already been selected
                break
            end
        end
        RawFileName{i}=[tempPathName,tempFileName];   %#ok
        
        %user input scanner parameters
        [t_exp, t_wait, t_read]=input_scanner_parameters(settings);
        settings.exposure_time(i)=t_exp;
        settings.waiting_time(i) =t_wait;
        settings.reading_time(i)=t_read;
        
        choice = questdlg('Load more data files?',...
            'Load more','yes','no','no');
        if strcmp(choice,'no')
            break
        end
    end
    
end
settings.filename=RawFileName;
DataIP=cell(1,length(RawFileName));
for i=1:length(RawFileName)
    [~,~,file_format] = fileparts(RawFileName{i});
    switch file_format
        case '.tif'
            % Open the selected file
            %         RawDataInfo = imfinfo(RawFileName) ;
            
            RawDataIP = imread(RawFileName{i});  
            DataIP{i} = (convert2floatingpoint(RawDataIP)+0.5)*(2^16-1)/42948;
        case '.gel'
            RawDataIP = imread(RawFileName{i});
            DataIP{i} = convert2floatingpoint(RawDataIP).^2/42948;
        case '.dat'
            DataIP{i} = load(RawFileName{i});
        otherwise
            error('Fileformat not recognized')
    end
end
end


%% Decay correction
function int=correct_for_decay(int,settings)
global verbose
number_of_data_files=length(int);
% Set the constants for the decay function which defines the expected
% intensity decay as function of time.
if ischar(settings.decay_constants)
    switch settings.decay_constants
        case 'old'
            a0=0.36;
            t1=14.9; t2=434; %[min]
        case 'new'
            a0=0.3827;
            t1=8.3986; t2=324; %[min]
        case 'none'
            a0=1;
            t1=inf; t2=inf; %[min]
        otherwise
            error('Unknown string in ''decay_constants''')
    end
elseif isnumeric(settings.decay_constants)
    decay_par=settings.decay_constants;
    a0=decay_par(1);
    t1=decay_par(2);
    t2=decay_par(3);
else
    error('In decay correction: unknown type for "decay_constants"')
end
%the decay function as an
decay = @(t) a0.*exp(-t./t1)+(1-a0).*exp(-t./t2);

if verbose >2
    fig_decay=figure('Name','Decay Correction','NumberTitle','off'); hold on;
    xlabel('Time (min)')
    ylabel('Decay fraction')
end
for i=1:number_of_data_files
    %     convert the z-axis coordinates into a time scale
    tstart=settings.exposure_time(i)/2+...
        settings.waiting_time(i);
    tend=settings.exposure_time(i)/2+...
        settings.waiting_time(i)+...
        settings.reading_time(i); %time from start expoure to end development
    t=linspace(tstart,tend,size(int{i},2));
    
    int{i}=int{i}./repmat(decay(t),size(int{i},1),1);
    
    %Plot the results
    if verbose >2
        figure(fig_decay);
        plot(t,decay(t));
    end
end
end

%% Scanner position correction
function data_cor=correct_scanner_position(data,correction_function)

switch correction_function
    case 'nanna'
        pp1=load([ getpath '/spline_s3_px24k941.mat']);
        x_rec2true=  @(x) x; %#ok % from true value to recorded value    
        z_rec2true=  @(z) z-ppval(pp1,z);
        % this is approximate at best since the true recorded value should
        % go into ppval
        x_true2rec=  @(x) x; %  recorded value to true value
        z_true2rec=  @(z) z+ppval(pp1,z);
        
    case 'kasper'
        %cor.f will always return a collumn vector. Therefore z-must be a
        %collumn-vector
        cor=load([ getpath '/scanner_cor_spline_5em9.mat']);
        x_rec2true=  @(x) x; %#ok % from true value to recorded value
        z_rec2true=  @(z) z+cor.f(z);
        % this is approximate at best since the true recorded value should
        % go into the spline function f
        x_true2rec=  @(x) x; %  recorded value to true value
        z_true2rec=  @(z) z-cor.f(z);
    otherwise
        error('Unknown value for ''scanner_correction''')
end

nx=size(data,1); nz=size(data,2);
method=2;

%method 2 can be used if there is ever a need also to correct the position
% in X-direction.

switch method
    case 1
        z_rec=(1:nz)';
        z_true=z_rec2true(z_rec);
        data_cor=zeros(size(data));
        for i=1:nx
            data_cor(i,:)=interp1(z_true,data(i,:),z_rec);
        end
        data_cor=convert2floatingpoint(data_cor);
    case 2
                x_true=(1:nx)';  z_true=(1:nz)';
                x_rec=x_true2rec(x_true);
                z_rec=z_true2rec(z_true);
                [Z_rec, X_rec]=meshgrid(z_rec,x_rec);
                data_cor=convert2floatingpoint(interp2(data,Z_rec,X_rec));
% This method is not qute accurate but for 2D interpolation it is much, 
% much faster because the grid is regular. The two methods give almos 
% identical results
                
end

    function path = getpath()
        %get the path of this m-file
        [path , ~ , ~] = fileparts(mfilename('fullpath'));
    end
end

function [data, ip_position,rotation]=define_image_plate_graphic(data, settings)
global dfs
%this function works currently but the use of callbacks is a bit messy
% In this function the user selects the position of individual image image
% plates, from the scanned image.

% The shape of the image plate can be adjusted by dragging 8 points. Their
% number refers to their handle in variable h. The <,>,^,v indicate the
% directions the points can be dragged, i.e. 1,2,3,4,6,8 moves in the
% direction of the lines 4-3 and 1-2. The entire figure can be rotated
% around 4 by holding down "Alt" and draggind 2 or 3. 
%                 ^
% <4>-------------7------------ <3> (Alt-> rotate)
%  |              v             /
%  |                           /
% <8>                        <6>
%  |                         /
%  |           ^            /
% <1> ---------5----------<2>  (Alt-> rotate)
%              v

nx=size(data,1); nz=size(data,2);
x_limit=[1 nx];
z_limit=[1 nz];

%These functions are the output. Values are added to them in the
%subfunctions
ip_position=[];
rotation=0;

% Define the figure and user-interface
fig2d=figure('Outerposition',[ 0  50 dfs],  ...
    'Name','Select image plates',...
    'NumberTitle','off',...
    'Toolbar','none',...
    'menubar','none');

% Add menus with keyboard shortcuts ('Accelerators')
mymenu = uimenu('Parent',fig2d,'Label','Hot Keys');
uimenu('Parent',mymenu,'Label','Zoom on','Accelerator','a','Callback',@(src,evt)zoom(fig2d,'on'));
uimenu('Parent',mymenu,'Label','Zoom off','Accelerator','s','Callback',@(src,evt)zoom(fig2d,'off'));
uimenu('Parent',mymenu,'Label','Undo','Accelerator','d','Callback',@undo_ip);
uimenu('Parent',mymenu,'Label','Finish','Accelerator','f','Callback',@finish);
uimenu('Parent',mymenu,'Label','New IP','Accelerator','n','Callback',@select_ip_ui);


ax2d = axes('Parent',fig2d,...
    'units','normalized',...
    'position',[0.135  0.20 0.80   0.75]);


h1 = uibuttongroup('visible','on','units','pixels',...
    'Position',[5 5 360 45]);
h11 = uicontrol('Position',[10 10 50 20],...
    'parent',h1,'String','New IP',...
    'Callback',@select_ip_ui); %#ok
h12 = uicontrol('Position',[160 10 50 20],...
    'parent',h1,'String','Finish',...
    'Callback',@finish); %#ok
%Display the length of the IP. If edited; the drawn rectangle will be
%updated to march the new value
h14 = uicontrol('style','text','position',[240 25 50 20],...
    'parent',h1,'String','IP length'); %#ok
h13 = uicontrol('style','edit','position',[240 10 50 20],...
    'parent',h1,'String','0',...
    'callback',@update_length);

%Display the width of the IP. If edited; the drawn rectangle will be
%updated to march the new value
h16 = uicontrol('style','text','position',[300 25 50 20],...
    'parent',h1,'String','IP width'); %#ok
h15 = uicontrol('style','edit','position',[300 10 50 20],...
    'parent',h1,'String','0',...
    'callback',@update_width);

h17 = uicontrol('position',[80 10 50 20],...
    'parent',h1,'String','Undo IP','callback',@undo_ip); %#ok

%Create fields to set lower and upper limit for color-axis
[h_clim_min, h_clim_max]=add_colorlimit_fields(fig2d,ax2d,[385 5 120 45]);

h2 = uibuttongroup('visible','on','units','pixels',...
    'Position',[515 5 180 45]);
h20 = uicontrol('style','text','position',[40 25 100 20],...
    'parent',h2,'String','Rotate data'); %#ok
h21 = uicontrol('Position',[5 10 50 20],...
    'parent',h2,'String','+90',...
    'Callback',@(src,evt) rotate_data(1)); %#ok
h22 = uicontrol('Position',[65 10 50 20],...
    'parent',h2,'String','-90',...
    'Callback',@(src,evt) rotate_data(-1)); %#ok
h23 = uicontrol('Position',[125 10 50 20],...
    'parent',h2,'String','180',...
    'Callback',@(src,evt) rotate_data(2)); %#ok

    function rotate_data(k_rotation)
        data = rot90(data,k_rotation);
        rotation=rotation+k_rotation;
        p_2d=plot_2d_data(p_2d,data);
    end

    function p_handle=plot_2d_data(varargin)
        if nargin==1
            data0=varargin{1};
            p_handle=imagesc_highspeed(log10(data0),settings.resolution_reduction,'parent',ax2d);
        elseif nargin==2
            p_handle=varargin{1};
            data0=varargin{2};
            p_handle=imagesc_highspeed(p_handle,log10(data0),settings.resolution_reduction,'parent',ax2d);
        end
        nx=size(data0,1); nz=size(data0,2);
        x_limit=[1 nx];
        z_limit=[1 nz];
        axis([z_limit(1)-200 z_limit(2)+200 x_limit(1)-200 x_limit(2)+200])
    end
%%
% plot the data
p_2d=plot_2d_data(data);
clim_values=get(ax2d,'Clim');
set(h_clim_min,'string',num2str(clim_values(1)));
set(h_clim_max,'string',num2str(clim_values(2)));
set(ax2d,'YDir','normal');
hold on;
xlabel('Z-direction (pixels)');
ylabel('X-direction (pixels)');

%
j=0; %image plate counter
h=[]; %handles for polygon corners that define the image plate
p_trapez=[]; %the temporary plot of the tetragon
xy1=[];
%Instruction to the user
questdlg(['Image plates must be selected in order of low angle' ...
    'to high angle'],'Select IP','OK','OK');
%Informs if a user-interactable object is drawn
is_ui_trapez_drawn=false;

%wait for user interaction to complete. Ends with push on "Finish"-button
uiwait(fig2d)

    function select_ip(~,~)
        %First the user draws a rectangle with mouse-drag. On release a
        %user-interactable polygon is created. The user can adjust the size
        %an position
        j=1+j;
        try
            while 1
                rect=getrect(fig2d);
                % avoid that rect is a single point so all dragable points
                % are plotted on top of each other
                if rect(3)>100 && rect(4)>100
                    break
                end
            end
            z0=rect(1); x0=rect(2);
            zw=rect(3); xw=rect(4);
            z_initial=[z0 z0+zw z0+zw z0]';
            x_initial=[x0 x0 x0+xw x0+xw]';
            
            %ensure that the rectangle is within the data-area
            x_initial(x_initial<x_limit(1))=x_limit(1); 
            x_initial(x_initial>x_limit(2))=x_limit(2);
            z_initial(z_initial<z_limit(1))=z_limit(1); 
            z_initial(z_initial>z_limit(2))=z_limit(2);
            
            % Plot a rectangle which the user can edit by dragging
%             [h, p_tetragon] =ui_tetragon([z_initial, x_initial]);
            [h, p_trapez] =ui_trapez([z_initial, x_initial]);
            is_ui_trapez_drawn=true; %used by "finish" to check
            %if an unfinished polygon is drawn
        catch err
            % If "Next IP" was pressed, without drawing a rectangle before
            % "finish" was pressed then a error is thrown because getrect
            % is interupted by closing the figure.
            if strcmp(err.identifier,'MATLAB:class:InvalidHandle')
                % catch the error if getrect is interrupted by closing
                % the figure
                %disp('Error caught in select_ip')
            else
                % rethrow(err)
            end
        end
    end

    function undo_ip(~,~)
        if is_ui_trapez_drawn
            if j>0
                delete(p_trapez);
                delete(h(1));delete(h(2)); delete(h(3)); delete(h(4));
                delete(h(5));delete(h(6)); delete(h(7)); delete(h(8));
                is_ui_trapez_drawn=false;
                j=j-1;
            end
        end
    end

    function select_ip_ui(obj,event)
        if is_ui_trapez_drawn
            finalize_polygon;
        end
        select_ip(obj,event)
    end

    function finish(~,~)
        if j>0
            if is_ui_trapez_drawn
                finalize_polygon;
            end
            uiresume(gcbf)
            close(fig2d)
        else
            fprintf('You have not selected any image plates yet\n')
        end
    end

    function update_length(obj,~)
        length_new=str2double(get(obj,'String'));
        length_old=sqrt(sum((xy1(2,:)-xy1(1,:)).^2));
        xy_move=(xy1(2,:)-xy1(1,:))./length_old*(length_new-length_old);
        h(2).setPosition(xy1(2,1)+xy_move(1),xy1(2,2)+xy_move(2));
        h(3).setPosition(xy1(3,1)+xy_move(1),xy1(3,2)+xy_move(2));
    end

    function update_width(obj,~)
        width_new=str2double(get(obj,'String'));
        
        width_old=distance_between_two_parallel_line(xy1(1,:),xy1(4,:),...
                                                    xy1(1,:)-xy1(2,:));
        move_dir=xy1(1,:)-xy1(4,:);
        move=(width_new-width_old)*move_dir/norm(move_dir);
        
        xy1(1,:)=xy1(1,:)+move;
        xy1(2,:)=xy1(2,:)+move;
        
        h(1).setPosition(xy1(1,1),xy1(1,2));
        h(2).setPosition(xy1(2,1),xy1(2,2));
    end

    function d=distance_between_two_parallel_line(x1,x2,r1)
        x1=x1(:); x2=x2(:); r1=r1(:);
        
        r=x2-x1;
        
        r_proj=dot(r,r1)/norm(r1);
        
        d=sqrt(norm(r).^2-r_proj.^2);
        
    end

    function finalize_polygon()
        %get positions of corners
        final_position=xy1(1:4,:);
        %plot closed polygon
        xy_closed=[final_position; final_position(1,:)];
        plot(xy_closed(:,1),xy_closed(:,2),'b','linewidth',1);
        
        %clean up and save
        delete(h(1));delete(h(2)); delete(h(3)); delete(h(4));
        delete(h(5));delete(h(6)); delete(h(7)); delete(h(8));
        ip_position{j}=round([final_position(:,2) final_position(:,1)]);
        
        %Determine length of image plate, as used to determine the off-set
        %in Z-direction
        
        [ip_origin, ~]=get_ip_origin_and_length(ip_position{j});
        %mark the left and rigth corner of the image plate
        plot(ip_origin(2),ip_origin(1),'ok');
        
        is_ui_trapez_drawn=false;
        
        
    end

    function [h,poly_plot]=ui_trapez(xy0)
        
        xy1(1:4,:)=xy0;
        xy1(5,:)=(xy1(2,:)+xy1(1,:))/2;
        xy1(6,:)=(xy1(3,:)+xy1(2,:))/2;
        xy1(7,:)=(xy1(4,:)+xy1(3,:))/2;
        xy1(8,:)=(xy1(1,:)+xy1(4,:))/2;
        %Create as user-dragable polygon from a list of coordinates.
        hold on % do not overwrite
        r1=[]; r2=[];
        update_principal_directions;
        
        poly_plot=plot([xy0(:,1); xy0(1,1)],[xy0(:,2); xy0(1,2)],'r','linewidth',1);
        
        h(1) = impoint(gca,xy1(1,:));
        h(2) = impoint(gca,xy1(2,:));
        h(3) = impoint(gca,xy1(3,:));
        h(4) = impoint(gca,xy1(4,:));
        h(5) = impoint(gca,xy1(5,:));
        h(6) = impoint(gca,xy1(6,:));
        h(7) = impoint(gca,xy1(7,:));
        h(8) = impoint(gca,xy1(8,:));
        
        setColor(h(1),'blue');
        setColor(h(2),'red');
        setColor(h(4),'black');
        setColor(h(3),'red');

        
        %Add call back to every point
        %Corners
        addNewPositionCallback(h(1),@(p) move_h1); 
        addNewPositionCallback(h(2),@(p) move_h2); 
        addNewPositionCallback(h(3),@(p) move_h3); 
        addNewPositionCallback(h(4),@(p) move_h4); 
        %lines - midpoints
        addNewPositionCallback(h(5),@(p) move_h5); 
        addNewPositionCallback(h(7),@(p) move_h7); 
        addNewPositionCallback(h(6),@(p) move_h6);
        addNewPositionCallback(h(8),@(p) move_h8);

        
        %Now follows the individual callbacks for all points. Each callbac
        %defines how the trapez is moved. 
        function move_h1
            xy_old=xy1;
            xy_new=h(1).getPosition;
            d=xy_new-xy_old(1,:);
            
            d_proj=r1*dot(r1,d);
            xy1(1,:)=xy_old(1,:)+d_proj;
            
            update_mid_position;
            h(1).setPosition(xy1(1,:))
            h(5).setPosition(xy1(5,:))
            h(8).setPosition(xy1(8,:))   
            draw_tetragon;
        end
        function move_h2
            modifiers = get(gcf,'currentModifier');       %(Use an actual figure number if known)
            altIsPressed = ismember('alt',modifiers);
            if ~altIsPressed
                xy_old=xy1;
                xy_new=h(2).getPosition;
                d=xy_new-xy_old(2,:);
                
                d_proj=r1*dot(r1,d);
                xy1(2,:)=xy_old(2,:)+d_proj;
                
                update_mid_position;
                h(2).setPosition(xy1(2,:))
                h(1).setPosition(xy1(1,:))
                h(3).setPosition(xy1(3,:))
            else
                xy_old=xy1;
                xy_base=xy_old(4,:);
                old_dir=xy_old(2,:)-xy_base;
                old_dir=old_dir/norm(old_dir);
                xy_new=h(2).getPosition;
                new_dir=xy_new-xy_base;
                new_dir=new_dir/norm(new_dir);
                
                %find angle between vectors using cross product
                angle=get_angle(new_dir,old_dir);
                
                A0=rotation_matrix(angle);
                
                %rotate all points except the pivot point
                xy1(1,:)=rotate_point(A0,xy1(4,:),xy1(1,:));
                xy1(2,:)=rotate_point(A0,xy1(4,:),xy1(2,:));
                xy1(3,:)=rotate_point(A0,xy1(4,:),xy1(3,:));
                update_mid_position;
                update_principal_directions
                h(1).setPosition(xy1(1,:))
                h(2).setPosition(xy1(2,:))
                h(3).setPosition(xy1(3,:))
                
            end
            draw_tetragon;
        end
        function move_h3
            modifiers = get(gcf,'currentModifier');       %(Use an actual figure number if known)
            altIsPressed = ismember('alt',modifiers);
            if ~altIsPressed
                xy_old=xy1;
                xy_new=h(3).getPosition;
                d=xy_new-xy_old(3,:);
                
                d_proj=r1*dot(r1,d);
                xy1(3,:)=xy_old(3,:)+d_proj;
                
                update_mid_position;
                h(3).setPosition(xy1(3,:))
                h(2).setPosition(xy1(2,:))
                h(4).setPosition(xy1(4,:))
            else
                xy_old=xy1;
                xy_base=xy_old(4,:);
                old_dir=xy_old(3,:)-xy_base;
                old_dir=old_dir/norm(old_dir);
                xy_new=h(3).getPosition;
                new_dir=xy_new-xy_base;
                new_dir=new_dir/norm(new_dir);
                
                %find angle between vectors using cross product
                angle=asind(new_dir(1)*old_dir(2)-new_dir(2)*old_dir(1));
                
                A0=rotation_matrix(angle);
                
                %rotate all points except the pivot point
                xy1(3,:)=rotate_point(A0,xy1(4,:),xy1(3,:));
                xy1(1,:)=rotate_point(A0,xy1(4,:),xy1(1,:));
                xy1(2,:)=rotate_point(A0,xy1(4,:),xy1(2,:));
                update_mid_position;
                update_principal_directions
                h(3).setPosition(xy1(3,:))
                h(1).setPosition(xy1(1,:))
                h(2).setPosition(xy1(2,:))
                
            end
            draw_tetragon;
        end       
        function move_h4
            modifiers = get(gcf,'currentModifier');
            altIsPressed = ismember('alt',modifiers);
            if ~altIsPressed
                xy_old=xy1;
                xy_new=h(4).getPosition;
                d=xy_new-xy_old(4,:);
                
                d_proj=r1*dot(r1,d);
                xy1(4,:)=xy_old(4,:)+d_proj;
                
                update_mid_position;
                h(4).setPosition(xy1(4,:))
                h(1).setPosition(xy1(1,:))
                h(3).setPosition(xy1(3,:))
            else %move the entire trapez
                xy_old=xy1;
                xy_new=h(4).getPosition;
                d=xy_new-xy_old(4,:);
               
                xy1=xy_old+repmat(d,8,1);
                
                update_mid_position;
                
                h(1).setPosition(xy1(1,:))
                h(2).setPosition(xy1(2,:))
                h(3).setPosition(xy1(3,:))
                h(4).setPosition(xy1(4,:))
                h(5).setPosition(xy1(5,:))
                h(6).setPosition(xy1(6,:))
                h(7).setPosition(xy1(7,:))
                h(8).setPosition(xy1(8,:))
            end  
            draw_tetragon;
        end
        function move_h5
            xy_old=xy1;
            xy_new=h(5).getPosition;
            d=xy_new-xy_old(5,:);
            
            d_proj=r2*dot(r2,d);
            xy1(1,:)=xy_old(1,:)+d_proj;
            xy1(2,:)=xy_old(2,:)+d_proj;
            
            update_mid_position;
            
            h(1).setPosition(xy1(1,:))
            h(2).setPosition(xy1(2,:))
            h(6).setPosition(xy1(6,:))
            h(8).setPosition(xy1(8,:))
            
            draw_tetragon;
        end
        function move_h6
            xy_old=xy1;
            xy_new=h(6).getPosition;
            d=xy_new-xy_old(6,:);
            
            d_proj=r1*dot(r1,d);
            
            xy1(2,:)=xy_old(2,:)+d_proj;
            xy1(3,:)=xy_old(3,:)+d_proj;
            update_mid_position;
            h(2).setPosition(xy1(2,:))
            h(3).setPosition(xy1(3,:))
            h(5).setPosition(xy1(5,:))
            h(7).setPosition(xy1(7,:))
            
            draw_tetragon;
        end
        function move_h7
            xy_old=xy1;
            xy_new=h(7).getPosition;
            d=xy_new-xy_old(7,:);
            
            d_proj=r2*dot(r2,d);

            xy1(3,:)=xy_old(3,:)+d_proj;
            xy1(4,:)=xy_old(4,:)+d_proj;
            
            update_mid_position;
            h(3).setPosition(xy1(3,:))
            h(4).setPosition(xy1(4,:))
            h(6).setPosition(xy1(6,:))
            h(8).setPosition(xy1(8,:))
            
            draw_tetragon;
        end
        function move_h8
            xy_old=xy1;
            xy_new=h(8).getPosition;
            d=xy_new-xy_old(8,:);
            
            d_proj=r1*dot(r1,d);
            
            xy1(1,:)=xy_old(1,:)+d_proj;
            xy1(4,:)=xy_old(4,:)+d_proj;
            update_mid_position;
            h(1).setPosition(xy1(1,:))
            h(4).setPosition(xy1(4,:))
            h(5).setPosition(xy1(5,:))
            h(7).setPosition(xy1(7,:))
            
            draw_tetragon;  
        end
        
        function update_mid_position()
            xy1(5,:)=(xy1(2,:)+xy1(1,:))/2;
            xy1(6,:)=(xy1(3,:)+xy1(2,:))/2;
            xy1(7,:)=(xy1(4,:)+xy1(3,:))/2;
            xy1(8,:)=(xy1(1,:)+xy1(4,:))/2;
        end
        
        function update_principal_directions
            r1=xy1(4,:)-xy1(3,:);
            r1=r1/norm(r1);
            r2=[-r1(2) r1(1)];
        end
        
        function draw_tetragon()
            %Get the current position from the impoint objects in "h"
            % and updates the coordinates of "pp" which is spanned by
            % those points
            
            xy_shape=[xy1(1:4,:); xy1(1,:)];
            set(poly_plot,'XData',xy_shape(:,1),'YData',xy_shape(:,2))

            ip_length=sqrt(sum((xy1(2,:)-xy1(1,:)).^2));
            set(h13,'string',num2str(round(ip_length)));
            ip_width=distance_between_two_parallel_line(xy1(1,:),xy1(4,:),...
                                                    xy1(1,:)-xy1(2,:));
            set(h15,'string',num2str(round(ip_width)));
        end
        
        function xy1_rot=rotate_point(A,xy0,xy1)
            %xy0 and xy1 are 1 x 2 matrices
            %Rotates the point xy1 around the point xy0. The rotation angle
            %is defined by the rotation matrix A
            r=xy1-xy0;
            r_rot=r*A;
            xy1_rot=xy0+r_rot;
        end
        function A=rotation_matrix(phi)
            A=[cosd(phi) -sind(phi);sind(phi) cosd(phi)];
        end
%         function xi=cyclic(xi,n)
%             xi=mod(xi-1,n)+1;
%         end
    end

end

function [data_cut, ip_positions_cut] =cut_raw_data(data,ip_positions)
%given a tetragon defining the image plate, this function cuts out the
% smallest surrounding rectangle from the entire scanned image.
data_cut=cell(length(ip_positions),1);
ip_positions_cut=cell(length(ip_positions),1);
nx=size(data,1);
nz=size(data,2);
x_limit=[1 nx];
z_limit=[1 nz];

for j=1:length(ip_positions)
    borders=ip_positions{j};
    xcut=round([min(borders(:,1)) max(borders(:,1))]);
    zcut=round([min(borders(:,2)) max(borders(:,2))]);
    xcut(xcut<x_limit(1))=x_limit(1); xcut(xcut>x_limit(2))=x_limit(2);
    zcut(zcut<z_limit(1))=z_limit(1); zcut(zcut>z_limit(2))=z_limit(2);
    data_cut{j}=data(xcut(1):xcut(2),zcut(1):zcut(2));
    ip_positions_cut{j}=[borders(:,1)-xcut(1), borders(:,2)-zcut(1)];
end
end

function [data_new, borders_new]=transform_to_rectangular_grid(varargin)
global verbose
if nargin ==2
    data=varargin{1};
    borders=varargin{2};
    unwarp_shape='twolines';
elseif nargin==3
    data=varargin{1};
    borders=varargin{2};
    unwarp_shape=varargin{3};
else
    error('Unknown number of arguments for ''transform_to_rectangular_grid''')
end
%Get the corners of the region defining the IP
[xi,zi] = get_bottom_left_corner(borders(:,1),borders(:,2));
xz0=[xi,zi]';
[xi,zi] = get_top_left_corner(borders(:,1),borders(:,2));
xz1=[xi,zi]';
[xi,zi] = get_bottom_right_corner(borders(:,1),borders(:,2));
xz2=[xi,zi]';
[xi,zi] = get_top_right_corner(borders(:,1),borders(:,2));
xz3=[xi,zi]';

% The IP is spanned by r1 and r2 starting from xy0

% xz1 ----------------- xz3
% |          r4 ->       |
% |                      |
% | r1 ^|          |^ r3 |
% |                      |
% |          r2 ->       |
% xz0-------------------xz2

r1=xz1-xz0;
r2=xz2-xz0;
r3=xz3-xz2;
r4=xz3-xz1;

if verbose > 2
    disp('The angles between the sides of the image plate are:')
    disp(get_angle(r1,r2))
    disp(get_angle(r2,r3))
    disp(get_angle(r3,r4))
    disp(get_angle(r4,r1))
    disp('degrees.')
end



switch unwarp_shape
    case 'parallelogram'
        xm=round(norm(r1));
        zm=round(norm(r2));
        x_new=(0:xm)';
        z_new=(0:zm)';
        r1n=r1/xm;
        r2n=r2/zm;
        
        [Z_new,X_new]=meshgrid(z_new,x_new);
        
        X_new_in_old=X_new*r1n(1)+Z_new*r2n(1)+xz0(1);
        Z_new_in_old=X_new*r1n(2)+Z_new*r2n(2)+xz0(2);
        
    case 'tetragon'
        xm=round(norm(r1));
        zm=round(norm(r2));
        x_new=(0:xm)';
        z_new=(0:zm)';
        r1n=r1/xm;
        r2n=r2/zm;
        r3n=r3/xm;
        
        [Z_new,X_new]=meshgrid(z_new,x_new);
        X_new_in_old=X_new.*(r1n(1).*(1-Z_new/zm)+r3n(1).*Z_new/zm)+...
            Z_new.* r2n(1)+...
            xz0(1);
        Z_new_in_old=X_new.*(r1n(2).*(1-Z_new./zm)+r3n(2).*Z_new./zm)+...
            Z_new.* r2n(2)+...
            xz0(2);
    case {'twoline', 'twolines'}
        r2o=[-r2(2), r2(1),];
        r4o=[-r4(2), r4(1),];
        r_vert=(r2o+r4o)/2;
        r_vert=r_vert/norm(r_vert);
       
        %intersection between vertical lines and top horizontal lines
        [~, xz0proj]=vector_lines_intersect(xz0,r2,xz1,r_vert);
        [~, xz1proj]=vector_lines_intersect(xz1,r4,xz0,r_vert);
        [~, xz2proj]=vector_lines_intersect(xz2,r2,xz3,r_vert);
        [~, xz3proj]=vector_lines_intersect(xz3,r4,xz2,r_vert);
        
        
        
        xz0_org2proj=xz0proj-xz0;
        xz1_org2proj=xz1proj-xz1;
        xz2_org2proj=xz2proj-xz2;
        xz3_org2proj=xz3proj-xz3;
        
        if verbose > 2
            fig_test=figure; hold on;
            plot([xz0(2) xz1(2) xz3(2) xz2(2) xz0(2) ],...
                [xz0(1) xz1(1)  xz3(1) xz2(1) xz0(1)],'-r')
            plot([xz0proj(2) xz1proj(2)  xz3proj(2) xz2proj(2) xz0proj(2)],...
                [xz0proj(1) xz1proj(1)  xz3proj(1) xz2proj(1) xz0proj(1)],'-b')
            
            plot([xz0(2) xz1(2) xz2(2) xz3(2)],...
                [xz0(1) xz1(1) xz2(1) xz3(1)],'or')
            plot([xz0proj(2) xz1proj(2) xz2proj(2) xz3proj(2)],...
                [xz0proj(1) xz1proj(1) xz2proj(1) xz3proj(1)],'xb')
            
            legend('Original','Projected')
        end
        if dot(xz0_org2proj,r2)<0;
            xz0=xz0proj;
        elseif dot(xz1_org2proj,r4)<0;
            xz1=xz1proj;
        else
            error('')
        end
        
        
        if dot(xz2_org2proj,r2)>0;
            xz2=xz2proj;
        elseif dot(xz3_org2proj,r4)<0;
            xz3=xz3proj;
        end
        
        if verbose >2
            figure(fig_test)
            plot([xz0(2) xz1(2) xz3(2) xz2(2) xz0(2) ],...
                [xz0(1) xz1(1)  xz3(1) xz2(1) xz0(1)],'-k')
        end
        
        r1=xz1-xz0;
        r2=xz2-xz0;
        r3=xz3-xz2;
        xm=round(norm(r1));
        zm=round(norm(r2));
        x_new=(0:xm)';
        z_new=(0:zm)';
        
        r1n=r1/xm;
        r2n=r2/zm;
        r3n=r3/xm;
        
        [Z_new,X_new]=meshgrid(z_new,x_new);

        x_new2old= @(x,z) x.*(r1n(1).*(1-z/zm)+r3n(1).*z/zm)+...
            z.* r2n(1)+...
            xz0(1);
        z_new2old= @(x,z) x.*(r1n(2).*(1-z./zm)+r3n(2).*z./zm)+...
            z.* r2n(2)+...
            xz0(2);
        X_new_in_old=x_new2old(X_new,Z_new);
        Z_new_in_old=z_new2old(X_new,Z_new);
        
        
    case {'oneline', 'onelines'}
        %only the direction of top horizontal line is used to unwarp the 
        %IP. Same as above if the two "horizontal" lines are parallel
        r_vert=[-r4(2), r4(1),];
        
        %intersection between vertical lines and top horizontal lines
        [~, xz0proj]=vector_lines_intersect(xz0,r2,xz1,r_vert);
        [~, xz1proj]=vector_lines_intersect(xz1,r4,xz0,r_vert);
        [~, xz2proj]=vector_lines_intersect(xz2,r2,xz3,r_vert);
        [~, xz3proj]=vector_lines_intersect(xz3,r4,xz2,r_vert);
        
        %vector from original point to projected point
        xz0_org2proj=xz0proj-xz0;
        xz1_org2proj=xz1proj-xz1;
        xz2_org2proj=xz2proj-xz2;
        xz3_org2proj=xz3proj-xz3;
        
        if verbose > 2
            fig_test=figure; hold on;
            plot([xz0(2) xz1(2) xz3(2) xz2(2) xz0(2) ],...
                [xz0(1) xz1(1)  xz3(1) xz2(1) xz0(1)],'-r')
            plot([xz0proj(2) xz1proj(2)  xz3proj(2) xz2proj(2) xz0proj(2)],...
                [xz0proj(1) xz1proj(1)  xz3proj(1) xz2proj(1) xz0proj(1)],'-b')
            
            plot([xz0(2) xz1(2) xz2(2) xz3(2)],...
                [xz0(1) xz1(1) xz2(1) xz3(1)],'or')
            plot([xz0proj(2) xz1proj(2) xz2proj(2) xz3proj(2)],...
                [xz0proj(1) xz1proj(1) xz2proj(1) xz3proj(1)],'xb')
            
            legend('Original','Projected')
        end
        
        %ensure that both projected and original points are included
        if dot(xz0_org2proj,r2)<0;
            xz0=xz0proj;
        elseif dot(xz1_org2proj,r4)<0;
            xz1=xz1proj;
        else
            error('')
        end
        if dot(xz2_org2proj,r2)>0;
            xz2=xz2proj;
        elseif dot(xz3_org2proj,r4)<0;
            xz3=xz3proj;
        end
        
        if verbose >2
            figure(fig_test)
            plot([xz0(2) xz1(2) xz3(2) xz2(2) xz0(2) ],...
                [xz0(1) xz1(1)  xz3(1) xz2(1) xz0(1)],'-k')
        end
        
        r1=xz1-xz0;
        r2=xz2-xz0;
        xm=round(norm(r1));
        zm=round(norm(r2));
        x_new=(0:xm)';
        z_new=(0:zm)';
        
        r1n=r1/xm;
        r2n=r2/zm;
        
        [Z_new,X_new]=meshgrid(z_new,x_new);
        
        x_new2old= @(x,z) x.*r1n(1)+...
            z.* r2n(1)+...
            xz0(1);
        z_new2old= @(x,z) x.*r1n(2)+...
            z.* r2n(2)+...
            xz0(2);
        
        X_new_in_old=x_new2old(X_new,Z_new);
        Z_new_in_old=z_new2old(X_new,Z_new);
        
        %         Transform the ip-borders to the new coordinates
        B=[r1n r2n];
        r_new=@(r_old) B\(r_old-xz0);
        borders_new=zeros(size(borders));
        for i=1:4
            borders_new(i,:)=(r_new(borders(i,:)'))';
        end
        if verbose >2
            figure(fig_test)
            plot(borders_new(:,2),borders_new(:,1),'-k')
        end
    otherwise
        error('Unknown ''unwarp_shape'' requested')
end

if ~exist('borders_new','var')
disp(['The IP-border coordinates has not been transformed to the new'...
    ' coordinate system. Maybe it has not been implemented for the '...
    'selected unwarp method. The IP-borders are only used for the '...
    'optional initial mask']);
borders_new=borders;
end

% interp2 assumes that the default grid points cover the rectangular
%region, X=1:n and Y=1:m, where [m,n] = size(V)
data_new=interp2(data,Z_new_in_old,X_new_in_old,'linear',0);


end

%% Define mask
function [mask_poly, mask_inout, mask_2d]=define_mask(data,borders,settings)
%User interface to mask the image plates before integration.
global dfs

current_number=0; %used for undo/redo feature

%The masks are defined by these variables
mask_poly={};  %dim [n_masks, 1]
mask_inout=[]; %dim [n_masks, 1]
mask=[];       %dim [n_masks, 1]
mask_plot=[];
nx=size(data,1); nz=size(data,2);
temp_obj=[];
if isempty(borders)
borders=[1 1; 1 nz;nx nz;nx 1];
end

fig_mask = figure('PaperUnits',get(0,'defaultfigurePaperUnits'),...
                  'Name','Mask image plate',...
                  'NumberTitle','off',...
                  'Position',[50 50 dfs], ...
                  'Visible','on',...
                  'menubar','none',...
                  'toolbar','figure');

ax_mask = axes('Parent',fig_mask,...
               'units','normalized',...
               'position',[0.1 0.18 0.85 0.77 ]);

% Create the button group,
h1 = uibuttongroup('visible','on', ...
                   'units','pixels', ...
                   'Position',[5 20 215 40]);
% which contains 4 buttons:
% Finish: Press this button when all masks have been drawn.
% Undo: Undo the previous mask
% Redo:
% Hide: Hide or unhide all masks

h11 = uicontrol(...
    'Parent',h1,...
    'Callback',@finish,...
    'Position',[5 7 50 22],...
    'String','Finish' );%#ok
h12 = uicontrol(...
    'Parent',h1,...
    'Callback',@undo_callback,...
    'Position',[60 7 40 22],...
    'String','Undo');%#ok
h13 = uicontrol(...
    'Parent',h1,...
    'Callback',@redo_callback,...
    'Position',[105 7 40 22],...
    'String','Redo');%#ok

h14 = uicontrol(...
    'Parent',h1,...
    'Callback',@hide_callback,...
    'Position',[150 7 40 22],...
    'String','Hide');


% Create a group with 2 radio buttons: Keep and remove. Keep
% means that everything outside the selected area will be masked. Remove
% the selected area is masked. The mask will be colored accordingly.
h2 = uibuttongroup('visible','on','units','pixels','Position',[230 20 220 40]);
h21 = uicontrol('Style','Radio','String','Remove','Units','pixels', ...
    'pos',[10 7 90 23],'parent',h2,'HandleVisibility','on');%#ok
h22 = uicontrol('Style','Radio','String','Keep','Units','pixels', ...
    'pos',[120 7 90 23],'parent',h2,'HandleVisibility','on');
% Initialize some button group properties.
set(h2,'SelectedObject',h22);  % No selection


% Create a group with 3 buttons to select the shape of the mask. When one
% shape is selected the buttons are inactive/hidden. To select another
% shape, first press "escape". This is also written in the userinterface.
% While one shape is selected it is still possible to switch between
% Keep and remove mask using the radio-buttons.
h3 = uibuttongroup('visible','on','units','pixels',...
    'Position',[460 20 240 40]);
h31 = uicontrol('style','pushbutton','Parent',h3,'Position',[5 7 70 22],...
    'String','Rectangle','callback', @draw_mask );%#ok
h32 = uicontrol('style','pushbutton','Parent',h3,'Position',[80 7 70 22],...
    'String','Circle' ,'callback', @draw_mask );%#ok
h33 = uicontrol('style','pushbutton','Parent',h3,'Position',[155 7 70 22],...
    'String','Polygon','callback', @draw_mask  );%#ok

%This object is shown instead of shape selection buttons, when a specific
%shape is selected
h4 = uicontrol('style','text',...
    'visible','off',...
    'string',{' ' ,'Press Escape to select another mask shape'},... % cell is a work-around to make a newline.
    'units','pixels',...
    'Position',[460 20 240 40]);


[h_clim_min, h_clim_max]=add_colorlimit_fields(fig_mask,ax_mask,[720 20 120 45]);

h31 = uicontrol('style','pushbutton','Parent',fig_mask,'Position',[860 27 40 20],...
    'String','Save IP','callback', @save_data );%#ok
h32 = uicontrol('style','pushbutton','Parent',fig_mask,'Position',[910 27 40 20],...
    'String','Save mask','callback', @save_mask );%#ok

ef_height=0.05;
ef_width=0.05;
h_edge_left = uicontrol('style','edit',...
    'units','normalized',...
    'position',[0.1-ef_width 0.18 ef_width ef_height],...
    'parent',fig_mask,'String','0',...
    'callback',(@(src,evt) update_edge_mask(borders)));
h_edge_right = uicontrol('style','edit',...
    'units','normalized',...
    'position',[0.95 0.9 ef_width ef_height],...
    'parent',fig_mask,'String','0',...
    'callback',(@(src,evt) update_edge_mask(borders)));
h_edge_top = uicontrol('style','edit',...
    'units','normalized',...
    'position',[0.9 0.95 ef_width ef_height],...
    'parent',fig_mask,'String','0',...
    'callback',(@(src,evt) update_edge_mask(borders)));
h_edge_bottom = uicontrol('style','edit',...
    'units','normalized',...
    'position',[0.1 0.18-ef_height ef_width ef_height],...
    'parent',fig_mask,'String','0',...
    'callback',(@(src,evt) update_edge_mask(borders)));

%Deactivate shape-selection buttons and show h4 when a call-back is running


callbackrunning = false;
set(h3,'userdata',callbackrunning);
continuedrawingmask=false;
set(ax_mask,'ButtonDownFcn',@draw_mask)




%plot the data
imagesc_highspeed(log10(data),...
    settings.resolution_reduction,...
    'Parent',ax_mask);
set(ax_mask,'YDir','normal');
        xlabel('Z-direction (pixels)')
        ylabel('X-direction (pixels)')
hold on
clim_values=get(ax_mask,'Clim');
set(h_clim_min,'string',num2str(clim_values(1)));
set(h_clim_max,'string',num2str(clim_values(2)));


update_edge_mask(borders)

%The program waits here so the user can interact with the figures, until
%the "Finish"-button is pressed.
waitfor(fig_mask,'Userdata','finished')

%User interaction must have finished. Cleaning up and preparing output.
%Delete undon masks
if current_number > 0
    mask_poly=mask_poly(1:current_number);
    mask=mask(:,:,1:current_number);
    mask_inout=mask_inout(1:current_number);
else % if no masks have been made, make a dummy unity mask.
    mask=ones(size(data));
    mask_poly=[];
    mask_inout=1;
end


% combine all the individual masks into a single mask.
mask_2d=single(prod(mask,3));
close(fig_mask);

%Nested functions
    function finish(~,~)
        %presses the escape button to end the imrect,impoly, and imellipse
        send_escape;
        %Mark that user-interaction has finished. This breaks the waitfor
        set(fig_mask,'Userdata','finished')
    end

    function save_data(~,~)
        [filename,path] = uiputfile('*.mat','Save file name','');
        if isequal(filename,0)
            output_file_name=0;
        else
            output_file_name=[path filename];
        end
        if output_file_name~=0
            save(output_file_name,'data')
        end
    end

    function save_mask(~,~)
        [filename,path] = uiputfile('*.mat','Save file name','');
        if isequal(filename,0)
            output_file_name=0;
        else
            output_file_name=[path filename];
        end
        if output_file_name~=0
            if current_number > 0
            mask2d=single(prod(mask(:,:,current_number),3)); %#ok
            else % if no masks have been made, make a dummy unity mask.
            mask2d=ones(size(data)); %#ok
            end
            save(output_file_name,'mask2d')
        end
    end

    function undo_callback(~,~)
        %undo a mask. Every mask has a number. Masks lower of equal to
        %current_number are shown
        if current_number>0
            set(mask_plot(current_number),'Visible','off')
            current_number=current_number-1;
        else
            disp('Cannot undo mask')
            beep;
        end
    end

    function redo_callback(~,~)
        %redo a mask which was undone, by undo_callback
        if current_number<length(mask_poly)
            current_number=current_number+1;
            set(mask_plot(current_number),'Visible','on')
        else
            disp('Cannot redo mask')
            beep;
        end
    end

    function hide_callback(~,~)
        %Hide or reshow all masks
        if current_number>0
            switch get(mask_plot(1),'Visible')
                case 'on'
                    set(mask_plot(1:current_number),'Visible','off')
                    set(h14,'string','Show')
                case 'off'
                    set(mask_plot(1:current_number),'Visible','on')
                    set(h14,'string','Hide')
            end
        else
            disp('No mask to hide')
            beep;
        end
    end

    function hide_buttons(~,~)
        %deactivate shape-selection buttons
        set(h3,'visible','off')
        set(h4,'visible','on') %Shows info to user
    end

    function show_buttons(~,~)
        % Reactivate shape-selection buttons
        set(h3,'visible','on')
        set(h4,'visible','off')
    end

    function update_edge_mask(edge0)
        %Get the corners of the region defining the IP
        [xi,zi] = get_bottom_left_corner(edge0(:,1),edge0(:,2));
        xz0 =[xi,zi]';
        [xi,zi] = get_top_left_corner(edge0(:,1),edge0(:,2));
        xz1=[xi,zi]';
        [xi,zi] = get_bottom_right_corner(edge0(:,1),edge0(:,2));
        xz2=[xi,zi]';
        [xi,zi] = get_top_right_corner(edge0(:,1),edge0(:,2));
        xz3=[xi,zi]';
        r1=xz1-xz0; r1=r1/norm(r1); r1c=[-r1(2); r1(1)];
        r2=xz2-xz0; r2=r2/norm(r2); r2c=[r2(2); -r2(1)];
        r3=xz3-xz2; r3=r3/norm(r3); r3c=[r3(2); -r3(1)];
        r4=xz3-xz1; r4=r4/norm(r4); r4c=[-r4(2); r4(1)];
        
        d1c=str2double(get(h_edge_left,'string'));
        d2c=str2double(get(h_edge_bottom,'string'));
        d3c=str2double(get(h_edge_right,'string'));
        d4c=str2double(get(h_edge_top,'string'));
        
        [~, xz0i]=vector_lines_intersect(xz0+d1c*r1c,r1,xz0+d2c*r2c,r2);
        [~, xz1i]=vector_lines_intersect(xz0+d1c*r1c,r1,xz1+d4c*r4c,r4);
        [~, xz2i]=vector_lines_intersect(xz2+d3c*r3c,r3,xz0+d2c*r2c,r2);
        [~, xz3i]=vector_lines_intersect(xz1+d4c*r4c,r4,xz2+d3c*r3c,r3);
        
        edge_new=[xz0i xz1i xz3i xz2i ]';
        
        mask0=edge_new([1:end 1],:);
        draw_edge_mask(mask0)
        
        % The IP is spanned by r1 and r2 starting from xy0
        
        %                 r4 ->
        %     xz1 ----------------- xz3
        %      |         r4c v        |
        %      |                      |
        %r1 |^ | r1c ->        <- r3c | ^| r3
        %      |                      |
        %      |         r2c|^        |
        %      xz0-------------------xz2
        %                r2 ->
        
        % r1c is orthogonal to r1 and is pointing the center of the tetragon
    end

    function draw_mask(obj,~)
        %checks if another callback is currently running
        if get(h3,'userdata')
            return
        end
        continuedrawingmask=true;
        %         determine which button was pressed last
        switch get(obj,'String')
            case 'Rectangle'
                mask_rectangle;
            case 'Circle'
                mask_circle;
            case 'Polygon'
                mask_polygon;
            otherwise
                disp('shape not selected')
        end
        
        %Draw new mask of same shape unless continuedrawingmask was set to
        %false. This happens when call to mask_* has been ended by
        %escape-press by user.
        
        if continuedrawingmask
            draw_mask(obj,[]);
        end
    end

    function mask_rectangle(~,~)
        %Make of rectangular shape
        
        %Ensure that another mask-shape is not chosen while current
        %callback is running.
        set(h3, 'userdata',true);  %callbackrunning
        hide_buttons;
        %User draws object
        temp_obj=imrect(ax_mask);
        if isempty(temp_obj) % Check is imrect returned empty because
            % escape was pressed
            %make it possible to select another mask shape again
            set(h3, 'userdata',false);  %callbackrunning
            show_buttons;
            continuedrawingmask=false;
            return;
        end
        
        %list of 4 coordinates that defines the corners of the rectangle
        poly=rect2poly(temp_obj.getPosition);
        temp_obj.delete;
        %check if the mask should be a remove or keep mask.
        mask_type=string2inout(get(get(h2,'SelectedObject'),'string'));
        draw_poly_mask(poly,mask_type);
        set(h3, 'userdata',false);  %callbackrunning
        %         uiwait(fig_mask)
    end



    function mask_circle(~,~)
        %equivalent to mask_rectangle, but draws a circle.
        
        set(h3, 'userdata',true);  %callbackrunning
        hide_buttons;
        temp_obj=imellipse(ax_mask);
        if isempty(temp_obj)
            show_buttons;
            set(h3, 'userdata',false);  %callbackrunning
            continuedrawingmask=false;
            return;
        end
        temp=ellipse2polygon(temp_obj.getPosition,64);
        poly=[temp(:,2) temp(:,1)];
        temp_obj.delete;
        mask_type=string2inout(get(get(h2,'SelectedObject'),'string'));
        draw_poly_mask(poly,mask_type);
        set(h3, 'userdata',false);  %callbackrunning
        
        %         uiwait(fig_mask)
    end

    function poly=rect2poly(rect)
        %rect contains the origo, width, and height of a rectangle. poly contains
        %the coordinates of the 4 corners instead.
        z0=rect(1); x0=rect(2);
        zw=rect(3); xw=rect(4);
        poly=[x0 x0 x0+xw x0+xw; z0 z0+zw z0+zw z0]';
    end

    function poly=ellipse2polygon(el,N)
        %Approximate an ellipse as a N-gon, i.e., polygon with N-corners)
        
        %el defines the bounding rectangle of the ellipse
        xmin=el(1);
        ymin=el(2);
        xw=el(3);
        yw=el(4);
        
        % create circle
        t = linspace(0, 2*pi, N+1)';
        x = xmin+xw/2 + xw/2 * cos(t);
        y = ymin+yw/2 + yw/2 * sin(t);
        poly=[x y];
    end

    function mask_polygon(~,~)
        %equivalent to mask_rectangle, but draws a polygon.
        set(h3, 'userdata',true);  %callbackrunning
        hide_buttons;
        temp_obj=impoly(ax_mask);
        if isempty(temp_obj)
            set(h3, 'userdata',false);  %callbackrunning
            show_buttons;
            continuedrawingmask=false;
            return;
        end
        temp=temp_obj.getPosition;
        temp_obj.delete;
        
        poly=[temp(:,2) temp(:,1)];
        mask_type=string2inout(get(get(h2,'SelectedObject'),'string'));
        draw_poly_mask(poly,mask_type);
        set(h3, 'userdata',false);  %callbackrunning

    end

    function draw_poly_mask(poly0,mask_type0)
        if current_number>0
        mn=current_number+1;
        else %number 1 is reserved for the edge mask
            mn=2; 
        end
        
        %save all information about the created mask in the common
        %variables
        poly0=round(poly0);
        mask_plot(mn)=fill_inout(poly0(:,2),poly0(:,1),mask_type0);
        mask(:,:,mn)=(poly2mask(poly0(:,2),poly0(:,1), nx, nz)==mask_type0);
        mask_poly{mn}=poly0;
        mask_inout(mn)=mask_type0;
        current_number=mn;
    end

    function draw_edge_mask(poly0)
        %Hides the previous edge mask
     try
         delete(mask_plot(1))
     catch
         %do nothing
     end
        mn=1;
        
        %save all information about the created mask in the common
        %variables
        poly0=round(poly0);
        mask_plot(mn)=fill_inout(poly0(:,2),poly0(:,1),1);
        mask(:,:,mn)=(poly2mask(poly0(:,2),poly0(:,1), nx, nz)==1);
        mask_poly{mn}=poly0;
        mask_inout(mn)=1;
        current_number=mn;
    end

    function p=fill_inout(x,y,in_or_outside)
        % Colours outside the polygon (defined by [x,y]) it 
        % in_or_outside=1, i.e. keep mask. Colours inside the polygon if 
        % in_or_outside ~=1, i.e. remove mask
        if in_or_outside
            %Color the area outside the polygon
            [~, i]=min(x);
            
            xi=0; xe=nz;
            yi=0; ye=nx;
            
            x=x(:);
            y=y(:);
            x=[x(i:end)' x(1:i-1)' x(i)];
            y=[y(i:end)' y(1:i-1)' y(i)];
            
            if xi>=min(x); xi=min(x)-1;end
            if yi>=min(y); yi=min(y)-1;end
            if xe<=max(x); xe=max(x)+1;end
            if ye<=max(y); ye=max(y)+1;end
            
            x=[xi   xi xe xe xi xi   x(1) x];
            y=[y(1) ye ye yi yi y(1) y(1) y];
            
            p=fill(x,y,[0 0 1]);
        else
            %Color the area inside the polygon
            p=fill(x,y,[1 0 0]);
        end
        set(p,'edgecolor','none');
    end
    function val=string2inout(str)
        %convert the string obtained from Keep/remove radiobuttons to a
        %numerical value
        switch str
            case 'Remove'
                val=0;
            case 'Keep'
                val=1;
        end
    end
    function send_escape(~,~)
        % Press the escape key. This is okay for ending the figure, but
        % cannot be used to abort one mask-shape and choose another. In
        % this case the user will have to manually press escape as he or
        % she is instructed to.
        % disp('Fire escape key to cancel previous imrect');
        % Create robot that fires an escape key
        robot = java.awt.Robot;
        robot.keyPress    (java.awt.event.KeyEvent.VK_DELETE);
        robot.keyRelease  (java.awt.event.KeyEvent.VK_DELETE);
        %         disp('Finished keypress')
    end
end



function mask=define_mask_from_polygons(data,mask_poly,type)
% create a maks of 0's and 1's of same size as data{j} based on masks
% given as polygons. "type" determines if mask is keep (=1) or
% remove(=0).
nx=size(data,1); nz=size(data,2);
mask=single(true(nx,nz));
for j=1:length(mask_poly)
    mask_i = poly2mask(mask_poly{j}(:,2),mask_poly{j}(:,1), nx, nz);
    %mask_i has size [nx,nz], and is 1 for all points within the polygon
    %and 0 outside
    mask=mask.*(mask_i==type(j));
end
mask=single(mask);

end

%% UI edit color limits
function [hb3, hb4]=add_colorlimit_fields(parent_handle,axis_handle,button_position)
hbg = uibuttongroup('parent',parent_handle, 'visible','on','units','pixels',...
    'Position',button_position);
hb1 = uicontrol('style','text','position',[5 25 50 20],...
    'parent',hbg,'String','Cmin'); %#ok
hb2 = uicontrol('style','text','position',[65 25 50 20],...
    'parent',hbg,'String','Cmax'); %#ok
hb3 = uicontrol('style','edit','position',[5 10 50 20],...
    'parent',hbg,'String','0',...
    'callback',@set_clim);
hb4 = uicontrol('style','edit','position',[65 10 50 20],...
    'parent',hbg,'String','0',...
    'callback',@set_clim);
    function set_clim(~,~)
        minVal = str2double(get(hb3,'string'));
        maxVal = str2double(get(hb4,'string'));
        if minVal>maxVal
            minVal=maxVal;
            set(hb3,'string',num2str(minVal));
        end
        set(axis_handle,'Clim',[minVal maxVal]);
    end
end

function [ip_origin, ip_length]=get_ip_origin_and_length(ip_position)

%get bottom left corner
[x_bl,z_bl]=get_bottom_left_corner(ip_position(:,1),...
    ip_position(:,2));
ip_origin=[x_bl-1,z_bl-1];
%get bottom right corner
[x_br,z_br]=get_bottom_right_corner(ip_position(:,1),...
    ip_position(:,2));

ip_length=sqrt((z_br-z_bl).^2+(x_br-x_bl).^2);



end

function p=imagesc_highspeed(varargin)
% For plotting very large images. It reduces the memory consumption and
% thus increases the speed of plotting and user interaction. Works by
% downsampling the plotted image
if isa(varargin{1},'matlab.graphics.primitive.Image')
    p=varargin{1};
    varargin=varargin(2:end);   
else
    p=0;
end
data=varargin{1}; % 2D-array to be plotted
step=varargin{2}; % downsampling. This is a positive integer. step = 1
% means no downsampling, i.e. plot full resolution

% define the new axes for the downsampled data
nx=size(data,1);      nz=size(data,2);
x_axis=1:step:nx-step;  z_axis=1:step:nz-step;
hs_data=zeros(length(x_axis),length(z_axis));
for ix=1:step
    for iz=1:step
        hs_data=hs_data+data(x_axis+(ix-1),z_axis+(iz-1));
    end
end
hs_data=hs_data/step^2;

if p==0
    % if more than 2 arguments are given, these are passed on to "imagesc"
    if length(varargin)>2
        p=imagesc(z_axis,x_axis,hs_data,varargin{3:end});
    else
        p=imagesc(z_axis,x_axis,hs_data);
    end
else
    % Could not figure out how to update the image correctly. Something
    % changes when dimension of the CData changes. THe axis mus be reset
    % somehow
    %So for now i am replotting
    delete(p)
    if length(varargin)>2
        p=imagesc(z_axis,x_axis,hs_data,varargin{3:end});
    else
        p=imagesc(z_axis,x_axis,hs_data);
    end
    
%     if length(varargin)>2
%         set(p,'cdata',hs_data,varargin{3:end});
%     else
%         set(p,'Cdata',hs_data);
%     end
end

end

%% Determine beam center
function [x_ui, z_ui]=get_beam_center_from_user_graphical(data)
global dfs resred
%Creates a user interface where the user can select the approximate beam
%center position with a mouse-click. The coordinates can also be input into
%the 2 editable fields.
fig_2d=figure('Outerposition',[ 0  50   dfs], ...
              'Name','Find beam center',...
              'NumberTitle','off',...
              'menubar','none',...
              'toolbar','figure');

ax_2d   =axes('Parent',fig_2d,...
              'units','normalized',...
              'position',[0.1 0.19 0.85 0.77 ]);

imagesc_highspeed(log10(data),resred,'parent',ax_2d);
set(ax_2d,'YDir','normal');
xlabel('Z-direction (pixels)')
ylabel('X-direction (pixels)')


h=uibuttongroup('visible','on',...
              'units','pixels',...
              'Position',[5 5 350 50]);

input_z=uicontrol('parent',h,...
              'Style', 'edit',...
              'String', '',...
              'Position', [120 10 80 20]);
input_x=uicontrol('parent',h,...
              'Style', 'edit',...
              'String', '',...
              'Position', [210 10 80 20]);
text_z=uicontrol('parent',h,...
              'Style', 'text',...
              'String', 'Z',...
              'Position', [120 30 80 15]); %#ok
text_x=uicontrol('parent',h,...
              'Style', 'text',...
              'String', 'X',...
              'Position', [210 30 80 15]); %#ok
select_button=uicontrol('parent',h,...
              'Style', 'pushbutton',...
              'String', 'Select with mouse',...
              'Position', [10 10 100 20],...
              'Callback',@select_with_mouse_callback); %#ok
ok_button=uicontrol('parent',h,...
              'Style', 'pushbutton',...
              'String', 'OK',...
              'Position', [300 10 40 20],...
              'Callback',@finish); %#ok

%waits here until the OK-button is pressed
uiwait(fig_2d) 

z_ui =str2double(get(input_z,'string'));
x_ui =str2double(get(input_x,'string'));

close(fig_2d);


    function select_with_mouse_callback(~,~)
        % Select a single point with the mouse
        [z,x] = ginput(1);
        % and input these coordinates into the editable fields
        set(input_z,'String',num2str(round(z)));
        set(input_x,'String',num2str(round(x)));
    end

    function finish(~,~)
       if isempty(get(input_z,'string')) || isempty(get(input_x,'string'))
           %ignore the push on OK button because a beam center has not been
           %selected
       else
           uiresume(gcbf)
       end
    end

end

%% Refine Beam center
function bc=refine_beam_center(data,bc_guess)
global verbose

xguess=double(bc_guess(1)); zguess=double(bc_guess(2));
nx=size(data,1);            nz=size(data,2);

% Cutout the data around the beamposition
indX=(-100:100)'+round(xguess);
indZ=(-101:101)'+round(zguess);
indX(indX<1 | indX>nx)=[];
indZ(indZ<1 | indZ>nz)=[];
IPBeam=data(indX,indZ);


% Fit the cross section
xcs=sum(IPBeam,2);
zcs=sum(IPBeam,1)';
% Fit a double gaussian to the beam center
ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+a2*exp(-((x-b1)/c2)^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
% Order of parameters: (amplitude1,amplitude2,position1,sigma1,sigma2,x)
opts.Lower = [0 0 0 0 0];
opts.Upper = [Inf Inf Inf Inf Inf];

opts.StartPoint = [max(xcs) max(xcs)/10 xguess 9 54];
[cfx, ~]=fit(indX,double(xcs),ft,opts);
opts.StartPoint = [max(zcs) max(zcs)/10 zguess 9 54];
[cfz, ~]=fit(indZ,double(zcs),ft,opts);

if verbose > 2
    % 1D figure
    figure;
    gauss = @(x,a,b,c) a.*exp(-((x-b)/c).^2);
    subplot(1,2,1); hold on;
    plot(indX,xcs, 'xb')
    plot(indX,cfx(indX), '-r')
    plot(indX,gauss(indX,cfx.a1,cfx.b1,cfx.c1),'--k')
    plot(indX,gauss(indX,cfx.a2,cfx.b1,cfx.c2),':k')
    legend('Data','Fit','Gauss1','Gauss2')
    axis tight
    xlabel('pixel X-direction')
    
    subplot(1,2,2); hold on;
    plot(indZ,zcs, 'xb')
    plot(indZ,cfz(indZ), '-r')
    plot(indZ,gauss(indZ,cfz.a1,cfz.b1,cfz.c1),'--k')
    plot(indZ,gauss(indZ,cfz.a2,cfz.b1,cfz.c2),':k')
    axis tight
    xlabel('pixel Z-direction')
end
if verbose > 0
    fprintf('The 1D fit beam center is: (%.2f,%.2f)\n',cfx.b1,cfz.b1)
end

bc=[cfx.b1 cfz.b1];
%Beamcenter can be fitted with 2D gaussians, however, it is a bit slow and
%it only changes the optimized value of the beam center minutely.
fit_2d_beam_center=0;
if fit_2d_beam_center
    % Use 1D results as starting point for 2D fit
    [X,Z]=ndgrid(indX,indZ);
    XZ(:,:,1)=X;
    XZ(:,:,2)=Z;
    
    if verbose > 0
        % Use 1D results as starting point for 2D fit
        fprintf('Determining beam position using 2D functions...\n')
    end
    % Obtain starting guess from 1D-fit.
    xg=cfx.b1;      zg=cfz.b1;
    sigxg1=cfx.c1;   sigzg1=cfz.c1;
    sigxg2=cfx.c2;   sigzg2=cfz.c2;
    % Define 2D functions
    % Parameters are: [x0, y0, Amplitude1, sigmax1, sigmay1, Amplitude2, sigmax2, sigmay2, constant_bg]
    x0 = [xg,zg, 2.6232e+03,sigxg1,sigzg1,300,sigxg2,sigzg2,0]; %Inital guess parameters
    xlower=[min(indX),  min(indZ), 0 , 0 ,0 , 0 , 20 ,20 ,0];
    xupper=[max(indX),max(indZ), Inf, (max(indX)-min(indX)),(max(indZ)-min(indZ)), ...
        Inf, (max(indX)-min(indX)),(max(indZ)-min(indZ)), Inf];
    DoubleGauss2D = @(x, xdata)  x(3)*exp(   -(((xdata(:,:,1)-x(1))/x(4)).^2 + ...
        ((xdata(:,:,2)-x(2))/x(5)).^2) )+ ...
        x(6)*exp(   -(((xdata(:,:,1)-x(1))/x(7)).^2 + ...
        ((xdata(:,:,2)-x(2))/x(8)).^2) )+x(9);
    BC_F=DoubleGauss2D;
    % Fit to 2D data
    [x,resnorm,residual,exitflag] = lsqcurvefit(BC_F,x0,XZ,IPBeam,xlower,xupper);
    
    % Print results
    fprintf('The 2D fit beam center is: (%.2f,%.2f)\n',x(1:2))
    
    % Plot results
    % Limits for plots of beam center
    lim = [x(1)-4*sigxg1 x(1)+4*sigxg1 x(2)-4*sigzg1 x(2)+4*sigzg1];
    if verbose > 2
        figure;
        % Subplot 1
        subplot(1,3,1); hold on;
        surf(X,Z,IPBeam);
        plot3(x(1), x(2), max(max(IPBeam))+1, 'marker', '+')
        shading flat; view([0 90]); title('Raw data'); axis(lim); colorbar;
        % Subplot 2
        subplot(1,3,2); hold on
        surf(X,Z,BC_F(x,XZ))
        plot3(x(1), x(2), max(max(BC_F(x,XZ)))+1, 'marker', '+')
        shading flat; view([0 90]); title(['Fit - ' beamshape],'interpreter','none')
        axis(lim); colorbar;
        % Subplot 3
        subplot(1,3,3); hold on
        surf(X,Z,residual)
        plot3(x(1), x(2), max(max(residual))+1, 'marker', '+')
        shading flat; view([0 90]); title('Residual'); axis(lim); colorbar;
    end
    % Choose beam center for integration
    
    bc=[x(1) x(2)];
end


end

%% Integrate image plate
function int_data=integrate_data(int,mask,settings)
global verbose pixsize
% bring commandwindow to the front so the user sees that integration 
% has started
commandwindow;

number_of_data_files=length(int);
int_data=cell(number_of_data_files,1);

for i=1:number_of_data_files
    number_of_image_plates=length(int{i});
    int_data{i}=cell(number_of_image_plates,1);
    
    for j=1:number_of_image_plates
        
        %relative beam center of this imageplate
        ip_bc=settings.image_plate_beam_center{i}(j,:);
        %coordination transformations
        [pix2tth, tth2pix, tth2mm, ...
         x_pix2mm, x_mm2pix, ...
         z_pix2mm, z_mm2pix ]=coordinate_transformations(ip_bc);  %#ok
        
        %The integration path as function of x and 2theta
        curve_func=get_curve_func(ip_bc,settings.path_shape);
        
        %setup the 2theta axis
        tth_step=pix2tth(2)-pix2tth(1);
        minz=1;
        tth_min=pix2tth(minz);
        if tth_min<0.01
            tth_min=pix2tth(ceil(tth2pix(0.01))); %FIXLATER
        end
        tth_max=pix2tth(size(int{i}{j},2));
        tth_values=tth_min:tth_step:tth_max;
        
        nx=size(int{i}{j},1);
        x=1:nx;
        
        %initiate variables for integrated result
        int_mean=zeros(size(tth_values));
        int_esd=zeros(size(tth_values));
        
        if verbose > 0
            fprintf('Starting curve integration\n')
        end
        integration_method='block_interp2';
        switch integration_method
            case 'block_interp2'
                block_size=1000;
                dx=1;
                x_sampling=(x(1):dx:x(end))';
                n_sampling=length(x_sampling);
                X_sampling=repmat(x_sampling,1,block_size);
                int_temp=[int{i}{j} int{i}{j}(:,end)];
                mask_temp=[mask{i}{j} mask{i}{j}(:,end)];
                %could this loop could be optimized with parfor?
                for k =1:block_size:length(tth_values)
                    % take out a "block" of data. This limits the memory
                    % consumption while retaining the benefits of vectorization.
                    try
                        tth_block=tth_values(k:(k+block_size-1));
                    catch err   %#ok
                        tth_block=tth_values(k:end);
                        block_size=length(tth_block);
                        X_sampling=repmat(x_sampling,1,block_size);
                    end
                    
                    %determine the range needed to contain all curved paths
                    Z_curve_test=curve_func(x_sampling,tth_block(1));
                    Z_min=floor(min(Z_curve_test));
                    Z_curve_test=curve_func(x_sampling,tth_block(end));
                    Z_max=ceil(max(Z_curve_test));
                    
                    if Z_min<=0; Z_min=1; end
                    
                    z_block_ext=Z_min:Z_max;
                    
                    if verbose >2 
                        fprintf('Block length: %.0f\n',length(z_block_ext))
                    end
                    % necessary to include the extra line because numerical errors
                    % might make sampling values slightly too large
                    
                    int_block=int_temp(:,z_block_ext);
                    mask_block=mask_temp(:,z_block_ext);
                    
                    % calculate the curved paths
                    Z_curve=curve_func(X_sampling,repmat(tth_block,n_sampling,1));
                    ind_keep=isreal(Z_curve);
                    Z_curve(~ind_keep)=0;
                    Z_curve=real(Z_curve);
                    
                    %get interpolated intensities along curved paths
                    [Z_block, X_block]=meshgrid(z_block_ext, x);
                    int_curve  =interp2(Z_block,X_block,int_block,...
                        Z_curve,X_sampling,'linear',0);
                    mask_curve =interp2(Z_block,X_block,mask_block,...
                        Z_curve,X_sampling,'linear',0);
%                     int_curve  =interp2(Z_block,X_block,int_block,...
%                         Z_curve,X_sampling,'nearest',0);
%                     mask_curve =interp2(Z_block,X_block,mask_block,...
%                         Z_curve,X_sampling,'nearest',0);
                    
                    if sum(sum(isnan(mask_curve)))>0
                        fprintf(['NaN found during integration of pixel'...
                                 ' %.0f to %.0f\n'],k, k+block_size-1)
                    end
                    
                    % calculate weights for integration
                    switch settings.weighting_method
                        case 'curve_length' 
                            %this looks like method used in old-script
                            dZ=diff(Z_curve,1);
                            dr=sqrt((dZ*pixsize(2)).^2+(dx*pixsize(1)).^2); %mm
                            %make symmetric around the point
                            dra=(dr(2:end,:)+dr(1:end-1,:))/2; 
                            %pad edges of matrix
                            weight=[dra(1,:); dra; dra(end,:)];                             
                        case 'area' %calculate the area of a parallelogram.
                            dZ=diff(Z_curve,1,2);
                            if size(dZ,2)~=1
                                dZa=(dZ(:,2:end)+dZ(:,1:end-1))/2;
                            else
                                dZa=dZ; %very rare occurence
                            end
                            
                            dA=abs((dZa.*pixsize(2)).*(dx.*pixsize(1)));  %dZa is the extra area
                            weight=[dA(:,1) dA dA(:,end)]; %pad edges
                            
                        case 'unity'
                            weight=ones(size(int_curve));
                            
                        otherwise
                            error('Unknown weighting scheme')
                    end
                    
                    %apply mask
                    weight=weight.*mask_curve;
                    % calculate integrated intensity
                    mean_curve=sum(int_curve.*weight,1)./sum(weight,1);
                    %the calculation of esd is problematic when the pixels 
                    %are "over-sampled"
                    if dx==1
%                         esd_curve=sum(((repmat(mean_curve,n_sampling,1)-int_curve).*weight).^2,1)./(sum(weight,1)).^2;
                        esd_curve=sum(((repmat(mean_curve,n_sampling,1)-int_curve).*mask_curve).^2,1)./(sum(mask_curve,1)).^2;
                    else
                        %oversampling for integration
                        error('Do not know what to do, due to oversampling')
                    end
                    
                    int_mean(k:(k+block_size-1))=mean_curve;
                    int_esd(k:(k+block_size-1))=esd_curve;
                end
                
            case 'tine'
                %           x=[1:1:xpmax];  %defines pixels in x-direction
                % for r=minz:1:zpmax;  %integration from chosen z and onwards
                % u(r)=abs(r-bc(1,1)+cor)*pixsize/R*180/pi;     %u is equal to 2theta
                % y=real(sqrt((r+cor-bc(1,1))^2-(x-bc(1,2)).^2)+bc(1,1)-cor);
                % % The summation along the circular arc
                % ds(1)=0;F(1)=0;
                % for n=2:xpmax
                %     F(1)=0;
                %     m1=floor(y(n));
                %     temp=y(n)-m1-0.5;
                % if temp>0 f1=(data(m1,n)+data(m1,n-1))/2;
                %    m2=m1+1; f2=(data(m2,n)+data(m2,n-1))/2;
                %    F(n)=f1+(f2-f1)*temp;
                % else f1=(data(m1,n)+data(m1,n-1))/2;
                %      m2=m1-1;f2=(data(m2,n)+data(m2,n-1))/2;
                %      F(n)=f1-(f2-f1)*temp;
                % end
                %     ds(n)=sqrt((y(n)-y(n-1))^2+1);
                %     Ids(n)=ds(n)*F(n);
                % end
                %     pathds(r)=sum(ds);
                %     Int(r)=sum(Ids);
                %     Int_av(r)=sum(Ids)/(n-1); %Gennemsnit af int.
                %     sig_st_af(r)=sqrt(sum((Ids-Int_av(r)).^2)/(n-2)); %n har n-1 elementer. -1 for at få standard afvigelsen
                %     error_v(r)=sig_st_af(r)/sqrt(n-1);
                %
                % end
                % pathsrel=pathds(end)./pathds; %weights acoording to path lenght
                % Intr=Int.*pathsrel; %Intensity normalized to same solid angle
                % %Intp=Intr./(con*(xpmax-1));  %converts intensity from IP to photons
                % Intp=Intr./(con);  %converts intensity from IP to photons
                % error=sqrt(Intp);
                % error_vc=error_v./con.*pathsrel*(xpmax-1); %xpmax-1: fordi vi skal have for den integrerede int. ikke gennemsnittet
                %
        end
        
        int_data{i}{j}=[tth_values', int_mean', int_esd'];
        
        
        if verbose > 1
            
            int_masked=int_temp;
            int_masked(mask_temp==0)=0;
            plot_and_analyse_2d(int_masked,settings, i,j);
        end
    end
end


end

function [cf, cf_bc]=get_curve_func(bc,path_shape)
global R pixsize



[pix2tth, tth2pix, tth2mm, ...
 x_pix2mm, x_mm2pix,...
 z_pix2mm, z_mm2pix ]=coordinate_transformations(bc); %#ok

switch path_shape
    case 'line'
        cf=@(x_pix,tth) ones(size(x_pix)).*tth2pix(tth);
        cf_bc=@(x_pix,tth,bc) ones(size(x_pix)).*tth2pix(tth);
    case 'circle'
        cf=@(x,tth) z_mm2pix(abs((tth2mm(tth)).^2-(x_pix2mm(x)).^2).^.5);
        cf_bc=@(x,tth,bc) (abs((tth/180*pi*R).^2-((x-bc(1))*pixsize(1)).^2).^.5)/pixsize(2)+bc(2);
    case 'debye_cone'
        r  = @(x_mm) sqrt(R.^2 + x_mm.^2); % mm
        z_mm = @(x_mm,tth) R/180*pi.*atand(tand(tth).*cosd(abs(asind(-x_mm./sind(tth)./r(x_mm)))));
        cf= @(x_pix,tth) z_mm2pix(z_mm(x_pix2mm(x_pix), tth)); %pixel
        
       %CHECK THAT THIS IS CORRECT. PIXELSIZE IN X AND Z DIRECTION
        r_pix  = @(x_pix) sqrt(R.^2 + (x_pix.*pixsize(1)).^2)./pixsize(1); % pix
        cf_bc = @(x_pix,tth,bc) R/180*pi/pixsize(2).*atand(tand(tth).*cosd(abs(asind(-(x_pix-bc(1))./(sind(tth).*r_pix(x_pix-bc(1)))))))+bc(2);
    otherwise
        error('Unknown name for path_shape')
end
end

%% Various coordinate transformations
function [pix2tth, ...
    tth2pix, ...
    tth2mm,...
    x_pix2mm, ...
    x_mm2pix, ...
    z_pix2mm, ...
    z_mm2pix]=coordinate_transformations(bc)

global pixsize R
x_pix2mm = @(x) (x-bc(1))*pixsize(1);
z_pix2mm = @(z) (z-bc(2))*pixsize(2);
x_mm2pix = @(x_mm) x_mm./pixsize(1)+bc(1);
z_mm2pix = @(z_mm) z_mm./pixsize(2)+bc(2);

pix2tth=@(zpix) z_pix2mm(zpix)./R*180/pi;  % from number of pixels to radians
tth2mm= @(tth) tth/180*pi*R;
tth2pix= @(tth) z_mm2pix(tth2mm(tth));


end

%% Post integration analysis of the 2D data
function plot_and_analyse_2d(int_mask,settings,file_number,ip_number)
global dfs

ip_bc=settings.image_plate_beam_center{file_number}(ip_number,:);
% origin=settings.image_plate_origin{file_number}(ip_number,:);

pix2tth=coordinate_transformations(ip_bc);
x_sampling=1:size(int_mask,1);
curve_func=get_curve_func(ip_bc,settings.path_shape);

fig=figure('Outerposition',[ 0  50   dfs],...
    'Name',sprintf('Scan file no. %.0f - Image plate no. %.0f',...
    file_number,ip_number),...
    'NumberTitle','off');

tgroup = uitabgroup('Parent', fig);
tab1 = uitab('Parent', tgroup, 'Title', '2D data');
tab2 = uitab('Parent', tgroup, 'Title', 'Fit results');

ax_pix = axes('Parent',tab1,...
              'units','normalized',...
              'position',[0.1 0.17 0.85 0.75 ],...
              'fontsize',12);
imagesc(log10(int_mask));

set(ax_pix,'YDir','normal',...
    'YLim',[1, size(int_mask,1)],...
    'XLim',[1, size(int_mask,2)],...
    'box','off');
hold on
%make an axis on top with 2theta angle
ax_tth = axes('Parent',tab1,...
    'Position',ax_pix.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'TickDir','out',...
    'box','off',...
    'fontsize',12);
set(ax_tth,'YTick',[])
set(ax_tth,'YColor','w')
set(ax_tth,'Xlim',pix2tth(ax_pix.XLim))

axes(ax_pix); %make pixel-axis active

%synchronize the pixel and 2theta axes.
addlistener( ax_pix, 'XLim', 'PostSet', ...
    @(src,evt) set(ax_tth,'XLim',pix2tth(ax_pix.XLim)) );
addlistener( ax_pix, 'YLim', 'PostSet', ...
    @(src,evt) set(ax_tth,'YLim',ax_pix.YLim));

%handle for curves to be plotted on top of data
ttl=[];

bg1 = uibuttongroup('parent',tab1, 'visible','on','units','pixels',...
    'Position',[5 5 310 40]);
pb1=uicontrol('Style', 'pushbutton', 'Parent',bg1,...
    'String', 'Analyse',...
    'Position', [10 10 50 20],...
    'Callback', {@ui_analyse_peak, int_mask, ip_bc});  %#ok
pb2=uicontrol('Style', 'pushbutton', 'Parent',bg1,...
    'String', 'Hide overlay',...
    'Position', [70 10 70 20],...
    'Callback', @overlay_visibility);
pb3=uicontrol('Style', 'pushbutton', 'Parent',bg1,...
    'String', 'Overlay path',...
    'Position', [150 10 70 20],...
    'Callback', {@overlay_path}); %#ok
pb4=uicontrol('Style', 'pushbutton', 'Parent',bg1,...
    'String', 'Fit beam center',...
    'Position', [230 10 70 20],...
    'Callback', {@fit_path}); %#ok

[h_clim_min, h_clim_max]=add_colorlimit_fields(tab1,ax_pix,[340 5 120 45]);

clim_values=get(ax_pix,'Clim');
set(h_clim_min,'string',num2str(clim_values(1)));
set(h_clim_max,'string',num2str(clim_values(2)));

% Create the column and row names in cell arrays 
cnames = {'2theta','X_BC',...
          'fit-2theta','esd',...
          'fit-X_BC','esd'};
% fit_data=zeros(1,8);
fit_data=[];
fit_table = uitable(tab2,'Data',fit_data,'ColumnWidth',{75},...
    'ColumnName',cnames,...
    'columnformat',{'short','short','short','short','short','short',},... 
    'units','normalized','position',[0.1 0.15 0.85 0.80 ]);

set( gcf, 'toolbar', 'figure' )
do_not_show_again=0; %do not show guide in fit_path-callback
disp('')

    function overlay_path(~,~)
        [z0,~] = ginput(1);
        tth0=pix2tth(z0);
        z0_curve=curve_func(x_sampling,tth0);
        ttl(end+1)=plot(z0_curve,x_sampling,'color',[0.1 0.1 0.1],'linewidth',1);
        
    end

    function fit_path(~,~)
        global R pixsize
        
%         [~, cf_bc]=get_curve_func(ip_bc,'circle');
        [~, cf_bc]=get_curve_func(ip_bc,settings.path_shape);
        cf_fit=@(p,x) cf_bc(x,p(1),[p(2) ip_bc(2)]);
        
        if ~do_not_show_again
            choice=questdlg(['Click on either side of a reflection to '...
                'determine the beam center from the path shape. Only '...
                'the Z-range between the two selected points is used '...
                'for the fitting'],...
                'Fit beam center','OK','Do not show again','OK');
            if strcmp(choice,'Do not show again')
                do_not_show_again=1;
            end
        end
        
        %get 2 points from user. One on either side of the reflection
        [z0,x0] = ginput(1);
        x0=round(x0);
        z0=round(z0);
        tth0=sqrt(((x0-ip_bc(1)).*pixsize(1)).^2+((z0-ip_bc(2)).*pixsize(2)).^2)/R*180/pi;
        
        x1=1:size(int_mask,1);
        z1=round(cf_fit([tth0 ip_bc(1)],x1));
        %determine the position of the peak maximum
        d1=100;
        zlimit=[z1'-d1 z1'+d1];
        
        
        
        
        
        zselect=zeros(length(x1),2*d1+1);
        for i=1:length(x1)
            %ensure that data-selection is within the range of data
            if zlimit(i,1)<1; zlimit(i,1)=1; end
            if zlimit(i,2)>size(int_mask,2); zlimit(i,2)=size(int_mask,2); end
            zselect(i,:)=zlimit(i,1):zlimit(i,2);
        end

        
        d2=20;
        z1=zeros(size(int_mask,1),d2*2+1);
        int_z1=zeros(size(int_mask,1),d2*2+1);
        z_peak=zeros(size(int_mask,1),1);
        for k=1:size(int_mask,1);
            [~,peak_index_i]=max(int_mask(k,zselect(k,:)),[],2);
            z_peak(k)=peak_index_i+zselect(k,1);
        
            z1(k,:)=z_peak(k)-d2:z_peak(k)+d2;
            int_z1(k,:)=int_mask(k,z1(k,:));
        end
        
        z_norm=sum(int_z1,2);
        z_m1=sum(int_z1.*z1,2)./z_norm;

        tth_guess=(z_m1(round(ip_bc(1)))-ip_bc(2))*pixsize(2)/R*180/pi;
        
        ind=~isnan(z_m1);
        x1=x1(ind);
        z_m1=z_m1(ind);
        z_peak=z_peak(ind);
        z_norm=z_norm(ind);
        
        %remove is outside image plate
        ind_outside=(z_norm<mean(z_norm)/2);
        x1=x1(~ind_outside);
        z_m1=z_m1(~ind_outside);
        z_peak=z_peak(~ind_outside);
        

        
        p_guess=[tth_guess ip_bc(1)];
        [p_fit,p_esd, ~]=fit_path_center(x1,z_m1,cf_fit,p_guess);
        ip_bc_fit=[p_fit(2) ip_bc(2)];
        tth_fit=p_fit(1);
        fprintf('Start guess:   beam center_x: %6.1f, 2theta: %6.4f\n',...
            ip_bc(1),tth_guess);
        fprintf(['Refine values: beam center_x: %6.1f (%6.1f),'...
            ' 2theta: %8.4f (%8.4f)\n'],...
            p_fit(2),p_esd(2),p_fit(1),p_esd(1));
        
        fit_data(end+1,:)=[tth_guess, ip_bc(1) , p_fit(1),p_esd(1),...
                           p_fit(2),p_esd(2)];
        set(fit_table,'Data',fit_data)
        
        ttl(end+1)=plot(z_peak,x1,'.r','linewidth',1);
        ttl(end+1)=plot(z_m1,x1,'.b','linewidth',1);
        %plot path with guessed and refined beam center
        
        z0_curve=cf_fit(p_guess,x_sampling);
        ttl(end+1)=plot(z0_curve,x_sampling,...
            '-','color',[0 0 0],'linewidth',1);

        z0_curve=cf_fit(p_fit,x_sampling);
        ttl(end+1)=plot(z0_curve,x_sampling,...
            '-','color',[1 0.1 0.1],'linewidth',1);
        
        legend(ttl(end-3:end),'max value','1^{st} moment', ...
                              'Initial','Fefined')
        
    end

    function [p_fit,esd,resnorm]=fit_path_center(xi,zi,z_pix_circle,p0)
        xi=xi(:); zi=zi(:);
        
        p_lower=[p0(1)-1,p0(2)-500]; p_upper=[p0(1)+2,p0(2)+500];
        options=optimoptions('lsqcurvefit','TolX',1e-9);
        [p_fit,resnorm,resid,~,~,~,J]=lsqcurvefit(z_pix_circle,p0,xi,zi,...
                                                  p_lower,p_upper,options);
        p_conf = nlparci(p_fit,resid,'jacobian',J);
        degrees_of_freedom=length(zi)-length(p0)-1;   %CHECK LATER is it correct
        esd=conf2esd(p_conf,degrees_of_freedom);
    end

    function overlay_visibility(~, ~)
        switch get(ttl(1),'visible')
            case 'off'
                set(ttl,'visible','on')
                set(pb2,'string','Hide overlay')
            case 'on'
                set(ttl,'visible','off')
                set(pb2,'string','Show overlay')
        end
    end

end

function ui_analyse_peak(~,~,data,ip_bc)

nx=size(data,1); nz=size(data,2);
xaxis=1:nx;      zaxis=1:nz;

rect=round(getrect(gcf));
z0=rect(1); x0=rect(2);
zw=rect(3); xw=rect(4);

%abort if the selected region has 0 area
if zw==0 || xw==0
    return
end

xmin=x0; xmax=x0+xw;
zmin=z0; zmax=z0+zw;

if xmin < 1;  xmin=1;  end; if zmin < 1;  zmin=1;  end
if xmax > nx; xmax=nx; end; if zmax > nz; zmax=nz; end

data_select=data(xmin:xmax,zmin:zmax);
xaxis_select=xaxis(xmin:xmax);
zaxis_select=zaxis(zmin:zmax);
analyse_peakshape(data_select,xaxis_select,zaxis_select,ip_bc)
end

% Anlyse the peakshape in a selected area
function analyse_peakshape(int,xaxis,zaxis,ip_bc)
% the input data must contain a section of the image plate with only a
% single image plate
fig1d2d=figure;   %#ok
subplot(2,2,1)
imagesc(zaxis,xaxis,int);
set(gca,'YDir','normal');
title('Radial')
hold on


% locate the diffraction line
x_index=round(size(int,1)/2);
x0=xaxis(x_index);
[~,z_index]=max(int(x_index,:));
z0=zaxis(z_index);
plot(z0,x0,'kx','markersize',20,'linewidth',2)

xbc=ip_bc(1); zbc=ip_bc(2);
curve_x= @(alpha,r) r*sind(alpha)+xbc;
curve_z= @(alpha,r) r*cosd(alpha)+zbc;
curve_r    =@(x,z) sqrt((x-xbc).^2+(z-zbc).^2);
%curve_alpha=@(x,z) atand((x-xbc)./(z-zbc));
r0=curve_r(x0,z0);

%find boudaries of alpha
x_min=xaxis(1);
x_max=xaxis(end);
alpha_xr=@(x,r) atand((x-xbc)./sqrt(r0.^2-(x-xbc).^2));

alpha_range=[alpha_xr(x_min,r0) alpha_xr(x_max,r0)];
alpha_sampling=linspace(alpha_range(1),alpha_range(2),100);
xc0=curve_x(alpha_sampling,r0);
zc0=curve_z(alpha_sampling,r0);

plot(zc0,xc0,'w','linewidth',2)
[X,Z]=meshgrid(xaxis,zaxis);

r_step=1;
wr=100;
r_sample=r0-wr:r_step:r0+wr;
xri=zeros(length(r_sample),length(alpha_sampling));
zri=zeros(length(r_sample),length(alpha_sampling));

for i = 1:length(alpha_sampling)
    xri(:,i)=curve_x(alpha_sampling(i),r_sample);
    zri(:,i)=curve_z(alpha_sampling(i),r_sample);
end
plot(zri(:,1:5:end),xri(:,1:5:end),'w')

% Moments along radial paths
int_radial=interp2(X,Z,int',xri,zri);
ind_nan=isnan(int_radial);
int_radial(ind_nan)=0;

% subtract a background
int_radial_nobg=subtract_bg(r_sample,int_radial);
R_sampling=repmat(r_sample',1,length(alpha_sampling));
norm_factor_radial=sum(int_radial_nobg,1);
r1=sum(int_radial_nobg.*R_sampling,1)./norm_factor_radial;


%redo using the first central moment as estimate of the peakposition

xri=zeros(length(r_sample),length(alpha_sampling));
zri=zeros(length(r_sample),length(alpha_sampling));
R_sampling=zeros(length(r_sample),length(alpha_sampling));

for i = 1:length(alpha_sampling)
    R_sampling(:,i)=r1(i)-wr:r_step:r1(i)+wr;
    xri(:,i)=curve_x(alpha_sampling(i),R_sampling(:,i));
    zri(:,i)=curve_z(alpha_sampling(i),R_sampling(:,i));
end
plot(zri(:,1:5:end),xri(:,1:5:end),'w')

% Moments along radial paths
int_radial=interp2(X,Z,int',xri,zri);
ind_nan=isnan(int_radial);
int_radial(ind_nan)=0;

% subtract a background
int_radial_nobg=subtract_bg(r_sample,int_radial);

R1=repmat(r1,length(r_sample),1);
r2=sum(int_radial_nobg.*(R_sampling-R1).^2,1)./norm_factor_radial;
r3=sum(int_radial_nobg.*(R_sampling-R1).^3,1)./norm_factor_radial;
r4=sum(int_radial_nobg.*(R_sampling-R1).^4,1)./norm_factor_radial;
r5=sum(int_radial_nobg.*(R_sampling-R1).^5,1)./norm_factor_radial;
r7=sum(int_radial_nobg.*(R_sampling-R1).^7,1)./norm_factor_radial;
% The third moment measures skewness, the lack of symmetry, while the
% fourth moment measures kurtosis, the degree to which the 
% distribution is peaked
s_r=sqrt(r2);
skewness_r=r3./s_r.^3;
kurtosis_r=r4./s_r.^4;

%Moments along "linear" paths
subplot(2,2,2)
imagesc(zaxis,xaxis,int);
set(gca,'YDir','normal');
title('Linear')
hold on
plot(z0,x0,'kx','markersize',20,'linewidth',2)

xli=zeros(length(r_sample),length(alpha_sampling));
zli=zeros(length(r_sample),length(alpha_sampling));

r1z=r1;
r1z(isnan(r1))=r0;

for i = 1:length(alpha_sampling)
    xli(:,i)=repmat(curve_x(alpha_sampling(i),r1z(i)),length(r_sample),1);
    zcenter=curve_z(alpha_sampling(i),r1z(i));
    zli(:,i)=(zcenter-wr:r_step:zcenter+wr)';
end
plot(zli(:,1:5:end),xli(:,1:5:end),'w')

int_line=interp2(X,Z,int',xli,zli);
ind_nan=isnan(int_line);
int_line(ind_nan)=0;


% subtract a background
int_line_nobg=subtract_bg(r_sample,int_line);

norm_factor_line=sum(int_line_nobg,1);
z1=sum(int_line_nobg.*(zli),1)./norm_factor_line;

wr=zeros(length(alpha_sampling),6);
wz=zeros(length(alpha_sampling),6);

parr=zeros(length(alpha_sampling),7);
parz=zeros(length(alpha_sampling),7);


for i = 5:length(alpha_sampling)-5
    
    [wr(i,1), wr(i,2)]=get_peak_width(r_sample',int_radial_nobg(:,i),0.5);
    [wr(i,3), wr(i,4)]=get_peak_width(r_sample',int_radial_nobg(:,i),0.05);
    [wr(i,5), wr(i,6)]=get_peak_width(r_sample',int_radial_nobg(:,i),0.01);
    
    sr=spec1d(r_sample',int_radial_nobg(:,i),ones(size(r_sample)));
    [sfr, pars]=fits(sr,'asympseudovoigt',[norm_factor_radial(i) r1(i) wr(i,1) 0 0.001 0 0],[1 1 1 1 1 1 0]);
    parr(i,:)=pars.pvals';
    
    [wz(i,1), wz(i,2)]=get_peak_width(zli(:,i),int_line_nobg(:,i),0.5);
    [wz(i,3), wz(i,4)]=get_peak_width(zli(:,i),int_line_nobg(:,i),0.05);
    [wz(i,5), wz(i,6)]=get_peak_width(zli(:,i),int_line_nobg(:,i),0.01);
    
    sz=spec1d(zli(:,i),int_line_nobg(:,i),ones(size(r_sample)));
    [sfz, pars]=fits(sz,'asympseudovoigt',[norm_factor_line(i) z1(i)  wz(i,1) 0 0.001 0 0],[1 1 1 1 1 1 0]);
    parz(i,:)=pars.pvals';
    
    
end




h=[];
ax3=subplot(2,2,3);
hold on
h(1 ,1)=plot(alpha_sampling,norm_factor_radial,'b');
h(2 ,1)=plot(alpha_sampling,r1,'b');
h(3 ,1)=plot(alpha_sampling,wr(:,1),'b');
h(4 ,1)=plot(alpha_sampling,wr(:,3),'b');
h(5 ,1)=plot(alpha_sampling,wr(:,5),'b');

h(6 ,1)=plot(alpha_sampling,wr(:,2),'b');
h(7 ,1)=plot(alpha_sampling,wr(:,4),'b');
h(8 ,1)=plot(alpha_sampling,wr(:,6),'b');

h(9 ,1)=plot(alpha_sampling,parr(:,4),'b');
h(10 ,1)=plot(alpha_sampling,parr(:,6),'b');


ylabel('-')
xlabel('alpha (degree)')

ax4=subplot(2,2,4);
hold on

h(1 ,2)=plot(alpha_sampling,norm_factor_radial,'b');
h(2 ,2)=plot(alpha_sampling,z1-ip_bc(2),'r');  %subtract the bc to make the
                                               % plot comparable to r1
                                               % plotted in h(2,1)
h(3 ,2)=plot(alpha_sampling,wz(:,1),'b');
h(4 ,2)=plot(alpha_sampling,wz(:,3),'b');
h(5 ,2)=plot(alpha_sampling,wz(:,5),'b');

h(6 ,2)=plot(alpha_sampling,wz(:,2),'b');
h(7 ,2)=plot(alpha_sampling,wz(:,4),'b');
h(8 ,2)=plot(alpha_sampling,wz(:,6),'b');

h(9 ,2)=plot(alpha_sampling,parz(:,4),'b');
h(10 ,2)=plot(alpha_sampling,parz(:,6),'b');






ylabel('-')
xlabel('alpha (degree)')
set(h(:),'visible','off')

linkaxes([ax3, ax4],'xy')
% Create pop-up menu
popup = uicontrol('Style', 'popup',...
    'String', {'norm factor','1st moment',...
    'width 50%' ,'width 5%','width 1%' ,...
    'asym 50%' ,'asym 5%','asym 1%' ,...
    'fit eta','fit asym1'},...
    'Position', [20 5 150 50],...
    'Callback', {@setplot,h});   %#ok


    function Y_nobg=subtract_bg(r,Y)
        Y_bg=Y([1:round(end*0.05) round(end*0.95):end],:);
        r_bg=r([1:round(end*0.05) round(end*0.95):end]);
        A=[r_bg' ones(size(r_bg))'*1000]; %multiplication with 1000 avoids numerical error
        Y_nobg=zeros(size(Y));
        for j = 1:size(Y,2)
            B=Y_bg(:,j);
            x=A\B; %linear least squares
            bg=r'*x(1)+1000*x(2);
            Y_nobg(:,j)=Y(:,j)-bg;
        end
    end
    function setplot(hObj,event,plot_handles) %#ok
        % Called when user activates popup menu
        val = get(hObj,'Value');
        all_str = get(hObj,'String');
        val_str=all_str{val};
        ax=get(plot_handles(1,:),'Parent');
        set(get(ax{1},'Ylabel'),'String',val_str)
        set(get(ax{2},'Ylabel'),'String',val_str)
        set(plot_handles(:),'visible','off')
        set(plot_handles(val,:),'visible','on')
        axis(ax{1},'auto');
        al1=axis(ax{1});
        axis(ax{2},'auto');
        al2=axis(ax{2});
        
        axis([al1(1:2) min([al1(3) al2(3)]) max([al1(4) al2(4)])] )
        
        
    end
end

%% Save integrated data
function save_integrated_data(data,settings,auto_mode)
number_of_data_files=length(data);

%make a "standard output filename" 
[path,name,~]=fileparts(settings.filename{1});
standard_output_file=[path '/' name '_int.dat'];

if isempty(settings.save_result)
    settings.save_result='';
end

if auto_mode
    switch settings.save_result
        case ''
            % Prompt user to input a filename for saving the data.
            choice=questdlg('Do you want to save the integrated data?',...
                'Saving file','yes','no','yes');
            if strcmp(choice,'yes')
                [filename,path] = uiputfile('*.dat','Save file name',...
                    ['./' standard_output_file]);
                if isequal(filename,0)
                    output_file_name=0;
                else
                    output_file_name=[path filename];
                end
            else
                output_file_name=0;
            end
        case '$'
            output_file_name=standard_output_file;
        otherwise
            output_file_name=settings.save_result;
    end
else
    [filename,path] = uiputfile('*.dat','Save file name',...
        ['./' standard_output_file]);
    if isequal(filename,0)
        output_file_name=0;
    else
        output_file_name=[path filename];
    end
end



if isequal(output_file_name,0)
    disp('The integrated intensity was not saved!')
else
    if number_of_data_files==1 && length(data{1})==1
        fid=fopen(output_file_name,'w');
        fprintf(fid,'%16s%16s%16s\r\n','2theta(degree)',...
            'intensity','esd');
        fprintf(fid,'%16.5f%16.6f%16.6f\r\n',data{1}{1}');
        fclose(fid);
    else
        %If there are more than one IP in total, each IP is given a bank
        %number which is increasing from low to high angle
        bank_number=0;
        %Each data bank is saved as a separate file. The number of the 
        % data bank is added to
        % the filename
        output_file_name0=output_file_name;
        [~,~,ext]=fileparts(output_file_name0);
        if isempty(ext)
            %if the name does not have a file type; one will be
            %provided
            ext='.dat';
            output_file_name0=[output_file_name0 ext];
        end
        for i=1:number_of_data_files;
            number_of_image_plates=length(data{i});
            for j=1:number_of_image_plates
                bank_number=bank_number+1;
                output_file_name=regexprep(output_file_name0,ext,...
                    ['_bank' num2str(bank_number) ext]);
                fid=fopen(output_file_name,'w');
                fprintf(fid,'%16s%16s%16s\r\n','2theta(degree)',...
                            'intensity','esd');
                fprintf(fid,'%16.5f%16.6f%16.6f\r\n',data{i}{j}');
                fclose(fid);
            end
        end
    end
end
end

%% Save configuration file
function save_settings_to_file(settings,auto_mode)
if isempty(settings.save_settings)
    settings.save_settings='';
end

if auto_mode
    switch settings.save_settings
        case ''
            choice = questdlg('Do you want to save the settings?',...
                'Saving file','yes','no','yes');
            if strcmp(choice,'yes')
                [filename,path] = uiputfile('*.ini','Save file name');
                if isequal(filename, 0)
                    output_file_name=0;
                else
                    output_file_name=[path filename];
                end
            else
                output_file_name=0;
            end
        case 0
            output_file_name=0;
        otherwise
            output_file_name=settings.save_settings;
    end
else
     [filename,path] = uiputfile('*.ini','Save file name');
    if isequal(filename,0)
        output_file_name=0;
    else
        output_file_name=[path filename];
    end
end

if isequal(output_file_name,0)
    disp('The settings were not saved!')
    settings.save_settings='';
else
    struct2ascii(settings,output_file_name)
    settings.save_settings=output_file_name;
end


end

% 2 functions ascii2struct and struct2ascii. Convert between the MATLAB 
% structure which contains all settings/parameters and an ascii. The ascii 
% file is written as "keyword = blabla". blabla can be a string  
% (optional to incase with '') of array of numbers of a cell array. An 
% array and cell array are written as one would define them in Matlab. An 
% array of numbers does not need to be enclosed with [].

function settings=ascii2struct(inputfile)
fid=fopen(inputfile);
if fid==-1
    error('Settings file (%s) not found',inputfile)
end 
settings_str=fread(fid,'*char')';

settings=[];
set=regexp(settings_str,'(^|\n)(?<!\%)\s*(\w+)\s*\=(.+?)($|[\r\n])','tokens');

%Explanation for regexp expression:
%(^|\n) beginning of string or newline: keyword must be at beginning of file
%(?<!%) must not be preceded by %-sign: comments are ignored
%\s*(\w+)\s*\= maybe some whitespaces, then a keyword, maybe some 
%whitespace then followed by =-sign 

for i=  1:length(set)
    field_name= set{i}{2};
    field_str = set{i}{3};
    
    
    %check number of occurences of keyword. ensure it only occures once
    number_of_occurences=length(regexp(settings_str,['(^|\n)(?<!\%)\s*' field_name '[ ]*\=(.+?)[\r\n^]+'],'match'));
    if number_of_occurences>1
        error('The keyword ''%s'' appeared multiple times in the input file',field_name)
    end
    
    %str2num(field_str) will return empty if field str does not contain a
    %string which can be interpreted as a number; however e.g. 'line' and
    %'rectangle' are names of MATLAB are matlab objects and they will be
    %interpreted as their numeric values
    %The work around below therefore first checks for letters (except e).
    a=regexp(field_str,'({.*})','tokens');
    if ~isempty(a) %is cell-array
        try
            eval(['field_value = ' field_str ';'])
        catch err
            error(fprintf('Failed to read keyword: %s',field_str))
            
        end
    else
        %does it contain letters
        if ~isempty(regexp(field_str,'[a-zA-Z]','once'));
            %could it be a number like '1e-234'
            if ~isempty(regexp(field_str,'^[0-9.,+ ]*[eE][0-9+-]*$','match'));
                %converting to number
                field_value=str2num(field_str);  %#ok
            else
                %saving as string
                field_value=regexprep(field_str,'["'' ]','');
            end
        else
            field_value=str2num(field_str);  %#ok
            if isempty(field_value) %could not convert to number
                %saving as string
                field_value=regexprep(field_str,'["'' ]','');
            end
        end
    end
    
    settings.(field_name)=field_value;
end
end

function struct2ascii(settings,outputfile)


fid=fopen(outputfile,'w');

all_fields=fieldnames(settings);
for i=1:length(all_fields)
    name =all_fields{i};
    value=settings.(name);
    outstr=ex_func(value,6);
    fprintf(fid,'%s = %s \r\n',name,outstr);
    
end

fclose(fid);

    function [ out ] = ex_func( in, prec )
        
        in_datatype = class(in);
        
        switch in_datatype
            case 'char'
                out =  in ;
            case {'single','double','logical'}
                out=matrix2str(double(in),prec);
                out=out(2:end-1); %remove the surrounding square brackets
            case 'cell'
                out=cell2str(in,prec);
            otherwise
                error('Unknown type');
        end
        
        function str=cell2str(in_cell,prec0)
            str='{';
            for j=1:length(in_cell)
                switch class(in_cell{j})
                    case 'double'
                        str=[str  matrix2str(in_cell{j},prec0) ' ,'];   %#ok
                    case 'char'
                        str=[str  '''' in_cell{j} '''' ' ,'];   %#ok
                    case 'cell'
                        str=[str  cell2str(in_cell{j},prec0) ' ,'];  %#ok
                end
            end
            
            str=str(1:end-1); %remove the last comma
            str=[str '}']; %close with curly bracket
            
        end
        
        function str=matrix2str(in_mat,prec0)
            str='[';
            for j=1:size(in_mat,1)
                str =[str num2str(in_mat(j,:),prec0) ';'];   %#ok
            end
            if length(str) > 1
                str=str(1:end-1); %remove the last ";"
            end
            str=[str ']']; %close
        end
        
        
    end
end

%% Small helper functions 
function esd=conf2esd(con,degrees_of_freedom)
width=abs(con(:,2)-con(:,1));
if nargin~=2
    degrees_of_freedom=Inf;
    disp('Assumes infinite degrees of freedom')
end
esd=width/(2*tinv(0.975,degrees_of_freedom));

end

function par_guess=make_peak_guess(tth,int)
%Estimate the intensity, position, width and background for a single
%gaussian on a constant backgroun
[position_guess, int_guess, ~]=get_peak_centre(tth,int);
FWHM_guess=get_peak_width(tth,int,0.5);
width_guess=FWHM_guess/2.3548;
%Background is estimated as the average of the endpoints
bg_guess=mean(int([1 end]));
par_guess=[int_guess, position_guess, width_guess, bg_guess];


end

function [FW, asym]=get_peak_width(x,y,fraction)
[~,yc, indc]=get_peak_centre(x,y);
cog=sum(x.*y)./sum(y);
%left side of peak
xl=x(1:indc-1); yl=y(1:indc-1);
y0=yc*fraction;
[~,indn]=minn(abs(yl-y0),2);
xn=xl(indn);
yn=yl(indn);
[xn, ind]=sort(xn);
yn=yn(ind);
xlf=(y0-yn(1))/(yn(2)-yn(1))*(xn(2)-xn(1))+xn(1);

%right side of peak
xr=x(indc+1:end); yr=y(indc+1:end);
[~,indn]=minn(abs(yr-y0),2);
xn=xr(indn);
yn=yr(indn);
[xn, ind]=sort(xn);
yn=yn(ind);
xrf=(y0-yn(1))/(yn(2)-yn(1))*(xn(2)-xn(1))+xn(1);
FW=xrf-xlf;

asym=((cog-xrf)-(xlf-cog))./FW; %approaches 0 for the symmetric peak
        

end

function [xc, yc, ind]=get_peak_centre(x,y)
[yc, ind]=max(y);
xc=x(ind);
end

function [xm, ind]=minn(x,n)
%get the n lowest numbers
x=x(:);
if length(x)<n
    error(['Number of requested values exceeds the number ' ...
        'of elements in input'])
end
xm=zeros(n,1);
ind=zeros(n,1);
for i=1:n
    [xm(i), ind(i)]=min(x);
    x(ind(i))=NaN;
end
end

function [xm, ind]=maxn(x,n)
%get the n highest numbers
x=x(:);
if length(x)<n
    error(['Number of requested values exceeds the number ' ...
        'of elements in input'])
end
xm=zeros(n,1);
ind=zeros(n,1);
for i=1:n
    [xm(i), ind(i)]=max(x);
    x(ind(i))=NaN;
end
end

function [t, xt1]=vector_lines_intersect(x1,r1,x2,r2)
        %line is defined by a point x1 and a direction r1
        % x(t)=x1+r1*t
        %this function returns the intersection between the two points
        %ensure all input are columns
        x1=x1(:); r1=r1(:);
        x2=x2(:); r2=r2(:);
        
        % Setup as a matrix equation Y=A*t
        Y=x2-x1;       
        A=[r1 -r2];
        %solve for t
        t=A\Y;
        
        %calculate the coordinates.  
        xt1=x1+r1*t(1);
        % xt2=x2+r2*t(2); %xt1 and xt2 should of course be identical
end

function [x0,z0]=get_bottom_right_corner(x,z)
%given the coordinates of an approximately rectangular shape,
%return the coordinates of the "bottom right corner"
%get two bottom points
[~, ind]=minn(x,2);
xb=x(ind); zb=z(ind);
%selct right-most point
[z0, ind]=max(zb);
x0=xb(ind);
end

function [x0,z0]=get_bottom_left_corner(x,z)
% same as above but returns the "bottom left corner"
[~, ind]=minn(x,2);
xb=x(ind); zb=z(ind);
[z0, ind]=min(zb);
x0=xb(ind);
end

function [x0,z0]=get_top_left_corner(x,z)
% same as above but returns the "bottom left corner"
[~, ind]=maxn(x,2);
xt=x(ind); zt=z(ind);
[z0, ind]=min(zt);
x0=xt(ind);
end

function [x0,z0]=get_top_right_corner(x,z)
% same as above but returns the "bottom left corner"
[~, ind]=maxn(x,2);
xt=x(ind); zt=z(ind);
[z0, ind]=max(zt);
x0=xt(ind);
end

function alpha=get_angle(r1,r2)
% alpha=acosd(dot(r1,r2)./(norm(r1)*norm(r2)));
alpha=asind((r1(1)*r2(2)-r1(2)*r2(1))./(norm(r1)*norm(r2)));
end

function out=convert2floatingpoint(in)
global precision
switch precision
    case 'single'
        out=single(in);
    case 'double'
        out=double(in);
    otherwise
        error('Unknown precision')
end

end


