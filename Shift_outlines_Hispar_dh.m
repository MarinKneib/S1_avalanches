%% Initialize
clear
close all

% This script calls functions from topotoolbox and ImGRAFT:
% Aslak Grinsted (2024). ImGRAFT (https://github.com/grinsted/ImGRAFT), GitHub. Retrieved February 26, 2024.
% Schwanghart, Wolfgang, and Nikolaus J. Kuhn. “TopoToolbox: A Set of Matlab Functions for Topographic Analysis.” Environmental Modelling & Software, vol. 25, no. 6, Elsevier BV, June 2010, pp. 770–81, doi:10.1016/j.envsoft.2009.12.002.

outline_path = 'path/to/outlines/';
angle_path = 'path/to/incidence/angles/';
dh_path = 'path/to/dh/Hugonnet2021/';
results_path = 'path/to/results/';
RGB_path = 'path/to/RGB/image/';

%% Define parameters

% satellite heading angles (to be modified)
heading_ASC = 349; % in degrees 
heading_DESC = 189;  % in degrees


%% Load data

% load Hugonnet data (m/yr)
[dh,x,y,~] = geoimread([dh_path,'N35E074_2000-01-01_2020-01-01_dhdt.tif']);
[X,Y] = meshgrid(x,y);

% filter abnormal values
dh(dh<-1000) = NaN;
dh(dh<-5) = -5;
dh(dh>5) = 5;

% reference raster (wanted extents)
[Raster,xr,yr,Ir] = geoimread([RGB_path,'image_name.tif']);
dgeot = geotiffinfo([RGB_path,'image_name.tif']);
[Xr,Yr] = meshgrid(xr,yr);

% load tif files with projected local incidence angle (obtained from
% pre-processing of initial images with ESA SNAP software)
[angle_ASC,xa,ya,~] = geoimread([angle_path,'ProjLocIncAngle_asc.tif']);
[angle_DESC,xd,yd,~] = geoimread([angle_path,'ProjLocIncAngle_desc.tif']);

[Xa,Ya] = meshgrid(xa,ya);
[Xd,Yd] = meshgrid(xd,yd);

% resample everything
dh = interp2(X,Y,dh,Xr,Yr,'bilinear');
angle_ASC = interp2(Xa,Ya,angle_ASC,Xr,Yr,'bilinear');
angle_DESC = interp2(Xd,Yd,angle_DESC,Xr,Yr,'bilinear');


%% Load & sort shapefiles by date

% find name of all ASC/DESC outlines
ASC_all = dir([outline_path,'*ASC*.shp']);
ASC_names = {ASC_all(:).name};

DESC_all = dir([outline_path,'*DESC*.shp']);
DESC_names = {DESC_all(:).name};

% initialize time
ASC_time1 = [];
ASC_time2 = [];
DESC_time1 = [];
DESC_time2 = [];

% extract date from image name
for ii = 1:length(ASC_all)
    folder = ASC_names{ii};
    ASC_time1 = [ASC_time1;folder(5:12)];
    ASC_time2 = [ASC_time2;folder(14:21)];
end
for ii = 1:length(DESC_all)
    folder = DESC_names{ii};
    DESC_time1 = [DESC_time1;folder(6:13)];
    DESC_time2 = [DESC_time2;folder(15:22)];
end

% convert time to datetime
ASC_time1_datetime = datetime(ASC_time1,'InputFormat','yyyyMMdd');
ASC_time2_datetime = datetime(ASC_time2,'InputFormat','yyyyMMdd');
DESC_time1_datetime = datetime(DESC_time1,'InputFormat','yyyyMMdd');
DESC_time2_datetime = datetime(DESC_time2,'InputFormat','yyyyMMdd');


%% For each date, shift outlines based on dh & angles

% remove NaNs in dh
dh(isnan(dh)) = 0;

% ASCENDING

% distance to shift (m/yr)
shift_dist = dh./tan(angle_ASC*pi/180);

for ii = 1:length(ASC_time1_datetime)
    disp(ii*100/length(ASC_time1_datetime))
    % load shapefile
    SHP = shaperead([outline_path,ASC_names{ii}]);
    SHP_upd = SHP;

    % convert m/yr to m
    tdiff = year(ASC_time1_datetime(ii))-2000; % Sentinel-1 images were pre-processed using 30m DEM from SRTM mission
    shift_dist_temp = shift_dist*tdiff;

    % for each polygon shift coordinates
    for jj = 1:length(SHP_upd)
        % extract local shift
        shift_local = nanmean(interp2(Xr,Yr,shift_dist_temp,SHP_upd(jj).X,SHP_upd(jj).Y,'nearest'));
        if ~isnan(shift_local)
            SHP_upd(jj).X = double(SHP_upd(jj).X+shift_local*cos((360-heading_ASC)*pi/180));
            SHP_upd(jj).Y = double(SHP_upd(jj).Y+shift_local*sin((360-heading_ASC)*pi/180));
        end
    end

    % export shapefile
    if ~isempty(SHP_upd)
        shapewrite(SHP_upd,[results_path,'dh_correct_',ASC_names{ii}]);
    end
end



% DESCENDING

% distance to shift (m/yr)
shift_dist = dh./tan(angle_DESC*pi/180);

for ii = 1:length(DESC_time1_datetime)
    disp(ii*100/length(DESC_time1_datetime))
    % load shapefile
    SHP = shaperead([outline_path,DESC_names{ii}]);
    SHP_upd = SHP;

    % convert m/yr to m
    tdiff = year(DESC_time1_datetime(ii))-2000; % from SRTM mission
    shift_dist_temp = shift_dist*tdiff;

    % for each polygon shift coordinates
    for jj = 1:length(SHP_upd)
        % extract local shift
        shift_local = nanmean(interp2(Xr,Yr,shift_dist_temp,SHP_upd(jj).X,SHP_upd(jj).Y,'nearest'));
        SHP_upd(jj).X = double(SHP_upd(jj).X+shift_local*cos((360-heading_DESC)*pi/180));
        SHP_upd(jj).Y = double(SHP_upd(jj).Y+shift_local*sin((360-heading_DESC)*pi/180));
    end

    % export shapefile
    if ~isempty(SHP_upd)
        shapewrite(SHP_upd,[results_path,'dh_correct_',DESC_names{ii}]);
    end
end







