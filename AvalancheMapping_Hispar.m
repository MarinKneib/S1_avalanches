%% Initialize
clear 
close all
 
% This script calls functions from topotoolbox and ImGRAFT:
% Aslak Grinsted (2024). ImGRAFT (https://github.com/grinsted/ImGRAFT), GitHub. Retrieved February 26, 2024.
% Schwanghart, Wolfgang, and Nikolaus J. Kuhn. “TopoToolbox: A Set of Matlab Functions for Topographic Analysis.” Environmental Modelling & Software, vol. 25, no. 6, Elsevier BV, June 2010, pp. 770–81, doi:10.1016/j.envsoft.2009.12.002.

DEM_path = '\path\to\AW3D_DEM\';
RGB_path = '\path\to\RGB\image\';
Brightness_path = '\path\to\brightness\';
Result_path = '\path\to\results\folder\';

%% Load raster data & resample as needed
[DEM,xdem,ydem,Idem] = geoimread([DEM_path,'Hispar_UTM43N.tif']);
[Xdem,Ydem] = meshgrid(xdem,ydem);
dgeot = geotiffinfo([DEM_path,'Hispar_UTM43N.tif']);
pixsz_dem = abs((xdem(2)-xdem(1))*(ydem(2)-ydem(1)));

[Im1,xdm1,ydm1,Idm1] = geoimread([RGB_path,'GEE_ASC_20170921_20171003.tif']);
[Xdm1,Ydm1] = meshgrid(xdm1,ydm1);
dgeot_dm1 = geotiffinfo([RGB_path,'GEE_ASC_20170921_20171003.tif']);
pixsz = abs((xdm1(2)-xdm1(1))*(ydm1(2)-ydm1(1)));

DEM_resamp = interp2(Xdem,Ydem,DEM,Xdm1,Ydm1,'cubic');

%% Slope calculation
DEM_gridobj = GRIDobj(xdm1,ydm1,DEM_resamp);
Slope = gradient8(DEM_gridobj,'deg');

%% Extract time information from geotiffs
ASC_all = dir([RGB_path,'*ASC*']);
ASC_names = {ASC_all(:).name};

DESC_all = dir([RGB_path,'*DESC*']);
DESC_names = {DESC_all(:).name};

ASC_time1 = [];
ASC_time2 = [];
DESC_time1 = [];
DESC_time2 = [];

for ii = 1:length(ASC_all)
    folder = ASC_names{ii};
    ASC_time1 = [ASC_time1;folder(9:16)];
    ASC_time2 = [ASC_time2;folder(18:25)];
end
for ii = 1:length(DESC_all)
    folder = DESC_names{ii};
    DESC_time1 = [DESC_time1;folder(10:17)];
    DESC_time2 = [DESC_time2;folder(19:26)];
end

ASC_time1_datetime = datetime(ASC_time1,'InputFormat','yyyyMMdd');
ASC_time2_datetime = datetime(ASC_time2,'InputFormat','yyyyMMdd');
DESC_time1_datetime = datetime(DESC_time1,'InputFormat','yyyyMMdd');
DESC_time2_datetime = datetime(DESC_time2,'InputFormat','yyyyMMdd');


%% On a bi-monthly basis, extract information of avalanche deposits VS non avalanche deposits


for year = 2017:2021
    % Initialize
    year1_str = string(year);
    year2_str = string(year+1);
    time_edges = datetime(char([year1_str{:},'-11-01'],[year2_str{:},'-05-01'],[year2_str{:},'-11-01']),'InputFormat','yyyy-MM-dd');
    
    % Define the threshold parameter values (different for orbits/seasons)
    threshS_range_desc = [0.3 0.32];
    threshV_range_desc = [0.65 0.4];
    threshO_range_desc = [0.04 0.04];
    thresh22_1_desc = [0.03 0.11];
    thresh22_2_desc = [0.06 0];
    thresh21_desc = [0.33 0.27];
    
    threshS_range_asc = [0.30 0.38];
    threshV_range_asc = [0.6 0.42];
    threshO_range_asc = [0.04 0.06];
    thresh22_1_asc = [0.06 0.13];
    thresh22_2_asc = [0.04 0.01];
    thresh21_asc = [0.38 0.26];
    
    % Ascending
    
    % load brightness
    [Brightness,xb,yb,Ib] = geoimread([Brightness_path,'Brightness_ASC.tif']);
    [Xb,Yb] = meshgrid(xb,yb);
    Brightness = interp2(Xb,Yb,Brightness,Xdm1,Ydm1,"linear");
    
    % parameters for region growing
    SE = strel('square',15);
    
    %parameters for smoothing
    filtSigma = 11;
    filtWidth = 2*ceil(2*filtSigma)+1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    imageFilter2=fspecial('gaussian',17,4);
    
    for jj = 1:length(time_edges)-1
        disp(jj)
        % select corresponding dates
        time_select1 = ASC_time1_datetime(ASC_time1_datetime>=time_edges(jj) & ASC_time1_datetime<time_edges(jj+1));
        time_select2 = ASC_time2_datetime(ASC_time1_datetime>=time_edges(jj) & ASC_time1_datetime<time_edges(jj+1));
    
        % mapping across all scenes
        for ii = 1:length(time_select1)
            
            % rewrite time
            t1string1 = string(datetime(time_select1(ii),'Format','yyyyMMdd'));
            t2string1 = string(datetime(time_select2(ii),'Format','yyyyMMdd'));
    
            % load RGB files 
            RGB1_name = [RGB_path,'GEE_ASC_',t1string1{:},'_',t2string1{:},'.tif'];
            [RGB1,x1,y1,~] = geoimread(RGB1_name);
            [X1,Y1] = meshgrid(x1,y1);
            RGB1_1 = interp2(X1,Y1,RGB1(:,:,1),Xdm1,Ydm1,"linear");
            RGB1_2 = interp2(X1,Y1,RGB1(:,:,2),Xdm1,Ydm1,"linear");
            RGB1 = cat(3,cat(3,RGB1_1,RGB1_2),RGB1_1);
    
            % convert to HSV
            HSV = rgb2hsv(RGB1);
            HSV(:,:,1) = HSV(:,:,1).*HSV(:,:,2)./HSV(:,:,2); % remove NaNs
    
            % remove image edge
            Edge = HSV(:,:,2)==1 & HSV(:,:,3)==1;
            Edge = imdilate(Edge,strel('square',30)); % add buffer
    
            % Normalize Saturation & brightness with mean saturation & brightness
            % of early november scene
            HSV(:,:,2) = HSV(:,:,2)*0.16/nanmean(nanmean(HSV(:,:,2)));
            HSV(:,:,3) = HSV(:,:,3)*0.56/nanmean(nanmean(HSV(:,:,3)));
    
            % remove steep slopes, small patches 
            if sum(sum(Edge>1000))
                Green_zones = HSV(:,:,1)<0.4 & HSV(:,:,2)>0.15 & HSV(:,:,3)>0.2 & Slope.Z<35 & Brightness<0.82 & Brightness>0.1 & ~Edge;
            else
                Green_zones = HSV(:,:,1)<0.4 & HSV(:,:,2)>0.15 & HSV(:,:,3)>0.2 & Slope.Z<35 & Brightness<0.82 & Brightness>0.1;
            end
            S = HSV(:,:,2);
            V = HSV(:,:,3);
            S(~Green_zones) = NaN;
            V(~Green_zones) = NaN;
    
            % Smooth images
            Img1 = RGB1(:,:,1);
            Img2 = RGB1(:,:,2);
            Img1(Brightness<0.1 | Brightness>0.82 | Slope.Z>35) = NaN;
            Img2(Brightness<0.1 | Brightness>0.82 | Slope.Z>35) = NaN;
            Smoothed1 = nanconv(Img1,imageFilter, 'nanout');
            Smoothed2 = nanconv(Img2,imageFilter, 'nanout');
            Smoothed22 = nanconv(Img2-Smoothed2,imageFilter2,'nanout');
            Smoothed21 = nanconv(Img2-Smoothed1,imageFilter2,'nanout');
    
            % threshold values
            threshS = threshS_range_asc(jj);
            threshV = threshV_range_asc(jj);
            threshO = threshO_range_asc(jj);
            thresh22_1 = thresh22_1_asc(jj);
            thresh22_2 = thresh22_2_asc(jj);
            thresh21 = thresh21_asc(jj);
    
            % Use the threshold values to classify raster data
            Classification = S>threshS & V>threshV;
    
            % remove patches smaller than 10 pixels
            Classification = bwareaopen(Classification,10,8); 
    
            % grow region withing given buffer
            Classification_buf = imdilate(Classification,SE);
            Classification = Classification_buf &...
                S>threshS-threshO & V>threshV-threshO &...
                Slope.Z<35 & Brightness<0.82 & Brightness>0.1;
    
            % remove patches smaller than 40 pixels
            Classification = bwareaopen(Classification,40,8);
    
            % additional filter based on smoothed backscatter images
            if jj == 1
                Classification = Classification & ((Img2-Smoothed2>thresh22_1) | ...
                        ((Img2-Smoothed2>thresh22_2) & (Img2-Smoothed1>thresh21)));
            else
                Classification = Classification & ((Smoothed22>thresh22_1) |...
                    ((Smoothed22>thresh22_2) & (Smoothed21>thresh21)));
            end
    
            % remove patches smaller than 40 pixels
            Classification = bwareaopen(Classification,40,8);
    
            % fill holes
            Classification = imfill(Classification,'holes');
    
            geotiffwrite([Result_path,'ASC_',t1string1{:},'-',t2string1{:},'.tif'],Classification,...
                dgeot_dm1.RefMatrix,'CoordRefSysCode',32643)
        end
    
    end
    
    
    
    
    
    
    % Descending
    
    % load brightness
    [Brightness,xb,yb,Ib] = geoimread([Brightness_path,'Brightness_DESC.tif']);
    [Xb,Yb] = meshgrid(xb,yb);
    Brightness = interp2(Xb,Yb,Brightness,Xdm1,Ydm1,"linear");
    
    % parameters for region growing
    SE = strel('square',15);
    
    %parameters for smoothing
    filtSigma = 11;
    filtWidth = 2*ceil(2*filtSigma)+1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    imageFilter2=fspecial('gaussian',17,4);
    
    for jj = 1:length(time_edges)-1
        disp(jj)
        % select corresponding dates
        time_select1 = DESC_time1_datetime(DESC_time1_datetime>=time_edges(jj) & DESC_time1_datetime<time_edges(jj+1));
        time_select2 = DESC_time2_datetime(DESC_time1_datetime>=time_edges(jj) & DESC_time1_datetime<time_edges(jj+1));
    
        % mapping across all scenes
        for ii = 1:length(time_select1)
            
            % rewrite time
            t1string1 = string(datetime(time_select1(ii),'Format','yyyyMMdd'));
            t2string1 = string(datetime(time_select2(ii),'Format','yyyyMMdd'));
    
            % load RGB files 
            RGB1_name = [RGB_path,'GEE_DESC_',t1string1{:},'_',t2string1{:},'.tif'];
            [RGB1,x1,y1,~] = geoimread(RGB1_name);
            [X1,Y1] = meshgrid(x1,y1);
            RGB1_1 = interp2(X1,Y1,RGB1(:,:,1),Xdm1,Ydm1,"linear");
            RGB1_2 = interp2(X1,Y1,RGB1(:,:,2),Xdm1,Ydm1,"linear");
            RGB1 = cat(3,cat(3,RGB1_1,RGB1_2),RGB1_1);
    
            % convert to HSV
            HSV = rgb2hsv(RGB1);
            HSV(:,:,1) = HSV(:,:,1).*HSV(:,:,2)./HSV(:,:,2); % remove NaNs
    
            % remove image edge
            Edge = HSV(:,:,2)==1 & HSV(:,:,3)==1;
            Edge = imdilate(Edge,strel('square',30)); % add buffer
    
            % Normalize Saturation & brightness with mean saturation & brightness
            % of early november scene
            HSV(:,:,2) = HSV(:,:,2)*0.16/nanmean(nanmean(HSV(:,:,2)));
            HSV(:,:,3) = HSV(:,:,3)*0.56/nanmean(nanmean(HSV(:,:,3)));
    
            % remove steep slopes, small patches 
            if sum(sum(Edge>1000))
                Green_zones = HSV(:,:,1)<0.4 & HSV(:,:,2)>0.15 & HSV(:,:,3)>0.2 & Slope.Z<35 & Brightness<0.82 & Brightness>0.1 & ~Edge;
            else
                Green_zones = HSV(:,:,1)<0.4 & HSV(:,:,2)>0.15 & HSV(:,:,3)>0.2 & Slope.Z<35 & Brightness<0.82 & Brightness>0.1;
            end
            S = HSV(:,:,2);
            V = HSV(:,:,3);
            S(~Green_zones) = NaN;
            V(~Green_zones) = NaN;
    
            % Smooth images
            Img1 = RGB1(:,:,1);
            Img2 = RGB1(:,:,2);
            Img1(Brightness<0.1 | Brightness>0.82 | Slope.Z>35) = NaN;
            Img2(Brightness<0.1 | Brightness>0.82 | Slope.Z>35) = NaN;
            Smoothed1 = nanconv(Img1,imageFilter, 'nanout');
            Smoothed2 = nanconv(Img2,imageFilter, 'nanout');
            Smoothed22 = nanconv(Img2-Smoothed2,imageFilter2,'nanout');
            Smoothed21 = nanconv(Img2-Smoothed1,imageFilter2,'nanout');
    
            % threshold values
            threshS = threshS_range_desc(jj);
            threshV = threshV_range_desc(jj);
            threshO = threshO_range_desc(jj);
            thresh22_1 = thresh22_1_desc(jj);
            thresh22_2 = thresh22_2_desc(jj);
            thresh21 = thresh21_desc(jj);
    
            % Use the threshold values to classify raster data
            Classification = S>threshS & V>threshV;
    
            % remove patches smaller than 10 pixels
            Classification = bwareaopen(Classification,10,8); 
    
            % grow region withing given buffer
            Classification_buf = imdilate(Classification,SE);
            Classification = Classification_buf &...
                S>threshS-threshO & V>threshV-threshO &...
                Slope.Z<35 & Brightness<0.82 & Brightness>0.1;

            % remove patches smaller than 40 pixels
            Classification = bwareaopen(Classification,40,8);
    
            % additional filter based on smoothed backscatter images
            if jj == 1
                Classification = Classification & ((Img2-Smoothed2>thresh22_1) | ...
                        ((Img2-Smoothed2>thresh22_2) & (Img2-Smoothed1>thresh21)));
            else
                Classification = Classification & ((Smoothed22>thresh22_1) |...
                    ((Smoothed22>thresh22_2) & (Smoothed21>thresh21)));
            end
    
            % remove patches smaller than 40 pixels
            Classification = bwareaopen(Classification,40,8);
    
            % fill holes
            Classification = imfill(Classification,'holes');
    
    
            geotiffwrite([Result_path,'DESC_',t1string1{:},'-',t2string1{:},'.tif'],Classification,...
                dgeot_dm1.RefMatrix,'CoordRefSysCode',32643)
        end
    
    end
    
end
