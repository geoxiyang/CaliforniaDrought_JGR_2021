%% Regrid and resample climatological data
%  This is the first step of the JGR paper
%  It regrids the PRISM data to the native resolution of MODIS data (MOD13A1.006_500m_aid0001.nc)

clear variables
clc

%% 1. load EVI data to use the lat/lon grid; here just need the lat lon

load('/Volumes/XiYangBackUp/Projects/6.CalDrought/CalEVI_2000_2019_500m_LConly.mat','LC')                        
modis_file                 = '/Volumes/XiYangBackUp/Data/2.SatelliteData/2.MODIS/MOD13A1.006_500m_aid0001.nc';
modlat                     = ncread(modis_file,'lat');
modlon                     = ncread(modis_file,'lon');

%% 2. PRISM datasets
%  all of them are zipped files, we need to unzip them

%  Precipitation (mm), temperature(mean), VPD max and min
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/PRISM.mat','tmean','ppt','vpdmax','vpdmin')

%  time
filenames =dir('/Volumes/XiYangBackUp/Data/28.PRISM/monthlystack/ppt/*_20*_bil.bil');
%  2/16/2020 change: move the first five names to the bottom -- they need
%  to be ordered by time
filenames1 = [filenames(6:end);filenames(1:5)];
yearnumber = nan(length(filenames1),1);
monthnumber= nan(length(filenames1),1);
for ii = 1:length(filenames)
   rawname          = filenames1(ii);
   tempname1        = strsplit(rawname.name,'_');
   tempname2        = tempname1{5};
   yearnumber(ii)   = str2double(tempname2(1:4));
   monthnumber(ii)  = str2double(tempname2(5:6));
end

%  latitude, longitude
NROWS       = 621;
NCOLS       = 1405;
ULXMAP      = -125.000000000000;
ULYMAP      = 49.916666666666;
XDIM        = .041666666667;
YDIM        = .041666666667;

latarray    = linspace(ULYMAP - NROWS*YDIM, ULYMAP,  NROWS);
lonarray    = linspace(ULXMAP, ULXMAP + NCOLS*XDIM,  NCOLS);

%% 3. Regridding

% modlat                  = lat(latroi_y);
% modlon                  = lon(lonroi_x);

[latmesh,lonmesh]       = meshgrid(latarray',lonarray');
[modlatmesh,modlonmesh] = meshgrid(modlat,modlon);
%[modlatmesh,modlonmesh] = meshgrid(lat,lon);

% Note hte regridded dataset is lon lat time
pptregrid               = nan(length(modlon),length(modlat),length(filenames1));
tmeanregrid             = nan(length(modlon),length(modlat),length(filenames1));
vpdregrid               = nan(length(modlon),length(modlat),length(filenames1));

for ii = 1:length(filenames1)
    
    pptfunction         = griddedInterpolant(lonmesh,latmesh,squeeze(ppt(:,:,ii)),'nearest');
    tmeanfunction       = griddedInterpolant(lonmesh,latmesh,squeeze(tmean(:,:,ii)),'nearest');
    vpdfunction         = griddedInterpolant(lonmesh,latmesh,squeeze(0.5*(vpdmax(:,:,ii)+vpdmin(:,:,ii))),'nearest');
    
    pptregrid(:,:,ii)   = pptfunction(modlonmesh,modlatmesh);
    tmeanregrid(:,:,ii) = tmeanfunction(modlonmesh,modlatmesh);
    vpdregrid(:,:,ii)   = vpdfunction(modlonmesh,modlatmesh);
    
end


save('/Volumes/XiYangBackUp/Projects/6.CalDrought/PRISMregrid_500m.mat')
























