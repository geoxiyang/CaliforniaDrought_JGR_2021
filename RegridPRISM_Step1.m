%% Regrid and resample climatological data
%  This is the first step of the JGR paper
%  It regrids the PRISM data to the native resolution of MODIS data (MOD13A1.006_500m_aid0001.nc)

clear variables
clc

%% 1. load EVI data to use the lat/lon grid; here just need the lat lon

modis_file                 = '/Volumes/XiYangBackUp/Data/2.SatelliteData/2.MODIS/MOD13A1.006_500m_aid0001.nc';
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/CalEVI_2000_2019_500m_LConly.mat','LC')                        
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

% figure
% calboundary         = shaperead('usastatehi', 'UseGeoCoords', true,...
%                                 'Selector',{@(name) strcmpi(name,'California'), 'Name'});
% [XX,YY]             = meshgrid(double(latarray),double(lonarray));
% in                  = inpolygon(XX,YY,calboundary.Lat,calboundary.Lon); 
% % === Mapping capacity===
% usamap('California')
% geoshow(calboundary,'DisplayType','polygon','FaceAlpha',0)
% 
% PPTmap              = squeeze(ppt(:,:,1)); %EVIannMean(:,:,1)
% %PPTmap(PPTmap==-9999)=NaN;
% PPTmap(~in)         = NaN;
% h1                  = pcolorm(latarray,lonarray,PPTmap');
% uistack(h1,'bottom')
% colormap(brewermap([],'Blues'))
% %  === Mapping capacity===


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


% %  === Mapping===
% figure
% subplot(2,1,1)
% calboundary         = shaperead('usastatehi', 'UseGeoCoords', true,...
%             'Selector',{@(name) strcmpi(name,'California'), 'Name'});
% usamap('California')
% geoshow(calboundary,'DisplayType','polygon','FaceAlpha',0)
% %in                  = inpolygon(modlatmesh,modlonmesh,calboundary.Lat,calboundary.Lon); 
% h1                  = pcolorm(modlatmesh,modlonmesh,squeeze(pptregrid(:,:,1))); %squeeze(tmeanregrid(:,:,1))
% uistack(h1,'bottom')
% caxis([0,100])
% colormap(brewermap([],'Blues'))
% colorbar
% subplot(2,1,2)
% calboundary         = shaperead('usastatehi', 'UseGeoCoords', true,...
%             'Selector',{@(name) strcmpi(name,'California'), 'Name'});
% usamap('California')
% geoshow(calboundary,'DisplayType','polygon','FaceAlpha',0)
% %in                  = inpolygon(modlatmesh,modlonmesh,calboundary.Lat,calboundary.Lon); 
% h1                  = pcolorm(latmesh,lonmesh,squeeze(ppt(:,:,1))); %squeeze(tmeanregrid(:,:,1))
% uistack(h1,'bottom')
% caxis([0,100])
% colormap(brewermap([],'Blues'))
% colorbar
% %  === Mapping capacity===


%% 4. Calculate anomaly of each variables for 2014, 2015 and 2016

% uniqyear                           = unique(yearnumber);
% 
% for ii = 1:length(modlon)
%     for jj = 1:length(modlat)
%         if ~ismember(LC(jj,ii,10),[1,4,5,6,7,8,9,10])
%             PPTltmean(ii,jj,:)     = NaN;
%             PPTltstd(ii,jj,:)      = NaN;
%             Tmltmean(ii,jj,:)      = NaN;
%             Tmltstd(ii,jj,:)       = NaN;
%             VPDltmean(ii,jj,:)     = NaN;
%             VPDltstd(ii,jj,:)      = NaN;            
%         else
%             % annual mean
%             for kk = 1:length(uniqyear)            
%                 PPTannMean(ii,jj,kk) = nansum(pptregrid(ii,jj,yearnumber==uniqyear(kk)));
%                 TmannMean(ii,jj,kk)  = nanmean(tmeanregrid(ii,jj,yearnumber==uniqyear(kk)));
%                 VPDannMean(ii,jj,kk) = nanmean(vpdregrid(ii,jj,yearnumber==uniqyear(kk)));
%             end
%             % long-term (2000-2010 -- pre-drought) EVI mean and variance
%             PPTltMean(ii,jj)         = nanmean(PPTannMean(ii,jj,uniqyear<=2010 & uniqyear>=2000));
%             PPTltstd(ii,jj)          = nanstd(PPTannMean(ii,jj,uniqyear<=2010 & uniqyear>=2000),0,3);
%             TmltMean(ii,jj)          = nanmean(TmannMean(ii,jj,uniqyear<=2010 & uniqyear>=2000));
%             Tmltstd(ii,jj)           = nanstd(TmannMean(ii,jj,uniqyear<=2010 & uniqyear>=2000),0,3);
%             VPDltmean(ii,jj)         = nanmean(VPDannMean(ii,jj,uniqyear<=2010 & uniqyear>=2000));
%             VPDltstd(ii,jj)          = nanstd(VPDannMean(ii,jj,uniqyear<=2010 & uniqyear>=2000),0,3);            
%         end
%     end 
% end


% % === Mapping capacity===
% calboundary         = shaperead('usastatehi', 'UseGeoCoords', true,...
%                                 'Selector',{@(name) strcmpi(name,'California'), 'Name'});
% usamap('California')
% geoshow(calboundary,'DisplayType','polygon','FaceAlpha',0)
% in                  = inpolygon(modlatmesh,modlonmesh,calboundary.Lat,calboundary.Lon);
% ppttest             = squeeze(PPTannMean(:,:,1));
% ppttest(~in)        = NaN;
% h1                  = pcolorm(modlatmesh,modlonmesh,ppttest); %squeeze(tmeanregrid(:,:,1))
% uistack(h1,'bottom')
% colormap(brewermap([],'Blues'))
% colorbar
% % === Mapping capacity===


%%  5. Annual anomaly of PRISM variables
% for ii = 1:length(modlon)
%     for jj = 1:length(modlat)
%         
%             PPTano2013(ii,jj) = (PPTannMean(ii,jj,14) - PPTltMean(ii,jj))/PPTltstd(ii,jj);
%             PPTano2014(ii,jj) = (PPTannMean(ii,jj,15) - PPTltMean(ii,jj))/PPTltstd(ii,jj);
%             PPTano2015(ii,jj) = (PPTannMean(ii,jj,16) - PPTltMean(ii,jj))/PPTltstd(ii,jj);
%             PPTano2016(ii,jj) = (PPTannMean(ii,jj,17) - PPTltMean(ii,jj))/PPTltstd(ii,jj);
%             
%             Tmano2013(ii,jj)  = (TmannMean(ii,jj,14) - TmltMean(ii,jj))/Tmltstd(ii,jj);
%             Tmano2014(ii,jj)  = (TmannMean(ii,jj,15) - TmltMean(ii,jj))/Tmltstd(ii,jj);
%             Tmano2015(ii,jj)  = (TmannMean(ii,jj,16) - TmltMean(ii,jj))/Tmltstd(ii,jj);
%             Tmano2016(ii,jj)  = (TmannMean(ii,jj,17) - TmltMean(ii,jj))/Tmltstd(ii,jj);
%             
%             VPDano2013(ii,jj) = (VPDannMean(ii,jj,14) - VPDltmean(ii,jj))/VPDltstd(ii,jj);
%             VPDano2014(ii,jj) = (VPDannMean(ii,jj,15) - VPDltmean(ii,jj))/VPDltstd(ii,jj);
%             VPDano2015(ii,jj) = (VPDannMean(ii,jj,16) - VPDltmean(ii,jj))/VPDltstd(ii,jj);
%             VPDano2016(ii,jj) = (VPDannMean(ii,jj,17) - VPDltmean(ii,jj))/VPDltstd(ii,jj);
%             
%             
%     end 
% end


save('/Volumes/XiYangBackUp/Projects/6.CalDrought/PRISMregrid_500m.mat')
























