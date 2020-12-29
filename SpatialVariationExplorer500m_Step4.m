%% SpatialVariationExplorer.m
%  Explore what drives the spatial variations in magnitude of breakpoints,
%  and the slope before and after the breakpoint
%  Final step to generate the data for spatial
%  500m

clear variables
clc

%% 1. Load all the data
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/PRISMregrid_500m.mat','pptregrid','tmeanregrid','vpdregrid','yearnumber','monthnumber',...
    'modlon','modlat','modlatmesh','modlonmesh','PPTseasonmean','PPTseasonstd','Tmseasonmean','Tmseasonstd','VPDseasonmean','VPDseasonstd')
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/CalEVI_2000_2019_500m_BP_all.mat','mag','slope','bp','lonroi','latroi','nobp','bpconf','sigslope')
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/CalEVI_2000_2019_500m_LConly.mat','LC')
modis_file              = '/Volumes/XiYangBackUp/Data/2.SatelliteData/2.MODIS/MOD13A1.006_500m_aid0001.nc';
lat                     = ncread(modis_file,'lat');
lon                     = ncread(modis_file,'lon');
doy                     = ncread(modis_file,'time');
datetime_doy            = datetime(doy+datenum(2000,1,1),'ConvertFrom','datenum');
years                   = year(datetime_doy);

load('/Volumes/XiYangBackUp/Projects/6.CalDrought/wtdregrid_500m.mat','wtdregrid')
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/demregrid_500m.mat','demregrid')
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/MonthlyAnomaly_500m.mat','EVIannMean','EVIltmean',...
    'PPTann_ano','Tmann_ano','VPDann_ano','PPTannMean','TmannMean','VPDannMean')
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/soilthickness_500m.mat','soilthickness');

load('/Volumes/XiYangBackUp/Projects/6.CalDrought/sawc_500m.mat','sawcregrid','soillatmesh','soilonmesh');
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/sawc200cm_500m.mat','sawcregrid200');

% BFAST results
BP1                 = squeeze(bp(:,:,1));
BP2                 = squeeze(bp(:,:,2));
slope1              = squeeze(slope(:,:,1));
slope2              = squeeze(slope(:,:,2));
slope3              = squeeze(slope(:,:,3));
sigslope1           = squeeze(sigslope(:,:,1));
sigslope2           = squeeze(sigslope(:,:,2));
sigslope3           = squeeze(sigslope(:,:,3));
mag1                = mag(:,:,1);
mag2                = mag(:,:,2);
bprange             = bpconf(:,:,3) - bpconf(:,:,1);

% rotate wtdregrid & demregrid
wtdregrid           = wtdregrid';
demregrid           = demregrid';

% fire pixels screening
load('/Volumes/XiYangBackUp/Projects/6.CalDrought/modisfire_500m.mat','modfire')
firepixel           = modfire > 0;
firepixel           = firepixel';

% land cover mask
%years               = year(doy+datenum(2000,1,1));
LC1                 = squeeze(LC(:,:,17));

% make sure these pixels are always within the group 2000, 2011, and 2018
LCmask_deforestation= (LC(:,:,1) == 1) & (LC(:,:,1) ~= LC(:,:,17));
LCmask              = ismember(squeeze(LC(:,:,1)),[1,4,5,6,7,8,9,10]) & ... 
                      ismember(squeeze(LC(:,:,10)),[1,4,5,6,7,8,9,10]) & ...
                      ismember(squeeze(LC(:,:,17)),[1,4,5,6,7,8,9,10]) & ...
                      ~LCmask_deforestation;
                  
Tmseasonmin = nanmin(Tmseasonmean,[],3);
Tmseasonmin1 = Tmseasonmin';

% %load SPEI12
% load('/Volumes/XiYangBackUp/Projects/6.CalDrought/SPEI12.mat','spei12')
% load('/Volumes/XiYangBackUp/Projects/6.CalDrought/SPEI06.mat','spei06')
% spei12year = reshape((repmat((2000:2016)',1,12))',1,204);
% spei12year = spei12year(1:202);
% spei12mon  = reshape((repmat((1:12)',1,17)),1,204);
% spei12mon  = spei12mon(1:202);
% spei12(~isfinite(spei12)) = NaN;
% spei06(~isfinite(spei06)) = NaN;

%% 2. GIS files
calboundary         = shaperead('usastatehi', 'UseGeoCoords', true,...
            'Selector',{@(name) strcmpi(name,'California'), 'Name'});
[XX,YY]             = meshgrid(double(lon),double(lat));
in                  = inpolygon(XX,YY,calboundary.Lon,calboundary.Lat); 

%% 3. Screening

BP1(BP1<=0 | BP1>=500)             = NaN;
BP2(BP2<=0 | BP2>=500)             = NaN;

EVItslength                        = length(doy)-228+1;
rel_mag1                           = mag1./EVIltmean;
slope1_relpercent                  = slope1.*(BP1*16/365)./EVIltmean;
slope2_relpercent                  = slope2.*((EVItslength-BP1)*16/365)./EVIltmean;

% Just make those non-BP pixels: if there is no BP, then the slope1 and
% slope2 are the same as slope1 spanning the entire length of time-series
slope1_relpercent(isnan(BP1))      = slope1(isnan(BP1)).*((EVItslength-228+1)*16/365)./EVIltmean(isnan(BP1));
slope1_relpercent(bprange>23)     = slope1(bprange>23).*((EVItslength-228+1)*16/365)./EVIltmean(bprange>23);
slope2_relpercent(isnan(BP1))      = slope1(isnan(BP1)).*((EVItslength-228+1)*16/365)./EVIltmean(isnan(BP1));
slope2_relpercent(bprange>23)     = slope1(bprange>23).*((EVItslength-228+1)*16/365)./EVIltmean(bprange>23);

LC1                                = uint8(LC1);
nobp(~in | firepixel | ~LCmask)    = NaN;
BP1(~in | firepixel | ~LCmask)     = NaN;
BP2(~in | firepixel | ~LCmask)     = NaN;
slope1(~in | firepixel | ~LCmask)  = NaN;
slope2(~in | firepixel | ~LCmask)  = NaN;
slope3(~in | firepixel | ~LCmask)  = NaN;
sigslope1(~in | firepixel | ~LCmask)  = NaN;
sigslope2(~in | firepixel | ~LCmask)  = NaN;
sigslope3(~in | firepixel | ~LCmask)  = NaN;
mag1(~in | firepixel | ~LCmask)    = NaN;
mag2(~in | firepixel | ~LCmask)    = NaN;
rel_mag1(~in | firepixel | ~LCmask)= NaN;
LC1(~in | firepixel  | ~LCmask)    = NaN;
LC1(LC1==0)                        = NaN;
wtdregrid(~in | firepixel | ~LCmask)= NaN;
demregrid(~in | firepixel | ~LCmask)= NaN;
Tmseasonmin1(~in | firepixel | ~LCmask)= NaN;
bprange(~in | firepixel | ~LCmask)    = NaN;


%% 4. Calculate anomaly

PPTano           = nan(length(modlon),length(modlat));
Tmano            = nan(length(modlon),length(modlat));
VPDano           = nan(length(modlon),length(modlat));
EVIbp            = nan(length(modlon),length(modlat));
EVIbp_ano        = nan(length(modlon),length(modlat));
% cwdbp            = nan(length(modlon),length(modlat));
% cuwdbp           = nan(length(modlon),length(modlat));
% cuwdbp2          = nan(length(modlon),length(modlat));
% cuwdbp3          = nan(length(modlon),length(modlat));
% cuwdbp4          = nan(length(modlon),length(modlat));
% cuwdbp5          = nan(length(modlon),length(modlat));
% cuwdbp6          = nan(length(modlon),length(modlat));

PPTano_2412      = nan(length(modlon),length(modlat));
Tano_2412        = nan(length(modlon),length(modlat));
VPDano_2412      = nan(length(modlon),length(modlat));
PPTano_1206      = nan(length(modlon),length(modlat));
Tano_1206        = nan(length(modlon),length(modlat));
VPDano_1206      = nan(length(modlon),length(modlat));
PPTano_0606      = nan(length(modlon),length(modlat));
Tano_0606        = nan(length(modlon),length(modlat));
VPDano_0606      = nan(length(modlon),length(modlat));
PPTano_0011      = nan(length(modlon),length(modlat));
Tano_0011        = nan(length(modlon),length(modlat));
VPDano_0011      = nan(length(modlon),length(modlat));
PPTano_1201      = nan(length(modlon),length(modlat));
Tano_1201        = nan(length(modlon),length(modlat));
VPDano_1201      = nan(length(modlon),length(modlat));

% SPEI12_BP        = nan(length(modlon),length(modlat));
% SPEI12_WINDOW    = nan(length(modlon),length(modlat));
% SPEI2412_BP      = nan(length(modlon),length(modlat));
% SPEI0012_BP      = nan(length(modlon),length(modlat));
% 
% SPEI06_BP        = nan(length(modlon),length(modlat));
% SPEI06_WINDOW    = nan(length(modlon),length(modlat));
% SPEI06_2412_BP   = nan(length(modlon),length(modlat));
% SPEI06_0012_BP   = nan(length(modlon),length(modlat));

% cwdvalue         = cwddata.cwdregrid;
% cwdyear          = cwddata.yearnumber;
% cwdmonth         = cwddata.monthnumber;
PPTannual        = nan(length(modlon),length(modlat));
Tannual          = nan(length(modlon),length(modlat));
VPDannual        = nan(length(modlon),length(modlat));
PPTano_duringBP  = nan(length(modlon),length(modlat));
Tano_duringBP    = nan(length(modlon),length(modlat));
VPDano_duringBP  = nan(length(modlon),length(modlat));
%cwd_duringdrought= nan(length(modlon),length(modlat));
PPTlt            = nan(length(modlon),length(modlat)); 
Tmlt             = nan(length(modlon),length(modlat)); 
VPDlt            = nan(length(modlon),length(modlat));
%CWDlt            = nan(length(modlon),length(modlat)); 
%cwdannmean       = nan(length(modlon),length(modlat),6);

PPTano_duringBP2  = nan(length(modlon),length(modlat));
Tano_duringBP2    = nan(length(modlon),length(modlat));
VPDano_duringBP2  = nan(length(modlon),length(modlat));
PPTano_duringBP3  = nan(length(modlon),length(modlat));
Tano_duringBP3    = nan(length(modlon),length(modlat));
VPDano_duringBP3  = nan(length(modlon),length(modlat));
PPTano_duringBP4  = nan(length(modlon),length(modlat));
Tano_duringBP4    = nan(length(modlon),length(modlat));
VPDano_duringBP4  = nan(length(modlon),length(modlat));

PPT_bp_3yr_ano    = nan(length(modlon),length(modlat));
Tmp_bp_3yr_ano    = nan(length(modlon),length(modlat));
VPD_bp_3yr_ano    = nan(length(modlon),length(modlat));

PPT_bp_3yr_dif    = nan(length(modlon),length(modlat));
Tmp_bp_3yr_dif    = nan(length(modlon),length(modlat));
VPD_bp_3yr_dif    = nan(length(modlon),length(modlat));

PPT_post_trend    = nan(length(modlon),length(modlat));
Tmp_post_trend    = nan(length(modlon),length(modlat));
VPD_post_trend    = nan(length(modlon),length(modlat));

PPT_post_trend_1  = nan(length(modlon),length(modlat));
Tmp_post_trend_1  = nan(length(modlon),length(modlat));
VPD_post_trend_1  = nan(length(modlon),length(modlat));

PPT_pre_ano       = nan(length(modlon),length(modlat));
Tmp_pre_ano       = nan(length(modlon),length(modlat));
VPD_pre_ano       = nan(length(modlon),length(modlat));

timelags = [1,3,6,12,24];
%timelags = [-12]; %-1,-3,-6,-12,
pcor_all = nan(7,25);
for i_lags = 4:4
    
    timelag = timelags(i_lags); %month
    for ii = 1:length(modlon)
        for jj = 1:length(modlat)

            if isnan(BP1(jj,ii))

                 PPTano(ii,jj)              = NaN;
                 Tmano(ii,jj)               = NaN;
                 VPDano(ii,jj)              = NaN;
                 EVIbp(ii,jj)               = NaN;
                 EVIbp_ano(ii,jj)           = NaN;
                 PPTano_duringBP(ii,jj)     = NaN;
                 Tano_duringBP(ii,jj)       = NaN;
                 VPDano_duringBP(ii,jj)     = NaN;
                 PPTlt(ii,jj)               = nanmean(PPTannMean(ii,jj,:));
                 Tmlt(ii,jj)                = nanmean(TmannMean(ii,jj,:));
                 VPDlt(ii,jj)               = nanmean(VPDannMean(ii,jj,:));
      
                 
                 
            else
                
                                
                 PPTlt(ii,jj)               = nanmean(PPTannMean(ii,jj,:));
                 Tmlt(ii,jj)                = nanmean(TmannMean(ii,jj,:));
                 VPDlt(ii,jj)               = nanmean(VPDannMean(ii,jj,:));
                
                 
                bp_tmp = doy(BP1(jj,ii)+228) + datenum(2000,1,1);
                bp_tmp = datetime(bp_tmp,'ConvertFrom','datenum');
                bp_year= year(bp_tmp);
                bp_mon = month(bp_tmp);

                time_index    = find(yearnumber==bp_year & monthnumber == bp_mon);
                if timelag > 0
                    bp_mon1 = unique(mod(((12+bp_mon-timelag):(12+bp_mon-1))-1,12))+1;
                    indx    = (time_index-timelag):(time_index-1);
%                    indx2   = (cwd_t_index-timelag):(cwd_t_index-1);
%                    indx2   = indx2(indx2>0);
%                    indx22  = (cwd_t_index-timelag):(cwd_t_index);
%                    indx22  = indx22(indx22>0);
                else
                    bp_mon1 = unique(mod(((12+bp_mon):(12+bp_mon-timelag))-1,12))+1;
                    indx    = time_index:(time_index-timelag);
%                    indx2   = cwd_t_index:(cwd_t_index-timelag);
                end

                for kk = 1:11
                    PPTtmp(kk)= nanmean(pptregrid(ii,jj,yearnumber==1999+kk & ismember(monthnumber,bp_mon1)));
                    Tmtmp(kk) = nanmean(tmeanregrid(ii,jj,yearnumber==1999+kk & ismember(monthnumber,bp_mon1)));
                    VPDtmp(kk) = nanmean(vpdregrid(ii,jj,yearnumber==1999+kk & ismember(monthnumber,bp_mon1)));
                end

                PPTano(ii,jj) = (nanmean(pptregrid(ii,jj,indx)) - nanmean(PPTseasonmean(ii,jj,bp_mon1)))/nanstd(PPTtmp,0,2);
                Tmano(ii,jj)  = (nanmean(tmeanregrid(ii,jj,indx)) - nanmean(Tmseasonmean(ii,jj,bp_mon1)))/nanstd(Tmtmp,0,2);
                VPDano(ii,jj) = (nanmean(vpdregrid(ii,jj,indx)) - nanmean(VPDseasonmean(ii,jj,bp_mon1)))/nanstd(VPDtmp,0,2);
                EVIbp(ii,jj)  = EVIannMean(ii,jj,bp_year-2000+1-1);
                EVIbp_ano(ii,jj) = (EVIannMean(ii,jj,bp_year-2000+1-1) - nanmean(EVIannMean(ii,jj,1:11)))/nanstd(EVIannMean(ii,jj,1:11),0,3);
                
                PPT_bp_3yr_ano(ii,jj) = mean(PPTann_ano(ii,jj,(bp_year-1999):(bp_year-1999+2)));
                Tmp_bp_3yr_ano(ii,jj) = mean(Tmann_ano(ii,jj,(bp_year-1999):(bp_year-1999+2)));
                VPD_bp_3yr_ano(ii,jj) = mean(VPDann_ano(ii,jj,(bp_year-1999):(bp_year-1999+2)));

                PPT_bp_3yr_dif(ii,jj) = PPTann_ano(ii,jj,(bp_year-1999+2)) - PPTann_ano(ii,jj,(bp_year-1999));
                Tmp_bp_3yr_dif(ii,jj) = Tmann_ano(ii,jj,(bp_year-1999+2)) - Tmann_ano(ii,jj,(bp_year-1999));
                VPD_bp_3yr_dif(ii,jj) = VPDann_ano(ii,jj,(bp_year-1999+2)) - VPDann_ano(ii,jj,(bp_year-1999));
                
                PPT_post_trend(ii,jj) = (PPTann_ano(ii,jj,20) - PPTann_ano(ii,jj,(bp_year-1999)))/(20-bp_year+1999);
                Tmp_post_trend(ii,jj) = (Tmann_ano(ii,jj,20) - Tmann_ano(ii,jj,(bp_year-1999)))/(20-bp_year+1999);
                VPD_post_trend(ii,jj) = (VPDann_ano(ii,jj,20) - VPDann_ano(ii,jj,(bp_year-1999)))/(20-bp_year+1999);

                PPT_post_trend_1(ii,jj) = (PPTannMean(ii,jj,20) - PPTannMean(ii,jj,(bp_year-1999)))/(20-bp_year+1999);
                Tmp_post_trend_1(ii,jj) = (TmannMean(ii,jj,20) - TmannMean(ii,jj,(bp_year-1999)))/(20-bp_year+1999);
                VPD_post_trend_1(ii,jj) = (VPDannMean(ii,jj,20) - VPDannMean(ii,jj,(bp_year-1999)))/(20-bp_year+1999);                
                
                
                PPT_pre_ano(ii,jj)    = mean(PPTann_ano(ii,jj,(bp_year-1999-1):(bp_year-1999)));
                Tmp_pre_ano(ii,jj)    = mean(Tmann_ano(ii,jj,(bp_year-1999-1):(bp_year-1999)));
                VPD_pre_ano(ii,jj)    = mean(VPDann_ano(ii,jj,(bp_year-1999-1):(bp_year-1999)));
%                cwdbp(ii,jj)  = nanmean(cwdvalue(ii,jj,indx2));
                
 %              WDi           = nan(length(indx22)-1,1);
%               WDi(1,1)      = nanmin(WD(ii,jj,indx22(1)),0);
%                for cwd_ii = 2:((length(indx22)))
%                    WDi(cwd_ii,1) = nanmin(WD(ii,jj,indx22(cwd_ii)) + WDi(cwd_ii-1,1),0);
%                end
%                cuwdbp(ii,jj) = WDi(length(indx22)-3,1);
%                cuwdbp2(ii,jj) = WDi(length(indx22)-6,1);
%                cuwdbp3(ii,jj) = WDi(length(indx22)-9,1);
%                cuwdbp4(ii,jj) = WDi(length(indx22),1);
%                if sum((length(indx22)-12)<0) >0 
%                    cuwdbp5(ii,jj) = NaN;
%                else
%                    cuwdbp5(ii,jj) = WDi(length(indx22)-12,1);
%                end
%                
%                if sum((length(indx22)-23)<0) >0 
%                    cuwdbp6(ii,jj) = NaN;
%                else
%                    cuwdbp6(ii,jj) = WDi(1,1);
%                end
%                
%                cwd_duringdrought(ii,jj) = nanmean(cwdvalue(ii,jj,cwdyear==bp_year));
                
                uniqyears   = unique(yearnumber);
                
                for year_i = 1:length(uniqyears)
                    
                    PPTannual(ii,jj,year_i) = nanmean(pptregrid(ii,jj,yearnumber==uniqyears(year_i)));
                    Tannual(ii,jj,year_i)   = nanmean(tmeanregrid(ii,jj,yearnumber==uniqyears(year_i)));
                    VPDannual(ii,jj,year_i) = nanmean(vpdregrid(ii,jj,yearnumber==uniqyears(year_i)));
                    
                end

                bp_mon2                = unique(mod(((12+bp_mon):(12+bp_mon+2))-1,12))+1;
                indx3                  = (time_index):(time_index+2);
                for kkk = 1:11
                    PPTtmp2(kkk)= nanmean(pptregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon2)));
                    Tmtmp2(kkk) = nanmean(tmeanregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon2)));
                    VPDtmp2(kkk) = nanmean(vpdregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon2)));
                end

                PPTano_duringBP(ii,jj) = (nanmean(pptregrid(ii,jj,indx3)) - nanmean(PPTseasonmean(ii,jj,bp_mon2)))/nanstd(PPTtmp2,0,2);
                Tano_duringBP(ii,jj)   = (nanmean(tmeanregrid(ii,jj,indx3)) - nanmean(Tmseasonmean(ii,jj,bp_mon2)))/nanstd(Tmtmp2,0,2);
                VPDano_duringBP(ii,jj) = (nanmean(vpdregrid(ii,jj,indx3)) - nanmean(VPDseasonmean(ii,jj,bp_mon2)))/nanstd(VPDtmp2,0,2);

                bp_mon3                = unique(mod(((12+bp_mon-3):(12+bp_mon+3))-1,12))+1;
                indx4                  = (time_index-3):(time_index+3);
                for kkk = 1:11
                    PPTtmp3(kkk)= nanmean(pptregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon3)));
                    Tmtmp3(kkk) = nanmean(tmeanregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon3)));
                    VPDtmp3(kkk) = nanmean(vpdregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon3)));
                end

                PPTano_duringBP2(ii,jj) = (nanmean(pptregrid(ii,jj,indx4)) - nanmean(PPTseasonmean(ii,jj,bp_mon3)))/nanstd(PPTtmp3,0,2);
                Tano_duringBP2(ii,jj)   = (nanmean(tmeanregrid(ii,jj,indx4)) - nanmean(Tmseasonmean(ii,jj,bp_mon3)))/nanstd(Tmtmp3,0,2);
                VPDano_duringBP2(ii,jj) = (nanmean(vpdregrid(ii,jj,indx4)) - nanmean(VPDseasonmean(ii,jj,bp_mon3)))/nanstd(VPDtmp3,0,2); 
                
                bp_mon4                = unique(mod(((12+bp_mon-5):(12+bp_mon+6))-1,12))+1;
                indx5                  = (time_index-5):(time_index+6);
                for kkk = 1:11
                    PPTtmp4(kkk)= nanmean(pptregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon4)));
                    Tmtmp4(kkk) = nanmean(tmeanregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon4)));
                    VPDtmp4(kkk) = nanmean(vpdregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon4)));
                end

                PPTano_duringBP3(ii,jj) = (nanmean(pptregrid(ii,jj,indx5)) - nanmean(PPTseasonmean(ii,jj,bp_mon4)))/nanstd(PPTtmp4,0,2);
                Tano_duringBP3(ii,jj)   = (nanmean(tmeanregrid(ii,jj,indx5)) - nanmean(Tmseasonmean(ii,jj,bp_mon4)))/nanstd(Tmtmp4,0,2);
                VPDano_duringBP3(ii,jj) = (nanmean(vpdregrid(ii,jj,indx5)) - nanmean(VPDseasonmean(ii,jj,bp_mon4)))/nanstd(VPDtmp4,0,2);                      
                
                PPTano_duringBP4(ii,jj) = (PPTannual(ii,jj,uniqyears==bp_year) - nanmean(PPTannual(ii,jj,:)))/nanstd(PPTannual(ii,jj,:),0,3);
                Tano_duringBP4(ii,jj)   = (Tannual(ii,jj,uniqyears==bp_year) - nanmean(Tannual(ii,jj,:)))/nanstd(Tannual(ii,jj,:),0,3);
                VPDano_duringBP4(ii,jj) = (VPDannual(ii,jj,uniqyears==bp_year) - nanmean(VPDannual(ii,jj,:)))/nanstd(VPDannual(ii,jj,:),0,3);
                
                
                % AVE_CLIM_ANO
                
                bp_mon5                = unique(mod(((12+bp_mon-24):(12+bp_mon-12))-1,12))+1;
                indx6                  = (time_index-24):(time_index-12);
                for kkk = 1:11
                    PPTtmp5(kkk)= nanmean(pptregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon5)));
                    Tmtmp5(kkk) = nanmean(tmeanregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon5)));
                    VPDtmp5(kkk) = nanmean(vpdregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon5)));
                end
                
                PPTano_2412(ii,jj)     = (nanmean(pptregrid(ii,jj,indx6)) - nanmean(PPTseasonmean(ii,jj,bp_mon5)))/nanstd(PPTtmp5,0,2);
                Tano_2412(ii,jj)       = (nanmean(tmeanregrid(ii,jj,indx6)) - nanmean(Tmseasonmean(ii,jj,bp_mon5)))/nanstd(Tmtmp5,0,2);
                VPDano_2412(ii,jj)     = (nanmean(vpdregrid(ii,jj,indx6)) - nanmean(VPDseasonmean(ii,jj,bp_mon5)))/nanstd(VPDtmp5,0,2);
                
                bp_mon6                = unique(mod(((12+bp_mon-12):(12+bp_mon-6))-1,12))+1;
                indx7                  = (time_index-12):(time_index-6);
                for kkk = 1:11
                    PPTtmp6(kkk)= nanmean(pptregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon6)));
                    Tmtmp6(kkk) = nanmean(tmeanregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon6)));
                    VPDtmp6(kkk) = nanmean(vpdregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon6)));
                end
                
                PPTano_1206(ii,jj)     = (nanmean(pptregrid(ii,jj,indx7)) - nanmean(PPTseasonmean(ii,jj,bp_mon6)))/nanstd(PPTtmp6,0,2);
                Tano_1206(ii,jj)       = (nanmean(tmeanregrid(ii,jj,indx7)) - nanmean(Tmseasonmean(ii,jj,bp_mon6)))/nanstd(Tmtmp6,0,2);
                VPDano_1206(ii,jj)     = (nanmean(vpdregrid(ii,jj,indx7)) - nanmean(VPDseasonmean(ii,jj,bp_mon6)))/nanstd(VPDtmp6,0,2);                
                
                bp_mon7                = unique(mod(((12+bp_mon-6):(12+bp_mon+5))-1,12))+1;
                indx8                  = (time_index-6):(time_index+5);
                for kkk = 1:11
                    PPTtmp7(kkk)= nanmean(pptregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon7)));
                    Tmtmp7(kkk) = nanmean(tmeanregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon7)));
                    VPDtmp7(kkk) = nanmean(vpdregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon7)));
                end
                
                PPTano_0606(ii,jj)     = (nanmean(pptregrid(ii,jj,indx8)) - nanmean(PPTseasonmean(ii,jj,bp_mon7)))/nanstd(PPTtmp7,0,2);
                Tano_0606(ii,jj)       = (nanmean(tmeanregrid(ii,jj,indx8)) - nanmean(Tmseasonmean(ii,jj,bp_mon7)))/nanstd(Tmtmp7,0,2);
                VPDano_0606(ii,jj)     = (nanmean(vpdregrid(ii,jj,indx8)) - nanmean(VPDseasonmean(ii,jj,bp_mon7)))/nanstd(VPDtmp7,0,2);   
                
                bp_mon8                = unique(mod(((12+bp_mon):(12+bp_mon+11))-1,12))+1;
                indx9                  = (time_index):(time_index+11);
                for kkk = 1:11
                    PPTtmp8(kkk)= nanmean(pptregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon8)));
                    Tmtmp8(kkk) = nanmean(tmeanregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon8)));
                    VPDtmp8(kkk) = nanmean(vpdregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon8)));
                end
                
                PPTano_0011(ii,jj)     = (nanmean(pptregrid(ii,jj,indx9)) - nanmean(PPTseasonmean(ii,jj,bp_mon8)))/nanstd(PPTtmp8,0,2);
                Tano_0011(ii,jj)       = (nanmean(tmeanregrid(ii,jj,indx9)) - nanmean(Tmseasonmean(ii,jj,bp_mon8)))/nanstd(Tmtmp8,0,2);
                VPDano_0011(ii,jj)     = (nanmean(vpdregrid(ii,jj,indx9)) - nanmean(VPDseasonmean(ii,jj,bp_mon8)))/nanstd(VPDtmp8,0,2);  
                
                bp_mon9                = unique(mod(((12+bp_mon-12):(12+bp_mon-1))-1,12))+1;
                indx10                  = (time_index-12):(time_index-1);
                for kkk = 1:11
                    PPTtmp9(kkk)= nanmean(pptregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon9)));
                    Tmtmp9(kkk) = nanmean(tmeanregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon9)));
                    VPDtmp9(kkk) = nanmean(vpdregrid(ii,jj,yearnumber==1999+kkk & ismember(monthnumber,bp_mon9)));
                end
                
                PPTano_1201(ii,jj)     = (nanmean(pptregrid(ii,jj,indx10)) - nanmean(PPTseasonmean(ii,jj,bp_mon9)))/nanstd(PPTtmp9,0,2);
                Tano_1201(ii,jj)       = (nanmean(tmeanregrid(ii,jj,indx10)) - nanmean(Tmseasonmean(ii,jj,bp_mon9)))/nanstd(Tmtmp9,0,2);
                VPDano_1201(ii,jj)     = (nanmean(vpdregrid(ii,jj,indx10)) - nanmean(VPDseasonmean(ii,jj,bp_mon9)))/nanstd(VPDtmp9,0,2); 

                
   
%                 % SPEI12 at the time of BP
%                 %spei12_ind             = knnsearch(time_spei,bp_tmp,'K',1);
%                 %SPEI12_BP(ii,jj)       = spei12regrid(ii,jj,spei12_ind);
%                 spei12_ind             = (spei12year == bp_year) & (spei12mon == bp_mon);
%                 SPEI12_BP(ii,jj)       = spei12(ii,jj,spei12_ind);
%                 SPEI12_WINDOW(ii,jj)   = nanmean(spei12(ii,jj,(find(spei12_ind)-11):(find(spei12_ind))));
%                 SPEI2412_BP(ii,jj)     = nanmean(spei12(ii,jj,(find(spei12_ind)-23):(find(spei12_ind)-11)));
%                 if find(spei12_ind)+12 > 202
%                     SPEI0012_BP(ii,jj)     = NaN;
%                     SPEI06_0012_BP(ii,jj)  = NaN;
%                 else
%                     SPEI0012_BP(ii,jj)     = nanmean(spei12(ii,jj,(find(spei12_ind)):(find(spei12_ind)+11)));
%                     SPEI06_0012_BP(ii,jj)  = nanmean(spei06(ii,jj,(find(spei12_ind)):(find(spei12_ind)+11)));
%                 end
%                 
%                 SPEI06_BP(ii,jj)       = spei06(ii,jj,spei12_ind);
%                 SPEI06_WINDOW(ii,jj)   = nanmean(spei06(ii,jj,(find(spei12_ind)-11):(find(spei12_ind))));
%                 SPEI06_2412_BP(ii,jj)  = nanmean(spei06(ii,jj,(find(spei12_ind)-23):(find(spei12_ind)-11)));                
                
            end
        end 
    end

%     PPTano(PPTano<-5 | PPTano >5)   = NaN;
%     Tmano(Tmano<-5 | Tmano >5)      = NaN;
%     VPDano(VPDano<-5 | VPDano >5)   = NaN;
% 
%     usamap('California')
%     h2                  = geoshow(calboundary,'DisplayType','polygon','FaceAlpha',0);   
%     h1                  = pcolorm(latroi,lonroi,rel_mag1);
%     caxis([-1.5, 1.5])
%     uistack(h1,'bottom')
%     colormap(brewermap([],'RdBu'))
%     %colormap(flipud(brewermap([],'RdBu')))
%     colorbar
%         
%     print('/Volumes/XiYangBackUp/Projects/6.CalDrought/SPEI12_BP_newcalculation_BP.png',gcf,'-dpng')
%     close(gcf)     

    
    %% 5. Partial regression between BFAST results and PRISM variables
    slope11    = slope1';
    relmag1    = rel_mag1';
    slope22    = slope2';
    nobp1      = nobp';
    wtdregrid  = wtdregrid';
    demregrid  = demregrid';
    BP11       = BP1';
    BP22       = BP2';
    Tmseasonmin1 = Tmseasonmin1';
    bprange1   = bprange';
    
    PPTano2016 = squeeze(PPTann_ano(:,:,17));
    Tmano2016 = squeeze(Tmann_ano(:,:,17));
    VPDano2016 = squeeze(VPDann_ano(:,:,17));

    PPTano2017 = squeeze(PPTann_ano(:,:,18));
    Tmano2017 = squeeze(Tmann_ano(:,:,18));
    VPDano2017 = squeeze(VPDann_ano(:,:,18));
    
    PPTano2018 = squeeze(PPTann_ano(:,:,19));
    Tmano2018 = squeeze(Tmann_ano(:,:,19));
    VPDano2018 = squeeze(VPDann_ano(:,:,19));
    
    PPTano1618 = squeeze(nanmean(PPTann_ano(:,:,17:19),3));
    Tmano1618 = squeeze(nanmean(Tmann_ano(:,:,17:19),3));
    VPDano1618 = squeeze(nanmean(VPDann_ano(:,:,17:19),3));
    
    nobp1(bprange1>23)                 = 0;
    total_relpercent                   = nan(size(slope1_relpercent));
    total_relpercent2016               = nan(size(slope1_relpercent));
    total_relpercent(nobp1'==0)        = slope1_relpercent(nobp1'==0);
    total_relpercent(nobp1'==1)        = slope1_relpercent(nobp1'==1) + rel_mag1(nobp1'==1) + slope2_relpercent(nobp1'==1);
    total_relpercent(nobp1'==2)        = NaN;

    total_relpercent2016(nobp1'==0)    = slope1_relpercent(nobp1'==0);
    total_relpercent2016(nobp1'==1)    = slope1_relpercent(nobp1'==1) + rel_mag1(nobp1'==1) + slope2_relpercent(nobp1'==1) .* (389-228-BP1(nobp1'==1))./(length(doy)-228-BP1(nobp1'==1)) * (16/365); %First data point in 2017 is 389
    total_relpercent2016(nobp1'==2)    = NaN;    
    
    
    slope1_relpercent(~in | firepixel | ~LCmask) = NaN;
    slope2_relpercent(~in | firepixel | ~LCmask) = NaN;
    total_relpercent(~in | firepixel | ~LCmask) = NaN;
    total_relpercent2016(~in | firepixel | ~LCmask) = NaN;
    
    slope1_relpercent1 = slope1_relpercent';
    slope2_relpercent1 = slope2_relpercent';
    total_relpercent1 = total_relpercent';
    total_relpercent2016_1 = total_relpercent2016';
 
       usamap('California')
    h2                  = geoshow(calboundary,'DisplayType','polygon','FaceAlpha',0);   
    h1                  = pcolorm(lat,lon,EVIltmean);
    caxis([0, 1])
    uistack(h1,'bottom')
    colormap(brewermap([],'RdBu'))
    %colormap(flipud(brewermap([],'RdBu')))
    colorbar 
    
    
%     %wtdregrid(~(nobp==1))=NaN;
%     
%     lc_type    = [[1,1];[6,7];[8,9];[10,10]];
%     
%     for lc_i   = 1:5
%         if lc_i == 1
%             [pcor,sig] = partialcorr([relmag1(:),PPTano(:),Tmano(:),VPDano(:),cwdbp(:),wtdregrid(:),EVIbp(:),EVIltmean(:)],'Rows','complete');
%         else
%             subsample  = LC1'==lc_type(lc_i-1,1)| LC1'==lc_type(lc_i-1,2);
%             [pcor,sig] = partialcorr([relmag1(subsample),PPTano(subsample),Tmano(subsample),VPDano(subsample),cwdbp(subsample),wtdregrid(subsample),EVIbp(subsample),EVIltmean(subsample)],'Rows','complete');
%         end
%         pcor_all(:,(lc_i-1)*5+i_lags) = pcor(2:8,1);
%     end
    %[pcor,sig] = partialcorr([slope11(:),PPTano(:),Tmano(:),VPDano(:),wtdregrid(:)],'Rows','complete')
    
    
    %% 6. Stepwise regression
    LC11        = LC1';
    nonnanidx   = ~isnan(PPTano(:)) & ~isnan(Tmano(:)) & ~isnan(VPDano(:)) & ~isnan(wtdregrid(:)) & ~isnan(LC11(:)) & ~isnan(BP11(:)) & ~isnan(EVIltmean(:)) & ~isnan(EVIbp(:)) & ~isnan(EVIbp_ano(:))...
                    & (slope1_relpercent1(:)>-5 & slope1_relpercent1(:)<5) ...
                    & (slope2_relpercent1(:)>-5 & slope2_relpercent1(:)<5) ...
                    & (relmag1(:)>-5 & relmag1(:)<5) ...
                    & (total_relpercent1(:)>-5 & total_relpercent1(:)<5); 
            
    EVIdeviation = PPTano2016; % just put a placeholder here for EVIdeviation.
    % SPEI12_BP(:),SPEI12_WINDOW(:),SPEI2412_BP(:),SPEI0012_BP(:),SPEI06_BP(:),SPEI06_WINDOW(:),SPEI06_2412_BP(:),SPEI06_0012_BP(:),,cuwdbp(:),cwd_duringdrought(:)
    XX          = [PPTano(:),Tmano(:),VPDano(:),PPTano_duringBP(:),Tano_duringBP(:),VPDano_duringBP(:),PPTano_duringBP3(:),Tano_duringBP3(:),VPDano_duringBP3(:),PPTano_duringBP4(:),Tano_duringBP4(:),VPDano_duringBP4(:),...
                    EVIbp_ano(:),double(LC11(:)),wtdregrid(:),EVIbp(:),EVIltmean(:),sawcregrid(:),sawcregrid200(:),demregrid(:),...
                    PPTlt(:),Tmlt(:),VPDlt(:),soilthickness(:), Tmseasonmin1(:),PPTano_duringBP2(:),Tano_duringBP2(:),VPDano_duringBP2(:),...
                    PPTano_2412(:),Tano_2412(:),VPDano_2412(:),PPTano_1206(:),Tano_1206(:),VPDano_1206(:),PPTano_0606(:),Tano_0606(:),VPDano_0606(:),PPTano_0011(:),Tano_0011(:),VPDano_0011(:),PPTano_1201(:),Tano_1201(:),VPDano_1201(:),...
                    EVIdeviation(:),PPTano2016(:),Tmano2016(:),VPDano2016(:),PPTano2017(:),Tmano2017(:),VPDano2017(:),PPTano2018(:),Tmano2018(:),...
                    VPDano2018(:),PPTano1618(:),Tmano1618(:),VPDano1618(:),PPT_bp_3yr_ano(:),Tmp_bp_3yr_ano(:),VPD_bp_3yr_ano(:), ...
                    PPT_bp_3yr_dif(:),Tmp_bp_3yr_dif(:),VPD_bp_3yr_dif(:),PPT_post_trend(:),Tmp_post_trend(:),VPD_post_trend(:), ...
                    PPT_pre_ano(:),Tmp_pre_ano(:),VPD_pre_ano(:),PPT_post_trend_1(:),Tmp_post_trend_1(:),VPD_post_trend_1(:)]; 
    YY          = relmag1(:);
    %All2 = [XX(nonnanidx,:),double(LC11(nonnanidx)),YY(nonnanidx)];

    All         = [XX(nonnanidx,:),YY(nonnanidx),slope11(nonnanidx),slope1_relpercent1(nonnanidx),slope2_relpercent1(nonnanidx), ...
                    double(doy(BP11(nonnanidx)+228) + datenum(2000,1,1)),double(doy(BP22(nonnanidx)+228) + datenum(2000,1,1)),bprange1(nonnanidx),nobp1(nonnanidx),total_relpercent1(nonnanidx), ...
                    sigslope1(nonnanidx),sigslope2(nonnanidx),sigslope3(nonnanidx)];
    All         = double(All);
    All1        = [XX(:,:),YY(:),slope11(:),slope1_relpercent1(:),slope2_relpercent1(:),bprange1(:),nobp1(:),total_relpercent1(:)];
    save('/Volumes/XiYangBackUp/Projects/6.CalDrought/spatialvariation_12192020.mat','-v7.3')
    save('/Volumes/XiYangBackUp/Projects/6.CalDrought/spatialvariation_12192020_1.mat','All')
    
    
    
end

















