%% Read MODIS data
%  This is the Step 0 for JGR paper
%  It reads EVI and land cover data in tiles to be processed in the BFAST
%  analysis (R code). The reason for tiles is that Matlab cannot handle the
%  entire matrix (~2600 x ~2600 x 439 double)

clear variables
clc

% MOD13A1 is VI
modis_file              = '/Volumes/XiYangBackUp/Data/2.SatelliteData/2.MODIS/MOD13A1.006_500m_aid0001.nc';
% MCD12Q1 is land cover
modis_file_2            = '/Volumes/XiYangBackUp/Data/2.SatelliteData/2.MODIS/MCD12Q1.006_500m_aid0001.nc';

lat                     = ncread(modis_file,'lat');
lon                     = ncread(modis_file,'lon');

% Here we break down the files into tiles
for ii = 1:ceil(length(lon)/100)
    for jj = 1:ceil(length(lat)/100)
            
            start_loc               = [(ii-1)*100+1,(jj-1)*100+1,1];         
        if ii == ceil(length(lon)/100) && jj == ceil(length(lat)/100)
            count                   = [length(lon)-start_loc(1),length(lat)-start_loc(2),Inf];
        elseif ii ~= ceil(length(lon)/100) && jj == ceil(length(lat)/100)
            count                   = [100,length(lat)-start_loc(2),Inf];
        elseif ii == ceil(length(lon)/100) && jj ~= ceil(length(lat)/100)
            count                   = [length(lon)-start_loc(1),100,Inf];
        else
            count                   = [100,100,Inf];        
        end
        % doy since 2000-01-01
        doy                     = ncread(modis_file,'_500m_16_days_composite_day_of_the_year',start_loc,count);
        EVI_tmp                 = ncread(modis_file,'_500m_16_days_EVI',start_loc,count);
        QA                      = ncread(modis_file,'_500m_16_days_pixel_reliability',start_loc,count);
        LC                      = ncread(modis_file_2,'LC_Type1',start_loc,count);%,start_loc,count);

        % switch lat and lon dimension for consistency with the 5.6 km data
        % so it is switched to lat, lon, time
        EVI                     = permute(EVI_tmp,[2,1,3]);
        LC                      = permute(LC,[2,1,3]);
        doy                     = permute(doy,[2,1,3]);
        
        % Note that for other use the file path on your computer
        save(['/Volumes/XiYangBackUp/Projects/6.CalDrought/subsets/CalEVI_2000_2019_500m_Subset_' num2str(ii) '_' num2str(jj) '.mat'], '-v7.3'); 
        
    end
end