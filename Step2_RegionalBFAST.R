## Processing MODIS EVI data
## HDF read from Matlab

rm(list=ls())

stime <- system.time({
  
library('R.matlab')
library('bfast')
library('parallel')
library('MASS')
library('iterators')
library('tictoc')
library('foreach')
library('doParallel')
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("rhdf5")
library('rhdf5')
  
numCores   <- detectCores()
registerDoParallel(numCores/2)
print(numCores)

filename         <- '/Volumes/XiYangBackUp/Projects/6.CalDrought/subsets/CalEVI_2000_2019_500m_Subset_8_8.mat' 
EVI              <- h5read(filename,'EVI')
doy              <- h5read(filename,'doy')
LC               <- h5read(filename,'LC')
EVI              <- EVI[,,228:456] #457 for CMG
doy              <- doy[,,228:456]
dims             <- dim(EVI)

#QA
EVI[EVI==0.0]    <-NA


iii        <- 1
allresult  <- list()
allresult <- foreach(ii=1:dims[1], .packages = c("bfast","foreach")) %dopar%  {
  allresult1 <- foreach(jj=1:dims[2], .packages = c("bfast","foreach")) %dopar%  {

    if(sum(is.na(EVI[ii,jj,]))>=110)
    {
      result     <- NA
    } else if(LC[ii,jj,10] %in% c(1,4,5,6,7,8,9,10))
      {
        time1      <- as.Date(doy[ii,jj,],origin="2000-01-01")
        years      <- unique(substring(time1,1,4))
        DF         <- data.frame(x=time1,y=EVI[ii,jj,])
        DF1        <- na.omit(DF)
        ts_example <- ts(data = DF1$y,frequency=23)
        result     <- bfast(ts_example,h=60/length(ts_example),season="harmonic",max.iter=2,breaks=1,type="OLS-MOSUM",hpc="foreach")
      } else
      {
        result     <- NA
      }
  }
}

})
print(stime)
save.image('/Volumes/XiYangBackUp/Projects/6.CalDrought/Caldorught_8_8.RData') #change the output filename
