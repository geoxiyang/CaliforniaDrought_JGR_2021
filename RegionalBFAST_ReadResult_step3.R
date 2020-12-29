## Read the BFAST result from RegionalBFAST_CMG.R 

rm(list=ls())

library('R.matlab')
library('bfast')

stime <- system.time({

  #change file name here. First -- the original EVI file; Second -- bfast result; Third -- result output
inputEVI_name   <- '/Volumes/XiYangBackUp/Projects/6.CalDrought/subsets/CalEVI_2000_2019_500m_Subset_8_8.mat'
bfast_name      <- '/Volumes/XiYangBackUp/Projects/6.CalDrought/Caldorught_8_8.RData'
output_filename <- '/Volumes/XiYangBackUp/Projects/6.CalDrought/Caldrought_output/CalEVI_2000_2019_500m_Subset_8_8_result.mat'
#  1. Read the lat/lon information from matlab file
EVI              <- h5read(inputEVI_name,'EVI')
dims             <- dim(EVI)

#  2. Read the BFAST result
load(bfast_name) #change the filename

#  3. Read BFAST result and put them in matrices
bp         <- array(NA, c(dims[1], dims[2], 2))
slope      <- array(NA, c(dims[1], dims[2], 3))
slopelm    <- array(NA, c(dims[1], dims[2], 3))
sigslope   <- array(NA, c(dims[1], dims[2], 3))
mag        <- array(NA, c(dims[1], dims[2], 2))
nobp       <- array(NA, c(dims[1], dims[2]))
slope1     <- array(NA, 3)
bpconf     <- array(NA, c(dims[1], dims[2], 6))
  
for (ii in c(1:dims[1])) { #
  for (jj in c(1:dims[2]))  { #1:dims[2]
    result       <- allresult[[ii]][[jj]]

    if (is.na(result)) {
      bp[ii,jj,]<-NA
      slope[ii,jj,]<- NA
      mag[ii,jj,]<-NA
      nobp[ii,jj]<-NA
      bpconf[ii,jj,] <- NA
    }
    else if (unname(result$nobp$Vt))
    {
      bp[ii,jj,]<-NA
      bpconf[ii,jj,] <- NA
      xx        <- c(1:length(result$Yt))
      yy        <- result$output[[length(result$output)]]$Tt
      templm    <- lm(yy~xx)
      slopelm[ii,jj,1] <- (365/16)*summary(templm)$coefficients[2,1]
      sigslope[ii,jj,1]<- summary(templm)$coefficients[2,4]
      slope[ii,jj,1]<- (365/16)*unname(result$output[[length(result$output)]]$Tt[length(result$Yt)] - result$output[[length(result$output)]]$Tt[1])/length(result$Yt)
      mag[ii,jj,]<-NA
      nobp[ii,jj]<-0
    } 
    else
    {
      bp[ii,jj,] <- unname(result$output[[length(result$output)]]$Vt.bp)
      bpconf[ii,jj,] <- unname(result$output[[length(result$output)]]$ci.Vt$confint)
      mag[ii,jj,] <- unname(result$Magnitude)
      if (bp[ii,jj,1] == bp[ii,jj,2])
      {
        xx1       <- c(1:bp[ii,jj,1])
        xx2       <- c((bp[ii,jj,1]+1):length(result$Yt))
        yy1       <- result$output[[length(result$output)]]$Tt[1:bp[ii,jj,1]]
        yy2       <- result$output[[length(result$output)]]$Tt[(bp[ii,jj,1]+1):length(result$Yt)]
        templm1   <- lm(yy1~xx1)
        templm2   <- lm(yy2~xx2) 
        slopelm[ii,jj,1] <- (365/16)*summary(templm1)$coefficients[2,1]
        slopelm[ii,jj,2] <- (365/16)*summary(templm2)$coefficients[2,1]
        sigslope[ii,jj,1]<- summary(templm1)$coefficients[2,4]
        sigslope[ii,jj,2]<- summary(templm2)$coefficients[2,4]
        slope1[1] <- (365/16)*unname(result$output[[length(result$output)]]$Tt[bp[ii,jj,1]] - result$output[[length(result$output)]]$Tt[1])/bp[ii,jj,1]
        slope1[2] <- (365/16)*unname(result$output[[length(result$output)]]$Tt[length(result$Yt)] - result$output[[length(result$output)]]$Tt[bp[ii,jj,1]+1])/(length(result$Yt) - bp[ii,jj,1])
        slope1[3] <- NA
        nobp[ii,jj]<-1
      } 
      else
      {
        xx1       <- c(1:bp[ii,jj,1])
        xx2       <- c((bp[ii,jj,1]+1):(bp[ii,jj,2]))
        xx3       <- c((bp[ii,jj,2]+1):length(result$Yt))
        yy1       <- result$output[[length(result$output)]]$Tt[1:bp[ii,jj,1]]
        yy2       <- result$output[[length(result$output)]]$Tt[(bp[ii,jj,1]+1):(bp[ii,jj,2])]
        yy3       <- result$output[[length(result$output)]]$Tt[(bp[ii,jj,2]+1):length(result$Yt)]
        templm1   <- lm(yy1~xx1)
        templm2   <- lm(yy2~xx2)
        templm3   <- lm(yy3~xx3)
        slopelm[ii,jj,1] <- (365/16)*summary(templm1)$coefficients[2,1]
        slopelm[ii,jj,2] <- (365/16)*summary(templm2)$coefficients[2,1]
        slopelm[ii,jj,3] <- (365/16)*summary(templm3)$coefficients[2,1]
        sigslope[ii,jj,1]<- summary(templm1)$coefficients[2,4]
        sigslope[ii,jj,2]<- summary(templm2)$coefficients[2,4]
        sigslope[ii,jj,3]<- summary(templm3)$coefficients[2,4]
        slope1[1] <- (365/16)*unname(result$output[[length(result$output)]]$Tt[bp[ii,jj,1]] - result$output[[length(result$output)]]$Tt[1])/bp[ii,jj,1]
        slope1[2]<- (365/16)*unname(result$output[[length(result$output)]]$Tt[bp[ii,jj,2]+1] - result$output[[length(result$output)]]$Tt[bp[ii,jj,1]+1])/(bp[ii,jj,2] - bp[ii,jj,1]+1)
        slope1[3]<-(365/16)*unname(result$output[[length(result$output)]]$Tt[-1] - result$output[[length(result$output)]]$Tt[bp[ii,jj,2]+1])/(length(result$Yt) - bp[ii,jj,2] -1)
        nobp[ii,jj]<-2
      }
      slope[ii,jj,]<- slope1
    }
  }
}
#save(file='E:/Caldrought/RegionalBFAST_result_QA01_04222020.RData',list=c('slope','bp','latroi','lonroi','mag','nobp','slopelm','sigslope','bpconf'))
writeMat(output_filename,slope=slope,bp=bp,latroi=latroi,lonroi=lonroi,mag=mag,nobp=nobp,slopelm=slopelm,sigslope=sigslope,bpconf=bpconf)

})
print(stime)