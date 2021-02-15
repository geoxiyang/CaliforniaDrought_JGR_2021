# CaliforniaDrought_JGR_2021
Source code for the JGR paper:

Recovery: Fast and Slow - vegetation response during the 2012-2016 California Drought.

Step0_ReadMODISnetcdf.m: read MODIS EVI and Land Cover data. For EVI data, they need to be cut into pieces as the entire dataset is too large to store in Matlab;<br/>
Step1_RegridPRISM.m: PRISM data are regridded to the same resolution;<br/>
Step2_RegionalBFAST.R: Run BFAST on each pixel. We ran the code for each subsection of the EVI data on a cluster;<br/>
Step3_RegionalBFAST_ReadResult.R: Read the result from Step 2 and put it into a format that is useful for random forest analysis;<br/>
Step4_SpatialVariationExplorer500m.m: Assemble all the data into one giant array.<br/>
Step5_Figures_And_RandomForestAnalysis.ipynb: Random forest analysis and figure making.<br/>

Note that you need R, Matlab, and Python to run all code.

Please contact xiyang@virginia.edu if you have any questions.
