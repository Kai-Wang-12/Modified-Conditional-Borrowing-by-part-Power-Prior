### Description
R codes to reproduce the results of leveraging external Gaussian endpoint with unknown variance from a single external control arm in the article: "Modified Conditional Borrowing-by-part Power Prior for the dynamic and parameter-specific information borrowing of Gaussian endpoint" by Kai Wang, Han Cao, Chen Yao. 


### Files and subfolders
*1  "Function_for_SingleEC.R" contains the R functions to be loaded.
*2  "Numerical_study_for_SingleEC.R" is the code to conduct the numerical study (the working directory must be       
      the directory of the file "Code_for_SingleEC").
*3  "Simulation_study_for_SingleEC.R" is the code to conduct the simulation study (the working directory must   
      be the directory of the file "Code_for_SingleEC").
*4  "Application.R" is the code to reproduce the application results in Section 6 of the article (the working    
      directory must be the directory of the file "Code_for_SingleEC").
*5  Folder "stan" contains the stan codes of full Bayesian dynamic borrowing methods for leveraging external Gaussian endpoint      
      with unknown variance.
*6  Folder "simulation_dataset" is used to save the simulated datasets
*7  Folder "rst_for_SingleEC" is used to save the results of the numerical study, simulation study and application.



### Platform
R version 4.1.3 runs on Linux  
Stan version 2.21.0, called via the R package rstan  

