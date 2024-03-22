### Description
R codes to reproduce the results of leveraging external Gaussian endpoint with unknown variance from a single external control arm in the article: "Modified Conditional Borrowing-by-part Power Prior for the dynamic and parameter-specific information borrowing of Gaussian endpoint" by Kai Wang, Han Cao, Chen Yao. 


### Files and subfolders
*1  "Function_for_MultipleECs.R" contains the R functions to be loaded.
*2  "Numerical_study_for_MultipleECs.R" is the code to conduct the numerical study (the working directory must be       
    the directory of the file "Code_for_MultipleECs").
*3  "Simulation_study_for_MultipleECs.R" is the code to conduct the simulation study (the working directory must   
    be the directory of the file "Code_for_MultipleECs").
*4  Folder "stan" contains the stan codes of full Bayesian dynamic borrowing methods for leveraging external Gaussian endpoint 
      with unknown variance from multiple (three) external control arms.
*5  Folder "simulation_dataset" is used to save the simulated datasets
*6  Folder "rst_for_MultipleECs" is used to save the results of the numerical and simulation study.



### Platform
R version 4.1.3 runs on Linux  
Stan version 2.21.0, called via the R package rstan  

