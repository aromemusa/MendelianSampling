# ```msvalues``` module
The main functions of ```msvalues```module.   
<br/><br/>

The following four functions import example empirical data to show users how to format data and use msvalues module.
### ```msvalues.inputfile1st()```
Imports example data (marker map and effects for a single trait).   
<br/>

### ```msvalues.inputfile1mt()```
Imports example data (marker map and effects for multiple traits).   
<br/>

### ```msvalues.inputfile2ss()```
Imports example data (Identification number, sex (single), and phased genotype data).   
<br/>

### ```msvalues.inputfile2ms()```
Imports example data (Identification number, sex (male and female), and phased genotype data).   
<br/>


### ```msvalues.fileschecks(inputfile1, inputfile2, no_traits, index_wt)```
Check the main data required for most functions for errors.
#### Parameters
* ```inputfile1```: pandas data frame containing marker map and marker effects
* ```inputfile2```: 2-dimensional string array containing identification numbers, sex, and phased genotypic data
* ```no_traits```: an integer of the number of traits
* ```index_wt```: an array of index weights; preferably NumPy array  
<br/><br/>

### ```msvalues.popcovmat(inputfile1, mposunit)```
Sets up population covariance matrix reflecting expected within-family linkage disequilibrium of markers.
#### Parameters
* ```inputfile1```: pandas data frame containing marker map and marker effects
* ```mposunit```: a string (either 'cM' or 'bp')  
<br/><br/>

### ```msvalues.MSvarcov(inputfile1, inputfile2, CovMat, index_wt, progress = True)```
Calculates Mendelian sampling co(variance) and aggregate breeding value
#### Parameters
* ```inputfile1```: pandas data frame containing marker map and marker effects
* ```inputfile2```: 2-dimensional string array containing identification numbers, sex, and phased genotypic data
* ```CovMat```: a list containing chromosome-wise population covariance matrices
* ```index_wt```: an array of index weights; preferably NumPy array
* ```progress```: if ```True``` the progress of the calculations will be printed to screen  
<br/><br/>

### ```msvalues.SelStratGam(selstrat, inputfile1, inputfile2, MSVMSC, index_wt, selinorthresh)```
Calculates selection criterium using the gametic approach
#### Parameters
* ```selstrat```: a string (any of 'GEBV', 'PBTI', and 'RPTA') 
* ```inputfile1```: pandas data frame containing marker map and marker effects
* ```inputfile2```: 2-dimensional string array containing identification numbers, sex, and phased genotypic data
* ```MSVMSC```: pandas data frame containing Mendelian sampling (co)variance and aggregate breeding value
* ```index_wt```: an array of index weights; preferably NumPy array
* ```selinorthresh```: a float variable indicating selection intensity if selstrat is RPTA or value used for threshold calculation if selstrat is PBTI   
<br/><br/>

### ```msvalues.SelStratZyg(selstrat, inputfile1, inputfile2, MSVMSC, index_wt, selinorthresh, selintmale, selintfemale)```
Calculates selection criterium using the zygotic approach
#### Parameters
* ```selstrat```: a string (any of 'GEBV', 'PBTI', and 'RPTA') 
* ```inputfile1```: pandas data frame containing marker map and marker effects
* ```inputfile2```: 2-dimensional string array containing identification numbers, sex, and phased genotypic data
* ```MSVMSC```: pandas data frame containing Mendelian sampling (co)variance and aggregate breeding value
* ```index_wt```: an array of index weights; preferably NumPy array
* ```selinorthresh```: a float variable indicating selection intensity if selstrat is RPTA or value used for threshold calculation if selstrat is PBTI
* ```selintmale```: a float variable indicating the percentage of top males to select
* ```selintfemale```: a float variable indicating the percentage of top females to select   
<br/><br/>

### ```msvalues.SimMatGam(inputfile1, inputfile2, index_wt, CovMat, chrinterest, stdsim = True, progress = True)```
Calculates similarity matrices based on Mendelian sampling values using the gametic approach
#### Parameters
* ```inputfile1```: pandas data frame containing marker map and marker effects
* ```inputfile2```: 2-dimensional string array containing identification numbers, sex, and phased genotypic data
* ```index_wt```: an array of index weights; preferably NumPy array
* ```CovMat```: a list containing population covariance matrix
* ```chrinterest```: a string of either 'all' or 'none', or a list with chromosome(s) of interest e.g., chrinterest = [4, 14] 
* ```stdsim```: if ```True``` standardized similarity matrices for chromosomes of interest will be saved to file and standardized similarity matrix based on aggregate genotype will be output; if ```False```, the same is true for unstandardized similarity matrices
* ```progress```: if ```True``` the progress of the calculations will be printed to screen   
<br/><br/>

### ```msvalues.SimMatZyg(inputfile1, inputfile2, index_wt, ZygoDF, CovMat, chrinterest, stdsim = True, progress = True)```
Calculates similarity matrices based on Mendelian sampling values using the zygotic approach
#### Parameters
* ```inputfile1```: pandas data frame containing marker map and marker effects
* ```inputfile2```: 2-dimensional string array containing identification numbers, sex, and phased genotypic data
* ```index_wt```: an array of index weights; preferably NumPy array
* ```ZygoDF```: pandas data frame containing selection criterium using the zygotic approach
* ```CovMat```: a list containing population covariance matrix
* ```chrinterest```: a string of either 'all' or 'none', or a list with chromosome(s) of interest e.g., chrinterest = [4, 14] 
* ```stdsim```: if ```True``` standardized similarity matrices for chromosomes of interest will be saved to file and standardized similarity matrix based on aggregate genotype will be output; if ```False```, the same is true for unstandardized similarity matrices
* ```progress```: if ```True``` the progress of the calculations will be printed to screen   
<br/><br/>

### ```msvalues.OptiMateAllocGam(data, Similaritymat, maxalloc, progress = True)```
Optimizes mate allocation for gametic approach
#### Parameters
* ```data```: pandas data frame containing selection criterium using the gametic approach
* ```Similaritymat```: 2-dimensional NumPy float array containing similarity matrix
* ```maxalloc```: a float variable indicating the maximum weight to be allotted to an individual. The maximum value is 1.0 
* ```progress```: if ```True``` the progress of the calculations will be printed to screen   
<br/><br/>


### ```msvalues.OptiMateAllocZyg(data, Similaritymat, maxalloc, progress = True)```
Optimizes mate allocation for zygotic approach
#### Parameters
* ```data```: pandas data frame containing selection criterium using the zygotic  approach
* ```Similaritymat```: 2-dimensional NumPy float array containing similarity matrix
* ```maxalloc```: a float variable indicating the maximum weight to be allotted to a mate-pair. The maximum value is 1.0 
* ```progress```: if ```True``` the progress of the calculations will be printed to screen
