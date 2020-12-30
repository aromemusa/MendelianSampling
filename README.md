# MSvaluesPy
MSvaluesPy is a python package for estimating Mendelian sampling-related quantities, such as variation, covariance, and similarity between potential parents. The package implements a computationally faster approach for estimating Mendelian sampling co(variance) and the novel similarity matrices based on Mendelian sampling values [1]. 

Also, the package estimates three selection criteria, namely genomic estimated breeding value (GEBV; [2]), relative predicted transmitting ability (RPTA; [3]) and the probability of breeding top-ranking individuals (PBTI; [1]). Finally, this package optimizes mate allocation and decisions based on a selection criterion using the novel similarity matrices following the efficient frontier. The efficient frontier provides breeders with optimized portfolios containing parents/matings that provide the highest possible genetic return with the lowest degree of haplotype similarity aiming to maintain diversity by keeping different heterozygous chromosome regions active in the population.

For simplicity, MSvaluesPy has only one module named [`msvalues`](docs/docs_msvalues.md): functions 

To use [`msvalues`](docs/docs_msvalues.md) module, include the following code in your python file:

`from MSvaluesPy import msvalues`


## Installation
MSvaluesPy can be installed using `pip` with the command:

`pip install MSvaluesPy`

To install the latest software version, the following command should be issued:

`pip install MSvaluesPy==1.3.2`

## Example
Examples are given as interactive python notebook (ipynb) files:
* [`demo_analysis.ipynb`](demo_analysis.ipynb): Analysis of example dataset showing the format of required input data and how to use the functions in the [`msvalues`](docs/docs_msvalues.md) module


## How to cite MSvaluesPy


## References
1. Musa AA, Reinsch N. Similarities of Mendelian sampling values between families and their application in hedging haplotype diversity. Genet Sel Evol. 2021;53.
2. Meuwissen THE, Hayes BJ, Goddard ME. Prediction of total genetic value using genome-wide dense marker maps. Genetics. 2001;157:1819–29.
3. Santos DJA, Cole JB, Lawlor TJ, VanRaden PM, Tonhati H, Ma L. Variance of gametic diversity and its application in selection programs. J Dairy Sci. 2019;102:5279–94.
