# Meta-analysis on global change impacts on soil biodiversity

This repository holds the code for the analysis presented in the paper 
"Global change and their environmental stressors have a significant impact on soil biodiversity – a meta-analysis".

The data for the analysis is also available: https://doi.org/10.5281/zenodo.6903152

## Overview

Scripts 01_, 02_ and 04_ (there is no 03_, I thought there would be, but it wasn't needed) prepare the dataset.
cleaning the data, harmonising species names, as well as calculating the effect sizes (amongst other things). 
Although the dataset is saved later, these scripts are what produced the dataset in the Zenodo link. The rest of the 
scripts do not change the data, just subset it.

Scripts 5_ and 6_x then perform the analyses for the different models (publication bias models are performed in each individual
script as well). The dataset from the Zenodo link is needed for these scripts.

Script 7_ then produces all the figures in the main manuscript as well the supplementary material. And script 8_ fixes the
dataset ready for exporting to the Zenodo link.

## Additional Information

All packages and additional functions are listed at the top of each script. 

## Contact info

Any questions or issues, please feel free to contact me:
helen[dot]phillips[at]smu[dot]ca

