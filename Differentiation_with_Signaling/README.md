"NDD_with_Signaling.m": The main m.file which produces the main results, FIG.3, and movies. (MATLAB 2016) 

"shoving.biomass": JAVA code to evaluate Shoving step.

"gene_expression": MATLAB function to evaluate Gillespie algorithm. 

"clust_coeff" : MATLAB function to stimate clustering coefficient. 

"simple_spectral_partitioning": MATLAB function to stimate comminity size. 


To run the code, COMSOL 5.2 should be installed. The "COMSOL_Installation_Instructions.pdf" file gives the installation instructions. It requires a license file. The installation would take few minutes. 

The code will generate a figure which shows cell population.   

Data and code associated with:

Safdari et al. Noise-driven Cell Differentiation and the Emergence of Organization. bioRxiv doi: https://doi.org/10.1101/220525

## How to run JAVA codes:

Instructor to use static path definition (instead of this dynamic path definition) to run JAVA codes.

 
1. In MATLAB command window  write,

       cd(prefdir)

to go to the preferences directory of MATLAB.

2. Then type,

      edit javaclasspath.txt

and this file will be created because it does not exist yet.

3. In the new (empty) file “javaclasspath.txt” write the path for the JAVA codes:

   FOLDER PATH\java\shoving\build\classes\

Then save the file.

4. Close MATLAB and COMSOL server and start again COMSOL with MATLAB.

5. Now the path to the shoving function will be static and you can check it by typing:

     javaclasspath('-static')

the path could be seen.
