This file explains useage of the Chaste bolt-on user project PlosOne_mRNA

This file is divided into four parts

1. Archive contents

2. Using the Chaste code

3. Supercomputing

4. Matlab files

%------------------------------------------------------------------------------
1. Archive contents

The Chaste bolt-on project you have downloaded should either be copied into the 'projects' folder of your Chaste folder or ideally given a symbolic link to the Chaste projects folder. Thus the directory structure should read
Chaste/projects/PlosOne_mRNA/test (etc).

The Chaste project contains six folders:

test: Contains the C++ file TestSensitivityAnalysisOHaraEndoLiteratePaper.hpp and the folder data, containing the experimental designs run in the paper.

src: Contains SensitivityDataStructure.cpp, SensitivityDataStructure.hpp, and ohara_rudy_2011.cellml

build: Currently empty. Once you compile TestSensitivityAnalysisOHaraEndo.hpp, binaries will go here.

matlab: Contains code used to analyse and plot the biomarkers in the study. Not used by Chaste.

output: Contains results obtained from our experiments.

pbs: Batch submission scripts used to simulate the model population on the Oxford Supercomputing Centre's machines.

%------------------------------------------------------------------------------
 
2. Using the Chaste code

This project is known to work with Chaste release 3.1. Note that you must have CVode working with Chaste for ths project to work.

To use the Chaste code, it is recommended you familiarise yourself with how Chaste works on the Chaste wiki.

The file to be compiled, TestSensitivityAnalysisOHaraEndoLiteratePaper.hpp, may be compiled from your Chaste folder as follows:

cd Chaste

scons ts=projects/PlosOne_mRNA/test/TestSensitivityAnalysisOHaraEndoLiteratePaper.hpp

scons will then compile and run the resulting executable. This will result in a request for the following arguments:
# 0 arguments supplied.
TestSensitivityAnalysis::Please input an argument
* --file  the path of an experimental design file
* --tag  tag for identifying output (optional)
* --store-data  record voltage and calcium transients (optional)
Passed
OK!

To add these arguments, we can use scons and run_time_flags (see below)

Arguments:
--file is a path to an experimental design file. This can be either absolute (ie relative to the root directory of the computer) or relative to the Chaste folder. These are formatted as follows:
1st line: number of experiments in file
2nd line: number of parameters to be modified followed by parameter mask (see src/SensitivityDataStructure.cpp for details)
3+ line: number of experiment, followed by parameter scaling for each parameter in mask.

Files used in this experiment are contained in the folder test/data/ 
uniform_30: The non-failing experiment.
uniform_3060: The failing experiment.
surf_30 and surf_3060: used for generating the surfaces in the supplement

--tag adds a string to the output in your Chaste test output folder to help you identify the output

--store-data stores the last two action potentials in each case in your chaste output folder for validation and plotting purposes. Note: this gives a lot of data very quickly, so be careful about how many you print at once! Consider also changing the printing time step in TestSensitivityAnalysisOHaraEndo.hpp

So, to run and tag all your output with test, load an example failing heart experiment (containing 25 experiments) and store the resulting final APs, do:

scons ts=projects/PlosOne_mRNA/test/TestSensitivityAnalysisOHaraEndo.hpp run_time_flags="--file projects/PlosOne_mRNA/test/data/input_uniform_3060/exp_design_uniform_3060_1.dat --tag test --store-data"

Alternatively, we may use the compile argument compile_only=1 (or co=1) to compile but not run Chaste, then once this is done we can run the executable with flags as normal. This will be useful when using supercomputers (see 4):

scons co=1 ts=projects/PlosOne_mRNA/test/TestSensitivityAnalysisOHaraEndo.hpp
projects/PlosOne_mRNA/build/debug/TestSensitivityAnalysisOHaraEndoRunner --file projects/PlosOne_mRNA/test/data/input_uniform_3060/exp_design_uniform_30_1.dat --tag test --store-data

Other useful arguments for scons:
cl=1 (compile chaste as libraries)
b=GccOpt (use optimised GCC compiler)
b=GccOptNative (use optimised GCC compiler native to this architecture)
b=Intel (use Intel C++ compiler, if available)

The output files have the following name: DYN_(experiment number).out, and appear in a folder in your Chaste test output directory called SensitivityAnalysis/OHara2011_endo_{tag}, where tag is supplied with --tag as above.
They are formatted as follows:
column 1: BCL
column 2: APD30 (beat one)
column 3: APD30 (beat two)
column 4: APD80 (beat one)
column 5: APD80 (beat two)
column 6: CaTD30 (beat one)
column 7: CaTD30 (beat two)
column 8: CaTD80 (beat one)
column 9: CaTD80 (beat two)
column 10: DAPCaT (beat one)
column 11: DAPCAT (beat two)
column 12: CaTmax (beat one)
column 13: CaTmax (beat two)


%------------------------------------------------------------------------------

3. Supercomputing

All files for supercomputing provided are provided for example only - we do not intend to provide support for this feature. Queries may be addressed to johnwalmsley@gmail.com, however. We use the pbs system and as such if your system is set up differently then the file format will be different to that shown. Similarly, all directories will need changing unless you happen to be a member of the Rodriguez group using HAL!

To run on the OSC system we performed the following:

1) Compile TestSensitivityAnalysisOHaraEndo on the compile node using the intel compiler:

scons co=1 b=Intel ts=projects/PlosOne_mRNA/test/TestSensitivityAnalysisOHaraEndo.hpp

2) Run the batch file for the appropriate number of experiments. Note that we loop in groups of eight for efficiency, as the OSC's HAL machine has 8 cores on each compute node:
(eg)
qsub -J 1-656:8 batchscript_batchscript_OHaraEndo_uniform_3060.sh

More information on PBS can be found at http://www.osc.ox.ac.uk/content/pbs

%------------------------------------------------------------------------------

4. Matlab files

All matlab files were written by John Walmsley, Gary Mirams and Jose-Felix Rodriguez.

These files are not supported by the Chaste team - queries may be sent to John Walmsley at johnwalmsley@gmail.com. If you wish to re-use this code, please either add John Walmsley as an acknowledgement and/or cite our mRNA expression paper in PLOS ONE, Walmsley et al 2013 as appropriate. We use the freely available matlab program export_fig to generate our figures as it offers vastly improved "what you see is what you get" over standard Matlab http://www.mathworks.co.uk/matlabcentral/fileexchange/23629-exportfig.

You will need to change the directories used in all these files to get them to work on your machine. Once the directories have been changed as appropriate (you may want to copy them out of your Chaste folder), the figures used in our paper should be generated.

The output of the files we obtained in our experiments can be found in the output folder. uniform_30 is the results for the non-failing populaiton, uniform_3060 is the results for the fialing population, and surface_30 and surface_60 are the results used to construct the parameter surfaces in the supplement.

You will need to change the directories at the head of the file so that they are appropriate for your machine.

Files:

exp_des_uniform.m 
Used to generate experimental designs using a uniform distribution.

exp_des_surface.m
Used to generate experimental designs for the parameter surface in the supplement.

post_proc.m 
This file performs post-processing of the biomarkers (detecting alternans, dividing through APD30 by APD80 to get APD3080, etc). 

statistics_Chaste.m
Generates the histograms seen in Fig. 3 as well as generating summary statistics from the distributions. To get the baseline and 60s files needed by the analysis program, you'll need to generate these files by running the Chaste executable in section 2 on the file test/data/BaselineAndExtremeCase.dat. You'll then need to rename the output appropriately! Using a sensible tag for your output is recommended to avoid overwriting other data.

statistics_correlations_Chaste.m
Generates the correlation plots in Fig. 4. Requires JW_plotmatrix_triangle.m, an overwritten form of the matlab plotmatrix command.

regressionanalysis.m
Performs linear regression analysis on the data provided, to generate the plots in Fig. 5. Also postprocesses the alternans data.

ParameterSurface.m 
Generates the parameter surface figures in the supplement.
 










