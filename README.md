# QUAK_instructions

## Fitting

Right now, the main fitting software is based in the following repository and branch: https://github.com/wpmccormack/CASEUtils/tree/add_covariance.  This is derived from my fork of https://github.com/case-team/CASEUtils.

To set this up, I suggest that you clone the CASEUtils repository (probably my fork and branch) into a pull of CMSSW.  I have my CASEUtil directory in a pull of CMSSW_10_2_13.  Specifically, I went into CMSSW_10_2_13/src before cloning CASEUtils.  Putting the code in CMSSW_10_2_13 should give you all of the correct versions of packages that you need, e.g. python, ROOT, numpy, etc.

Fitting is then performed in CASEUtils/fitting.  You'll need to pull in a few files from this github repo to steer the QUAK fitting procedure.

First, you'll need to copy QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights.py into CASEUtils/fitting.  This is an ugly macro that does the event selection and then runs the fitting procedure.

You'll also need to put together a directory inside CASEUtils/fitting to hold the data root files (or at least links to the data files) and the plots and fit results.  Please see below for a little section on making data root files.  Assuming you have the files (which in the current format are a separate background and signal file), you should do something like

```
mkdir data_file_folder
cd data_file_folder
```
I normally choose a slightly more descriptive name than data_file_folder, depending on what exactly I'm doing or testing.  The name does NOT have to be data_file_folder explicitly.
Then you need to either copy in the signal and background root files or make a symbolic link to them within data_file_folder.  I recommend the symbolic link to save space.  Here, the naming DOES matter.  The background file needs to be BKG.root, and the signal file needs to be myOutput.root.  Then you can cd up one directory to get back to the fitting directory.

To run the fitting code, you need to run with a command like:
```
python QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights.py -d Dec18_ptModel/xyy_100_normal/ -t sigTemplateMakerWithInterpolation_XToYYprimeTo4Q_MY400_MYprime400/ -l 3000 -r 300
```
Here, my data (simulation in this case) root files were in Dec18_ptModel/xyy_100_normal/.  As you can see, you can have the files in subdirectories if you want.  (Here the signal file was 100 XYY events, and I wanted to give the directory a name that made that obvious).  The directory with the root files is specified with the -d option.

The -t option specifies where the signal templates are.
IMPORTANT NOTES:
* You need to have the signal templates available locally (or at least links to the files).  The README here https://github.com/case-team/CASEUtils/tree/master/fitting has information on how to make the templates.  Also, I think everyone else has made templates for the various signals, so you can borrow them if you need to.  I have mine in /eos/uscms/store/user/wmccorma/CASE_SIGNAL_TEMPLATES/ for example.
* The template file names right now are in the format case_interpolation_M3600.0.root.  In older versions, sometimes the names is graviton_interpolation_M3600.0.root.  If the templates that you're using are in the legacy format, please just uncomment the line in QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights.py which reads self.templateHeader = 'graviton'.

The -l option specifies the mass hypothesis.  You can scan over a range of masses by specifying an upper bound of the scan with an -u option.  For example, if you did "-l 3000 -u 3400", then you would run fits for mass hypotheses of 3000, 3100, 3200, and 3300.  Right now the step size is 100 GeV, and the upper bound is not inclusive.  You need a template for each mass you try to fit for.

The -r option specifies the minimum number of events required to be selected in the signal mass window.  We are using an event selection algorithm that is in truth a HISTOGRAM BIN selection algorithm.  There is a procedure that loops over the least populated bins in the QUAK space for mass sidebands.  For each of these bins, we add the number of events in that bin from the signal mass window to a running total.  Once that value is equal to or greater than than the value specified by the -r option, then the loop terminates.  Note, technically for higher mass hypotheses, we select more events.  The true stopping point is the maximum of max(self.reqEvs, (self.reqEvs/1500.)*500.*np.power((self.masspoint/1000.)-3, 3)), where self.reqEvs is the value specified by the -r parameter.
VERY IMPORTANT NOTE:  For simulated samples, I use -r of 300, and for data, I use -r 1500.


### Other important options

* You can check the end of the QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights.py macro for more options that can be added when you run
* When you make a root file for actual data, the data itself goes in the BKG.root file, and there is no signal.  In this case, you should use the -g option.
* There are options to do various tests like the standard bias tests or tests of efficiency due to JME uncertainties
* You can use -v to make the fitting verbose
* When running on a batch system for getting limits, there are various names that need need to be specified to extract e.g. efficiency values from json files.  These names can also be specified with options


### Other fitting macros

I've included two other fit-steering macros in this repo.  They are both very similar to the main one that I've referenced above, but slightly modified for different purposes.  QUAK_Space_NNTester_LundWeights.py can be used for fitting events based on a neural-net event classifier.  This just avoids doing the event-selection procedure.  I've also used QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights_selectionRD.py to experiment with the event selection procedure itself.


## Making data ROOT files

To make data (and simulation) ROOT files, you'll want to use the macros in /uscms/home/sbrightt/nobackup/CASE/final_oldTraining/ on lpc.  I have some lightly edited versions of the macros: runQuakSpace_editJul25.py and quakSpaceMaker_v2_editJul25.py, which I typically use, and have included in this repo.  I also included versions that can be used to create root files that either don't have systematics, or can included pytorch models that I trained.

You use a command like
```
python runQuakSpace_editJul25.py -b BKG_mjjFlat -s M170-170 -s M170-400 -s M400-400 -s M80-170 -s M80-400 -s M80-80 -r signedL5 -i Wp3000_B400_T170 -d 1 -l lossEvals_pcaDecorr_boxCox_local_Jul24/ -n 1000 -m 0
```
where -s specifies the signal axis (or axes if you use more than one), -r specifies the reduction, -i specifies the signal name, the h5 files are in lossEvals_pcaDecorr_boxCox_local_Jul24, -n specifies the number of signal events to inject, and -m specifies the background type.  -m of 2 runs on actual data.  I've been using the h5 files that are in /eos/uscms/store/user/sbrightt/CASE/final_trainings_bugFix_July2023/lossEvals_pcaDecorr_boxCox/ on lpc.


## Running on batch systems

For getting limits for a large number of signals, you're going to want to run on a batch system.  I normally use LPC's.



PUT FILES IN EOS THAT MIGHT GO AWAY IN THIS GIT REPO!

ALSO specify how to get efficiencies?