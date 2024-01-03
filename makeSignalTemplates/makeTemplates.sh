#!/bin/bash
set -x
echo "pwd"
pwd
echo "ls"
ls

xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CMSSW_10_2_13_LundExplorer.tgz .
source /cvmfs/cms.cern.ch/cmsset_default.sh
tar -xf CMSSW_10_2_13_LundExplorer.tgz
rm CMSSW_10_2_13_LundExplorer.tgz
cd CMSSW_10_2_13/src
echo "pwd where am i"
pwd
echo "ls what's here"


scramv1 b -r ProjectRename # this handles linking the already compiled code - do NOT recompile
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

echo $CMSSW_BASE "is the CMSSW we have on the local worker node"

cd CASEUtils/fitting/

xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/sampleNames.json .
xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/jsonReader.py .

python jsonReader.py -k ${2}

ls | grep shapes

hostname
date
