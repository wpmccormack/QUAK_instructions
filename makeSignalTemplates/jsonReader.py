import os
import optparse
import json

parser = optparse.OptionParser()
parser.add_option("-k","--key", dest="key", type=int, default=-1, help="key of dict")


(options, args) = parser.parse_args()

with open('sampleNames.json', 'r') as openfile:
    json_object = json.load(openfile)

print(json_object[str(options.key)])


h5s = json_object[str(options.key)]["h5s"]

for h in h5s:
    #print(h)
    os.system("xrdcp root://eosuser.cern.ch//eos/cms/store/group/phys_b2g/CASE/h5_files/UL/merged/"+h+' .')


#python fit_signalshapes.py -M 2000 3000 5000 --inputFiles XToYYprimeTo4Q_MX2000_MY170_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5 XToYYprimeTo4Q_MX3000_MY170_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5 XToYYprimeTo4Q_MX5000_MY170_MYprime25_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5 --dcb-model --fitRange 1.0 -o XToYY_shapes/
os.system("python fit_signalshapes.py -M 2000 3000 5000 --inputFiles "+json_object[str(options.key)]["h5sString"]+" --dcb-model --fitRange 1.0 -o "+json_object[str(options.key)]["signalName"]+"_shapes")
#print("python fit_signalshapes.py -M 2000 3000 5000 --inputFiles "+json_object[str(options.key)]["h5sString"]+" --dcb-model --fitRange 1.0 -o "+json_object[str(options.key)]["signalName"]+"_shapes")




#python interpolation.py -i XToYY_shapes/full_fit.root -s case -o XToYY_shapes --masses 1600 1700 1800 1900 2100 2200 2300 2400 2500 2600 2700 2800 2900 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5100 5200 5300 5400 5500 5600 5700 5800 5900
os.system("python interpolation.py -i "+json_object[str(options.key)]["signalName"]+"_shapes"+"/full_fit.root -s case -o "+json_object[str(options.key)]["signalName"]+"_shapes"+" --masses 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5000 5100 5200 5300 5400 5500 5600 5700 5800 5900")
#print("python interpolation.py -i "+json_object[str(options.key)]["signalName"]+"_shapes"+"/full_fit.root -s case -o "+json_object[str(options.key)]["signalName"]+"_shapes"+" --masses 1600 1700 1800 1900 2100 2200 2300 2400 2500 2600 2700 2800 2900 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5100 5200 5300 5400 5500 5600 5700 5800 5900")


os.system("xrdfs root://cmseos.fnal.gov/ mkdir /store/user/wmccorma/CASE_SIGNAL_TEMPLATES/"+json_object[str(options.key)]["signalName"]+"_shapes")

for i in range(1600, 6000, 100):
    os.system("xrdcp -s "+json_object[str(options.key)]["signalName"]+"_shapes/"+"case_interpolation_M"+str(i)+".0.json root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_TEMPLATES/"+json_object[str(options.key)]["signalName"]+"_shapes/")
    os.system("xrdcp -s "+json_object[str(options.key)]["signalName"]+"_shapes/"+"case_interpolation_M"+str(i)+".0.root root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_TEMPLATES/"+json_object[str(options.key)]["signalName"]+"_shapes/")

os.system("xrdcp -s "+json_object[str(options.key)]["signalName"]+"_shapes/"+"case_interpolation.json root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_TEMPLATES/"+json_object[str(options.key)]["signalName"]+"_shapes/")
os.system("xrdcp -s "+json_object[str(options.key)]["signalName"]+"_shapes/"+"case_interpolation.root root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_TEMPLATES/"+json_object[str(options.key)]["signalName"]+"_shapes/")

os.system("ls "+json_object[str(options.key)]["signalName"]+"_shapes/")
