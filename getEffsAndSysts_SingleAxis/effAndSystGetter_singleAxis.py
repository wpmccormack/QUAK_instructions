import os
import optparse
import json

parser = optparse.OptionParser()
parser.add_option("-k","--key", dest="key", type=int, default=-1, help="key of dict")


(options, args) = parser.parse_args()

with open('sampleNamesEnhanced.json', 'r') as openfile:
    json_object = json.load(openfile)

print(json_object[str(options.key)])

templateTop = os.getcwd()

os.system("mkdir "+json_object[str(options.key)]["signalName"]+"_shapes")
os.chdir(json_object[str(options.key)]["signalName"]+"_shapes")


for i in range(1600, 6000, 100):
    os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_TEMPLATES/"+json_object[str(options.key)]["signalName"]+"_shapes/"+"case_interpolation_M"+str(i)+".0.json "+".")
    os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_TEMPLATES/"+json_object[str(options.key)]["signalName"]+"_shapes/"+"case_interpolation_M"+str(i)+".0.root "+".")

os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_TEMPLATES/"+json_object[str(options.key)]["signalName"]+"_shapes/"+"case_interpolation.json "+".")
os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_TEMPLATES/"+json_object[str(options.key)]["signalName"]+"_shapes/"+"case_interpolation.root "+".")

os.chdir(templateTop)


os.system("mkdir runningEffDir")
os.chdir("runningEffDir")


#os.system("xrdcp -s root://cmseos.fnal.gov//store/user/sbrightt/CASE/final_trainings_bugFix_July2023/injections/transform-None_reduce-signedL5/BKG_mjjFlat-M170-170-M170-400-M400-400-M80-170-M80-400-M80-80/inject-"+json_object[str(options.key)]["samSignalName"]+"_N1000.root ./myOutput.root")
#os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_EFFICIENCIES_AND_SYSTEMATICS/BKG.root ./BKG.root")
#os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_AUGUST30_SINGLEAXIS_ROOTFILES/"+json_object[str(options.key)]["samSignalName"]+"/inject-"+json_object[str(options.key)]["samSignalName"]+"_N1000.root ./myOutput.root")
os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_AUGUST30_SINGLEAXIS_ROOTFILES/"+json_object[str(options.key)]["samSignalName"]+"/inject-"+json_object[str(options.key)]["samSignalName"]+"_N200.root ./myOutput.root")
os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_AUGUST30_SINGLEAXIS_ROOTFILES/"+json_object[str(options.key)]["samSignalName"]+"/BKG.root ./BKG.root")
os.chdir(templateTop)

#python QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights.py -d runningEffDir/ -t sigTemplateMakerWithInterpolation_XToYYprimeTo4Q_MY80_MYprime170 -l 3000 -r 300 -f YhhTest -k 10 -q 1000 -i
#os.system("python QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights.py -d runningEffDir/ -t "+json_object[str(options.key)]["signalName"]+"_shapes"+" -l "+str(json_object[str(options.key)]["resMass"])+" -r 300 -f "+json_object[str(options.key)]["samSignalName"]+" -k "+str(options.key)+" -q "+str(1000)+" -i > /dev/null")
os.system("python QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights.py -d runningEffDir/ -t "+json_object[str(options.key)]["signalName"]+"_shapes"+" -l "+str(json_object[str(options.key)]["resMass"])+" -r 300 -f "+json_object[str(options.key)]["samSignalName"]+" -k "+str(options.key)+" -q "+str(200)+" -i > /dev/null")

"""
os.system("mkdir runningSystDir")
os.chdir("runningSystDir")


#os.system("xrdcp -s root://cmseos.fnal.gov//store/user/sbrightt/CASE/final_trainings_bugFix_July2023/injections/transform-None_reduce-signedL5/BKG_mjjFlat-M170-170-M170-400-M400-400-M80-170-M80-400-M80-80/inject-"+json_object[str(options.key)]["samSignalName"]+"_N4000.root ./myOutput.root")
#os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_EFFICIENCIES_AND_SYSTEMATICS/BKG.root ./BKG.root")
os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_AUGUST30_SINGLEAXIS_ROOTFILES/"+json_object[str(options.key)]["samSignalName"]+"/inject-"+json_object[str(options.key)]["samSignalName"]+"_N4000.root ./myOutput.root")
os.system("xrdcp -s root://cmseos.fnal.gov//store/user/wmccorma/CASE_AUGUST30_SINGLEAXIS_ROOTFILES/"+json_object[str(options.key)]["samSignalName"]+"/BKG.root ./BKG.root")
os.chdir(templateTop)
"""

#os.system("python QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights.py -d runningSystDir/ -t "+json_object[str(options.key)]["signalName"]+"_shapes"+" -l "+str(json_object[str(options.key)]["resMass"])+" -r 300 -f "+json_object[str(options.key)]["samSignalName"]+" -k "+str(options.key)+" -q "+str(4000)+" -x -w > /dev/null")
os.system("python QUAK_Space_SignalTester_BKGTemplater_Polar_OzFits_ClassBased_LundWeights.py -d runningEffDir/ -t "+json_object[str(options.key)]["signalName"]+"_shapes"+" -l "+str(json_object[str(options.key)]["resMass"])+" -r 300 -f "+json_object[str(options.key)]["samSignalName"]+" -k "+str(options.key)+" -q "+str(200)+" -x -w > /dev/null")


#os.system("xrdfs root://cmseos.fnal.gov/ mkdir /store/user/wmccorma/CASE_SIGNAL_EFFICIENCIES_AND_SYSTEMATICS/"+json_object[str(options.key)]["samSignalName"])
os.system("xrdfs root://cmseos.fnal.gov/ mkdir /store/user/wmccorma/CASE_SIGNAL_EFFICIENCIES_AND_SYSTEMATICS/"+json_object[str(options.key)]["samSignalName"]+"_singleAxis")

os.system("xrdcp -fs runningEffDir/"+json_object[str(options.key)]["samSignalName"]+"_selectionEfficiency.json root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_EFFICIENCIES_AND_SYSTEMATICS/"+json_object[str(options.key)]["samSignalName"]+"_singleAxis"+"/")
#os.system("xrdcp -fs runningSystDir/"+json_object[str(options.key)]["samSignalName"]+"_signalSystematic.json root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_EFFICIENCIES_AND_SYSTEMATICS/"+json_object[str(options.key)]["samSignalName"]+"_singleAxis"+"/")
os.system("xrdcp -fs runningEffDir/"+json_object[str(options.key)]["samSignalName"]+"_signalSystematic.json root://cmseos.fnal.gov//store/user/wmccorma/CASE_SIGNAL_EFFICIENCIES_AND_SYSTEMATICS/"+json_object[str(options.key)]["samSignalName"]+"_singleAxis"+"/")


os.system("ls runningEffDir/")
os.system("ls runningSystDir/")
