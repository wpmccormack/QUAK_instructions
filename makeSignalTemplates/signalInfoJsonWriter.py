import json

signals = {
    "QstarToQW": {
        "mainMass": "M_",
        "resonances": ["mW_"],
        "resonanceMasses": [[25, 80, 170, 400]],
        "tune": "TuneCP2_13TeV-pythia8_TIMBER",
        "samName": "Qstar",
        "samRes": ["W"],
        "samExtra": ""
    },
    "RSGravitonToGluonGluon_kMpl01": {
        "mainMass": "M_",
        "resonances": [],
        "resonanceMasses": [[]],
        "tune": "TuneCP5_13TeV_pythia8_TIMBER",
        "samName": "",
        "samRes": [],
        "samExtra": ""
    },
    "WkkToWRadionToWWW": {
        "mainMass": "M",
        "resonances": ["Mr"],
        "resonanceMasses": [[170, 400]],
        "tune": "TuneCP5_13TeV-madgraph-pythia8_TIMBER",
        "samName": "Wkk",
        "samRes": ["R"],
        "samExtra": ""
    },
    "WpToBpT": {
        "mainMass": "Wp",
        "resonances": ["Bp"],
        "resonanceMasses": [[25, 80, 170, 400]],
        "tune": "Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER",
        "samName": "Wp",
        "samRes": ["B"],
        "samExtra": "_T170"
    },
    "XToYYprimeTo4Q": {
        "mainMass": "MX",
        "resonances": ["MY", "MYprime"],
        "resonanceMasses": [[25, 80, 170, 400], [25, 80, 170, 400]],
        "tune": "narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER",
        "samName": "XYY",
        "samRes": ["Y", "Yp"],
        "samExtra": ""
    },
    "YtoHH_Htott": {
        "mainMass": "Y",
        "resonances": ["H"],
        "resonanceMasses": [[400]],
        "tune": "TuneCP5_13TeV-madgraph-pythia8_TIMBER",
        "samName": "YHH",
        "samRes": ["H"],
        "samExtra": ""
    },
    "ZpToTpTp": {
        "mainMass": "Zp",
        "resonances": ["Tp"],
        "resonanceMasses": [[400]],
        "tune": "TuneCP5_13TeV-madgraph-pythia8_TIMBER",
        "samName": "ZTT",
        "samRes": ["Tp"],
        "samExtra": ""
    }
}



outputDict = {}
outputDictEnhanced = {}
numKey = 0
numKeyEnhanced = 0
for s in signals:

    totalDecays = 1
    for r in range(len(signals[s]["resonances"])):
        totalDecays *= len(signals[s]["resonanceMasses"][r])

    combosout = ['']*totalDecays
    samCombos= ['']*totalDecays

    combosused = 1
    for r in range(len(signals[s]["resonances"])):
        for m in range(len(signals[s]["resonanceMasses"][r])):
            for x in range(totalDecays/len(signals[s]["resonanceMasses"][r])):
                combosout[ combosused*x + m*(totalDecays/len(signals[s]["resonanceMasses"][r]))/combosused ] += signals[s]["resonances"][r]+str(signals[s]["resonanceMasses"][r][m])
                samCombos[ combosused*x + m*(totalDecays/len(signals[s]["resonanceMasses"][r]))/combosused ] += signals[s]["samRes"][r]+str(signals[s]["resonanceMasses"][r][m])
                if(r != len(signals[s]["resonances"])-1):
                    combosout[ combosused*x + m*(totalDecays/len(signals[s]["resonanceMasses"][r]))/combosused] += '_'
                    samCombos[ combosused*x + m*(totalDecays/len(signals[s]["resonanceMasses"][r]))/combosused] += '_'
        combosused *= len(signals[s]["resonanceMasses"][r])

    #print(s,combosout)
    
    for c in range(len(combosout)):
        if(combosout[c]!=''):
            signalname = s+'_'+combosout[c]
        else:
            signalname = s
        h5sString = ''
        h5s = []
        for b in ["2000", "3000", "5000"]:
            if(combosout[c]!=''):
                h5sString += s+'_'+signals[s]["mainMass"]+b+'_'+combosout[c]+'_'+signals[s]["tune"]+'.h5'
                h5s.append(s+'_'+signals[s]["mainMass"]+b+'_'+combosout[c]+'_'+signals[s]["tune"]+'.h5')
                samSignalName = signals[s]["samName"]+str(b)+'_'+samCombos[c]+signals[s]["samExtra"]
                specificSignalName = s+'_'+signals[s]["mainMass"]+b+'_'+combosout[c]+'_'+signals[s]["tune"]
            else:
                h5sString += s+'_'+signals[s]["mainMass"]+b+'_'+signals[s]["tune"]+'.h5'
                h5s.append(s+'_'+signals[s]["mainMass"]+b+'_'+signals[s]["tune"]+'.h5')
                samSignalName = signals[s]["samName"]+str(b)+signals[s]["samExtra"]
                specificSignalName = s+'_'+signals[s]["mainMass"]+b+'_'+signals[s]["tune"]
            if(b != "5000"):
                h5sString += ' '

            outputDictEnhanced[numKeyEnhanced] = {"signalName": signalname, "samSignalName": samSignalName, "resMass": int(b), "specificSignalName": specificSignalName}
            numKeyEnhanced+=1

        outputDict[numKey] = {"signalName": signalname, "h5sString": h5sString, "h5s": h5s}

        numKey+=1

#print(outputDict)
with open("sampleNames.json", "w") as outfile:
    json.dump(outputDict, outfile)

with open("sampleNamesEnhanced.json", "w") as outfile:
    json.dump(outputDictEnhanced, outfile)
