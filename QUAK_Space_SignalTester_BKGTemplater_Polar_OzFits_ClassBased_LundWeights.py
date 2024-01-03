import ROOT as r,sys,math,array,os,sys
#from optparse import OptionParser
#from ROOT import std,RooDataHist
from array import array
import numpy as np
import h5py
import glob
import os
import optparse
import shutil
import numpy as np
#from Utils import load_h5_avgBkgLoss
#from Utils import load_h5_avgSigLoss
#import array as array

import uproot
import numpy as np
#import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import subprocess
import json
import scipy.stats as stats
#import mplhep as hep
#import awkward as ak


#theSysWeights = ['pdf_up', 'pdf_down', 'prefire_up', 'prefire_down', 'pileup_up', 'pileup_down', 'btag_up', 'btag_down', 'PS_ISR_up', 'PS_ISR_down', 'PS_FSR_up', 'PS_FSR_down', 'F_up', 'F_down', 'R_up', 'R_down', 'RF_up', 'RF_down', 'top_ptrw_up', 'top_ptrw_down']
theSysWeights = ['pdf_up', 'pdf_down', 'prefire_up', 'prefire_down', 'pileup_up', 'pileup_down', 'btag_up', 'btag_down', 'PS_ISR_up', 'PS_ISR_down', 'PS_FSR_up', 'PS_FSR_down', 'F_up', 'F_down', 'R_up', 'R_down', 'RF_up', 'RF_down', 'top_ptrw_up', 'top_ptrw_down', 'lund_plane_up', 'lund_plane_down', 'lund_bjet_up', 'lund_bjet_down']


minangle = -.1*np.pi
angularbincount = 25
maxRad = 10
radBins = 20

"""
minangle = -.1*np.pi
angularbincount = 50
maxRad = 10
radBins = 40
"""
#bkgCut = -0.5
bkgCut = -1.
#bkgCut = -3.

#windowOffset = 200
windowOffset = 100
#windowOffset = 150
#windowOffset = 0

sidebandsize = 500
#sidebandsize = 400

binDict = {
    "A0": [[], [1350,1650], []],
    "A1": [[1800, 1900], [1650, 2017], [.01]],
    "A2": [[2200, 2300], [2017, 2465], [.01]],
    "A3": [[2600, 2700, 2800], [2465, 3013], [.01]],
    "A4": [[3200, 3300, 3400, 3500], [3013, 3682], [.03]],
    "A5": [[3900, 4100, 4200, 4300], [3682, 4500], [.03]],
    "A6": [[4800, 4900, 5000, 5100, 5200], [4500, 5500], [.05]],
    "A7": [[], [5500, 8000], []],
    "B0": [[], [1492, 1824], []],
    "B1": [[2000, 2100], [1824, 2230], [.01]],
    "B2": [[2400, 2500], [2230, 2725], [.01]],
    "B3": [[2900, 3000, 3100], [2725, 3331], [.01]],
    "B4": [[3600, 3700, 3800], [3331, 4071], [.03]],
    "B5": [[4400, 4500, 4600, 4700], [4071, 4975], [.03]],
    "B6": [[5300, 5400, 5500, 5700, 5800], [4975, 6081], [.05]],
    "B7": [[], [6081, 8000], []],
}


def massArraySlicer(contents, binval, mass, loss1, loss2, masspoint, contentmax, binstoadd, useWSbins, wsbin, windowsize=500, nBins=50):
    rbins = np.linspace(0, maxRad, radBins)
    abins = np.linspace(minangle,(0.4 + 0.02*binstoadd)*np.pi, angularbincount+binstoadd)
    nullXBounds = []
    nullYBounds = []
    contentsmin = contents.min()
    newmax = contents.max() + 1
    
    if(contentsmin==contentmax):
        print('warning contentsmin')
        return 0,[],True,[]
        
    for a in range(angularbincount+binstoadd-1):
        for r in range(radBins-1):
            if(contents[a,r]==contentsmin):
                
                nullXBounds.append(a)
                nullYBounds.append(r)
                if(contents[a,r]>0):
                    print('first print bin with ',contents[a,r], 'at', abins[a], rbins[r], 'x =', rbins[r]*np.cos(abins[a]), 'y =', rbins[r]*np.sin(abins[a]))
                contents[a,r]=newmax
                
    angles = np.arctan(loss1/loss2)
    radii = np.sqrt(loss1*loss1 + loss2*loss2)
    thebool = []
    print('len(nullXBounds) = ',len(nullXBounds))
    for b in range(len(nullXBounds)):
        #print('loop over nullxbounds b =',b, 'nullXBounds[b] = ',nullXBounds[b], 'abins[nullXBounds[b]] = ', abins[nullXBounds[b]], 'abins[nullXBounds[b]+1] = ',abins[nullXBounds[b]+1], 'nullYBounds[b] =', nullYBounds[b], 'rbins[nullYBounds[b]] =', rbins[nullYBounds[b]], 'rbins[nullYBounds[b]+1] =',rbins[nullYBounds[b]+1])
        abool = (angles > abins[nullXBounds[b]]) & (angles < abins[nullXBounds[b]+1]) & (radii > rbins[nullYBounds[b]]) & (radii < rbins[nullYBounds[b]+1]) & (mass > (masspoint-(windowsize+windowOffset))) & (mass < (masspoint+(windowsize-windowOffset)))
        justTestMass = mass[abool]
        #if(len(justTestMass)>0):
        #    print('second bin with ', len(justTestMass), 'events in SR at ', rbins[nullYBounds[b]]*np.cos(abins[nullXBounds[b]]), ',', rbins[nullYBounds[b]]*np.sin(abins[nullXBounds[b]]))
        thebool.append( (angles > abins[nullXBounds[b]]) & (angles < abins[nullXBounds[b]+1]) & (radii > rbins[nullYBounds[b]]) & (radii < rbins[nullYBounds[b]+1]) )
    
    if(len(thebool) == 0):
        print('warning thebool len = 0')
        return 0,np.asarray([]),False, []
    
    thefinalbool = thebool[0]
    for b in range(1,len(nullXBounds)):
        thefinalbool = thefinalbool | thebool[b]

    masstoreturn = mass[thefinalbool]
    if(useWSbins):
        acceptMass = (masstoreturn > binDict[wsbin][1][0]) & (masstoreturn < binDict[wsbin][1][1])
    else:
        acceptMass = (masstoreturn > (masspoint-(windowsize+windowOffset))) & (masstoreturn < (masspoint+(windowsize-windowOffset)))
    massinwindow = masstoreturn[acceptMass]
    return len(massinwindow), masstoreturn, False, thefinalbool


class QUAKRunner:
    def __init__(self, options, massvalue):
        self.directory = options.directory
        self.template = options.template
        self.masspoint = massvalue
        self.useWSBins = options.useWSBins
        self.doBiasTest = options.doBiasTest
        self.doHighEta = options.doHighEta
        self.doBKGOnly = options.doBKGOnly
        self.doWeightingTest = options.doWeightingTest
        self.doJMETest = options.doJMETest

        self.reqEvs = options.reqEvs
        self.sigEvs = options.sigEvs
        self.evsInName = options.evsInName
        self.justPlots = options.justPlots
        self.allLimits = options.allLimits
        self.saveEffs = options.saveEffs
        self.saveSysts = options.saveSysts
        self.signalName = options.signalName
        self.samSignalName = options.samSignalName
        self.jsonKey = options.jsonKey
        
        self.binMin = options.binMin
        if(options.binMax<options.binMin):
            self.binMax = options.binMin + 100
        else:
            self.binMax = options.binMax
        self.sigOffset = options.sigOffset
        
        self.templateTop = os.getcwd()

        self.windowsize = 500
        self.windowsize = 300 #setting that got 3.3 sigma for w'
        #self.windowsize = 400
        self.histogrambins = 25
        self.verbose = options.verbose

        #self.tagname = 'BKG_Template_singleBin'
        self.tagname = 'BKG_Template_singleBin_'+str(self.masspoint)
        if(self.evsInName):
            self.tagname = 'BKG_Template_singleBin_req'+str(self.reqEvs)+'_'+str(self.masspoint)
        #self.tagname = 'fit_inputs_eff'+str(self.masspoint)

        self.WSbin = ''
        self.WSbinLOW = ''
        self.WSbinHIGH = ''
    
        self.templateHeader = 'case'
        #self.templateHeader = 'graviton'



    def makeH5s(self, jmevar=0):
        file1 = uproot.open('myOutput.root')
        tree1 = file1['output']
        
        branches = tree1.arrays()
        
        #mass = branches['mass']
        mass = branches['mjj']
        loss1 = branches['loss2'] #sig_loss
        loss2 = branches['loss1'] #bkg_loss
        label = branches['label']
        if('samp_label' in branches):
            samp_label = branches['samp_label']
        else:
            samp_label = np.zeros(len(branches['label']))
        if('mj1' in branches):
            mj1 = branches['mj1']
            mj2 = branches['mj2']
        else:
            mj1 = np.zeros(len(branches['label']))
            mj2 = np.zeros(len(branches['label']))
        if('evt_idx' in branches):
            evt_idx = branches['evt_idx']
        else:
            evt_idx = np.zeros(len(branches['label']))
        JME_flag = branches['JME_flag']
        
        loss1 = loss1.flatten()
        #loss1 = loss1 - self.sigOffset


        file2 = uproot.open('BKG.root')
        tree2 = file2['output']
        
        branches2 = tree2.arrays()
        
        mass_bkg = branches2['mjj']
        loss1_bkg = branches2['loss2'] #sig_loss
        loss1_bkg = loss1_bkg.flatten()
        loss2_bkg = branches2['loss1'] #bkg_loss
        label_bkg = branches2['label']
        JME_flag_bkg = np.asarray([-1]*len(mass_bkg))
        samp_label_bkg = np.zeros(len(mass_bkg))
        mj1_bkg = np.zeros(len(mass_bkg))
        mj2_bkg = np.zeros(len(mass_bkg))
        evt_idx_bkg = np.zeros(len(mass_bkg))
        

        mass = np.concatenate((mass, mass_bkg))
        loss1 = np.concatenate((loss1, loss1_bkg))
        loss2 = np.concatenate((loss2, loss2_bkg))
        label = np.concatenate((label, label_bkg))
        JME_flag = np.concatenate((JME_flag, JME_flag_bkg))
        samp_label = np.concatenate((samp_label, samp_label_bkg))
        mj1 = np.concatenate((mj1, mj1_bkg))
        mj2 = np.concatenate((mj2, mj2_bkg))
        evt_idx = np.concatenate((evt_idx, evt_idx_bkg))

        
        mass = mass[loss2>bkgCut]
        loss1 = loss1[loss2>bkgCut]
        label = label[loss2>bkgCut]
        samp_label = samp_label[loss2>bkgCut]
        mj1 = mj1[loss2>bkgCut]
        mj2 = mj2[loss2>bkgCut]
        evt_idx = evt_idx[loss2>bkgCut]
        JME_flag = JME_flag[loss2>bkgCut]
        loss2 = loss2[loss2>bkgCut]

        loss1 = loss1[mass>1460]
        label = label[mass>1460]
        samp_label = samp_label[mass>1460]
        mj1 = mj1[mass>1460]
        mj2 = mj2[mass>1460]
        evt_idx = evt_idx[mass>1460]
        JME_flag = JME_flag[mass>1460]
        loss2 = loss2[mass>1460]
        mass = mass[mass>1460]


        jmeselection = np.logical_or(JME_flag==-1, JME_flag==jmevar)
        mass = mass[jmeselection]
        loss1 = loss1[jmeselection]
        loss2 = loss2[jmeselection]
        label = label[jmeselection]
        samp_label = samp_label[jmeselection]
        mj1 = mj1[jmeselection]
        mj2 = mj2[jmeselection]
        evt_idx = evt_idx[jmeselection]
        JME_flag = JME_flag[jmeselection]

        print(mj1)
        print(mj2)

        
        print('the len of mass is ', len(mass), 'num sig = ', len(mass[label==1]), 'num bkg = ', len(mass[label==0]))

        """
        hf = h5py.File(self.tagname+"_allsignal.h5", 'w')
        hf.create_dataset('mjj', data=np.asarray(mass[(label==1)]))
        hf.close()
        hf = h5py.File(self.tagname+"_10percsignal_sigcut.h5", 'w')
        hf.create_dataset('mjj', data=np.asarray(mass[(label==1) & (loss1<-3.7)]))
        hf.close()
        hf = h5py.File(self.tagname+"_1percsignal_sigcut.h5", 'w')
        hf.create_dataset('mjj', data=np.asarray(mass[(label==1) & (loss1<-5)]))
        hf.close()
        hf = h5py.File(self.tagname+"_10percsignal_bkgcut.h5", 'w')
        hf.create_dataset('mjj', data=np.asarray(mass[(label==1) & (loss2>1.9)]))
        hf.close()
        hf = h5py.File(self.tagname+"_1percsignal_bkgcut.h5", 'w')
        hf.create_dataset('mjj', data=np.asarray(mass[(label==1) & (loss2>2.7)]))
        hf.close()
        """


        if(self.doBKGOnly):
            mass = mass[(label==0)]
            loss1 = loss1[(label==0)]
            loss2 = loss2[(label==0)]
            JME_flag = JME_flag[(label==0)]
            samp_label = samp_label[(label==0)]
            mj1 = mj1[(label==0)]
            mj2 = mj2[(label==0)]
            evt_idx = evt_idx[(label==0)]
            label = label[(label==0)]

            print('length of mass =', len(mass))

            """
            mj2copy = mj2
            mj1copy = mj1
            #squarebool = np.logical_and((mj1>280) & (mj1<340), (mj2>260) & (mj2<320))
            #squarebool = np.logical_and((mj1copy>250) & (mj1copy<360), (mj2copy>220) & (mj2copy<340))
            squarebool1 = np.logical_and((mj1>280) & (mj1<340), (mj2>260) & (mj2<320))
            squarebool2 = (mj1>160) & (mj1<200)
            squarebool = np.logical_or(squarebool1, squarebool2)
            squarebool = np.logical_not(squarebool)
            mass = mass[squarebool]
            loss1 = loss1[squarebool]
            loss2 = loss2[squarebool]
            JME_flag = JME_flag[squarebool]
            samp_label = samp_label[squarebool]
            label = label[squarebool]
            evt_idx = evt_idx[squarebool]
            mj1 = mj1[squarebool]
            mj2 = mj2[squarebool]
            """
            """
            mass = mass[((mj1<160) | (mj1>200))]
            loss1 = loss1[((mj1<160) | (mj1>200))]
            loss2 = loss2[((mj1<160) | (mj1>200))]
            JME_flag = JME_flag[((mj1<160) | (mj1>200))]
            samp_label = samp_label[((mj1<160) | (mj1>200))]
            mj2 = mj2[((mj1<160) | (mj1>200))]
            label = label[((mj1<160) | (mj1>200))]
            evt_idx = evt_idx[((mj1<160) | (mj1>200))]
            mj1 = mj1[((mj1<160) | (mj1>200))]
            
            
            mass = mass[((mj1<290) | (mj1>360)) & ((mj2<250) | (mj2>360))]
            loss1 = loss1[((mj1<290) | (mj1>360)) & ((mj2<250) | (mj2>360))]
            loss2 = loss2[((mj1<290) | (mj1>360)) & ((mj2<250) | (mj2>360))]
            JME_flag = JME_flag[((mj1<290) | (mj1>360)) & ((mj2<250) | (mj2>360))]
            samp_label = samp_label[((mj1<290) | (mj1>360)) & ((mj2<250) | (mj2>360))]
            label = label[((mj1<290) | (mj1>360)) & ((mj2<250) | (mj2>360))]
            evt_idx = evt_idx[((mj1<290) | (mj1>360)) & ((mj2<250) | (mj2>360))]
            mj2copy = mj2
            mj1copy = mj1
            mj1 = mj1[((mj1copy<290) | (mj1copy>360)) & ((mj2copy<250) | (mj2copy>360))]
            mj2 = mj2[((mj1copy<290) | (mj1copy>360)) & ((mj2copy<250) | (mj2copy>360))]
            #mj2 = mj2[((mj1<290) | (mj1>360))]
            #mj1 = mj1[((mj1<290) | (mj1>360))]
            #mj1 = mj1[((mj2<250) | (mj2>360))]
            #mj2 = mj2[((mj2<250) | (mj2>360))]
            """
            """
            mass = mass[(not((mj1>290) & (mj1<360) & (mj2>250) & (mj2<360)))]
            loss1 = loss1[(not((mj1>290) & (mj1<360) & (mj2>250) & (mj2<360)))]
            loss2 = loss2[(not((mj1>290) & (mj1<360) & (mj2>250) & (mj2<360)))]
            JME_flag = JME_flag[(not((mj1>290) & (mj1<360) & (mj2>250) & (mj2<360)))]
            samp_label = samp_label[(not((mj1>290) & (mj1<360) & (mj2>250) & (mj2<360)))]
            label = label[(not((mj1>290) & (mj1<360) & (mj2>250) & (mj2<360)))]
            evt_idx = evt_idx[(not((mj1>290) & (mj1<360) & (mj2>250) & (mj2<360)))]
            mj2copy = mj2
            mj1copy = mj1
            mj1 = mj1[(not((mj1copy>290) & (mj1copy<360) & (mj2copy>250) & (mj2copy<360)))]
            mj2 = mj2[(not((mj1copy>290) & (mj1copy<360) & (mj2copy>250) & (mj2copy<360)))]
            """

            """
            mass = mass[(samp_label==0)]
            loss1 = loss1[(samp_label==0)]
            loss2 = loss2[(samp_label==0)]
            JME_flag = JME_flag[(samp_label==0)]
            label = label[(samp_label==0)]
            samp_label = samp_label[(samp_label==0)]
            """
            """
            mass = np.concatenate((mass, mass[(samp_label==-2)], mass[(samp_label==-2)]))
            loss1 = np.concatenate((loss1, loss1[(samp_label==-2)], loss1[(samp_label==-2)]))
            loss2 = np.concatenate((loss2, loss2[(samp_label==-2)], loss2[(samp_label==-2)]))
            JME_flag = np.concatenate((JME_flag, JME_flag[(samp_label==-2)], JME_flag[(samp_label==-2)]))
            label = np.concatenate((label, label[(samp_label==-2)], label[(samp_label==-2)]))
            samp_label = np.concatenate((samp_label, samp_label[(samp_label==-2)], samp_label[(samp_label==-2)]))
            """
            """
            mass = np.concatenate((mass[(samp_label==-2)], mass[(samp_label==-2)]))
            loss1 = np.concatenate((loss1[(samp_label==-2)], loss1[(samp_label==-2)]))
            loss2 = np.concatenate((loss2[(samp_label==-2)], loss2[(samp_label==-2)]))
            JME_flag = np.concatenate((JME_flag[(samp_label==-2)], JME_flag[(samp_label==-2)]))
            label = np.concatenate((label[(samp_label==-2)], label[(samp_label==-2)]))
            samp_label = np.concatenate((samp_label[(samp_label==-2)], samp_label[(samp_label==-2)]))
            """
            """
            mass = mass[(samp_label==-2)]
            loss1 = loss1[(samp_label==-2)]
            loss2 = loss2[(samp_label==-2)]
            JME_flag = JME_flag[(samp_label==-2)]
            label = label[(samp_label==-2)]
            samp_label = samp_label[(samp_label==-2)]
            """
            """
            mass = np.concatenate((mass[(samp_label==-2)], mass[(samp_label==-2)], mass[(samp_label==-2)]))
            loss1 = np.concatenate((loss1[(samp_label==-2)], loss1[(samp_label==-2)], loss1[(samp_label==-2)]))
            loss2 = np.concatenate((loss2[(samp_label==-2)], loss2[(samp_label==-2)], loss2[(samp_label==-2)]))
            JME_flag = np.concatenate((JME_flag[(samp_label==-2)], JME_flag[(samp_label==-2)], JME_flag[(samp_label==-2)]))
            label = np.concatenate((label[(samp_label==-2)], label[(samp_label==-2)], label[(samp_label==-2)]))
            samp_label = np.concatenate((samp_label[(samp_label==-2)], samp_label[(samp_label==-2)], samp_label[(samp_label==-2)]))
            """

        print("sig loss average =", np.average(loss1))
        #self.sigOffset = np.average(loss1) - 9.5
        self.sigOffset = np.average(loss1) - 9
        #self.sigOffset = np.average(loss1) - 8.5
        loss1 = loss1 - self.sigOffset

        totalInWSBin = 0.
        minEventsInSR = 150.
        minEventsInSRv2 = 150.
        if(self.useWSBins):
            window_low = mass > binDict[self.WSbin][1][0]
            window_high = mass < binDict[self.WSbin][1][1]
            totalInWSBin = len(mass[window_low & window_high])
            totalInWSBinv2 = len(mass[(mass>(self.masspoint-300)) & (mass<(self.masspoint+300))])
            minEventsInSR = binDict[self.WSbin][2][0] * totalInWSBin
            minEventsInSRv2 = 0.01*totalInWSBinv2
            print('totalInWSBin =', totalInWSBin, 'minEventsInSR =', minEventsInSR, 'totalInWSBinv2 =', totalInWSBinv2, 'minEventsInSRv2 =', minEventsInSRv2)
        else:
            window_low = mass > (self.masspoint-(self.windowsize+windowOffset))
            window_high = mass < (self.masspoint+(self.windowsize-windowOffset))
        massWindow = window_low & window_high & (loss2 > bkgCut)
        for a in range(6):
            rbins = np.linspace(0, maxRad, radBins)
            abins = np.linspace(minangle,(0.4 + 0.02*a)*np.pi, angularbincount+a)
            loss2out = loss2[massWindow] - bkgCut
            loss1out = loss1[massWindow]
            angles = np.arctan(loss1out/loss2out)
            radii = np.sqrt(loss1out*loss1out + loss2out*loss2out)
            hist2, _, _ = np.histogram2d(list(angles), list(radii), bins=(abins, rbins))
            A, R = np.meshgrid(abins, rbins)
            plt.close()
            print(a, hist2.sum())
            if(hist2.sum()>5000):
                break
        binstoadd = a
        print('binstoadd = ',binstoadd)

        if(hist2.sum()<5000):
            #for so in [9., 8.5, 8., 7.5, 7., 6.5]:
            for so in [.5, 1., 1.5, 2., 2.5, 3.]:
                #self.sigOffset = self.sigOffset + so]
                self.sigOffset = self.sigOffset + 0.5
                print('np.average(loss1)', np.average(loss1))
                #loss1 = loss1 - so
                loss1 = loss1 - 0.5
                massWindow = window_low & window_high & (loss2 > bkgCut)
                rbins = np.linspace(0, maxRad, radBins)
                abins = np.linspace(minangle,(0.4 + 0.02*a)*np.pi, angularbincount+binstoadd)
                loss2out = loss2[massWindow] - bkgCut
                loss1out = loss1[massWindow]
                angles = np.arctan(loss1out/loss2out)
                radii = np.sqrt(loss1out*loss1out + loss2out*loss2out)
                hist2, _, _ = np.histogram2d(list(angles), list(radii), bins=(abins, rbins))
                A, R = np.meshgrid(abins, rbins)
                plt.close()
                print(so, hist2.sum())
                if(hist2.sum()>5000):
                    break

        massMin = mass > 1000
        #testWindow1 = mass < (self.masspoint-(self.windowsize+windowOffset))
        #testWindow2 = mass > (self.masspoint+(self.windowsize-windowOffset))
        #massWindow = massMin & (testWindow1 | testWindow2) & (loss2 > bkgCut)
        

        rbins = np.linspace(0, maxRad, radBins)
        abins = np.linspace(minangle,(0.4 + 0.02*binstoadd)*np.pi, angularbincount+binstoadd)

        massWindow_testIn = (mass>(self.masspoint-(self.windowsize+windowOffset))) & (mass<(self.masspoint+(self.windowsize-windowOffset)))
        massWindow_testIn_SIG = (mass>(self.masspoint-(self.windowsize+windowOffset))) & (mass<(self.masspoint+(self.windowsize-windowOffset))) & (label==1)
        loss2_testIn = loss2[massWindow_testIn] - bkgCut
        loss1_testIn = loss1[massWindow_testIn]
        loss2_testIn_SIG = loss2[massWindow_testIn_SIG] - bkgCut
        loss1_testIn_SIG = loss1[massWindow_testIn_SIG]

        print('number of events in window =', len(loss1_testIn), 'number of signal events in window =', len(loss1_testIn_SIG))
        
        angles_testIn = np.arctan(loss1_testIn/loss2_testIn)
        radii_testIn = np.sqrt(loss1_testIn*loss1_testIn + loss2_testIn*loss2_testIn)
        angles_testIn_SIG = np.arctan(loss1_testIn_SIG/loss2_testIn_SIG)
        radii_testIn_SIG = np.sqrt(loss1_testIn_SIG*loss1_testIn_SIG + loss2_testIn_SIG*loss2_testIn_SIG)

        hist_testIn, _, _ = np.histogram2d(list(angles_testIn), list(radii_testIn), bins=(abins, rbins))
        hist_testIn_SIG, _, _ = np.histogram2d(list(angles_testIn_SIG), list(radii_testIn_SIG), bins=(abins, rbins))

        print('number of events in histo =', np.sum(hist_testIn), 'number of signal events in histo =', np.sum(hist_testIn_SIG))

        numAll = 0.
        numSIG = 0.
        numAll_ycut = 0.
        numSIG_ycut = 0.
        for a in range(angularbincount+binstoadd-1):
            for r in range(radBins-1):
                if(hist_testIn[a,r]>0 and hist_testIn[a,r]<40):
                    print('bin angle =', abins[a], ', bin r =', rbins[r], 'x =', rbins[r]*np.cos(abins[a]), 'y =', rbins[r]*np.sin(abins[a]), 'regular pop =', hist_testIn[a,r], 'sig only =', hist_testIn_SIG[a,r], 'frac =', hist_testIn_SIG[a,r]/hist_testIn[a,r])
                    numAll+=hist_testIn[a,r]
                    numSIG+=hist_testIn_SIG[a,r]
                    if(rbins[r]*np.sin(abins[a])<5.7):
                        numAll_ycut+=hist_testIn[a,r]
                        numSIG_ycut+=hist_testIn_SIG[a,r]
            print('end of wedge '+str(abins[a])+'. numAll = ', numAll, 'numSIG =', numSIG, 'numAll_ycut =', numAll_ycut, 'numSIG_ycut =', numSIG_ycut)
                #elif(hist_testIn_SIG[a,r]>0):
                #    print('bin angle =', abins[a], ', bin r =', rbins[r], 'regular pop =', hist_testIn[a,r], 'sig only =', hist_testIn_SIG[a,r], 'PURE')
                #else:
                #    print('bin angle =', abins[a], ', bin r =', rbins[r], 'regular pop =', hist_testIn[a,r], 'sig only =', hist_testIn_SIG[a,r], 'nothing')










        #sidebandsize = 500
        if(self.useWSBins):
            testWindow1 = (mass > binDict[self.WSbinLOW][1][0]) & (mass < binDict[self.WSbinLOW][1][1])
            testWindow2 = (mass > binDict[self.WSbinHIGH][1][0]) & (mass < binDict[self.WSbinHIGH][1][1])
        else:
            testWindow1 = (mass > (self.masspoint-(self.windowsize+windowOffset+sidebandsize))) & (mass < (self.masspoint-(self.windowsize+windowOffset)))
            testWindow2 = (mass > (self.masspoint+(self.windowsize-windowOffset))) & (mass < (self.masspoint+(self.windowsize+windowOffset+sidebandsize)))
        #massWindow = massMin & (testWindow1 | testWindow2)
        massWindow = massMin & (testWindow1 | testWindow2) & (loss2 > bkgCut)
        print('total sideband events =', len(mass[massWindow]))

        rbins = np.linspace(0, maxRad, radBins)
        abins = np.linspace(minangle,(0.4 + 0.02*binstoadd)*np.pi, angularbincount+binstoadd)

        #calculate histogram
        loss2out = loss2[massWindow] - bkgCut
        loss1out = loss1[massWindow]
        
        angles = np.arctan(loss1out/loss2out)
        radii = np.sqrt(loss1out*loss1out + loss2out*loss2out)
        
        hist, _, _ = np.histogram2d(list(angles), list(radii), bins=(abins, rbins))
        A, R = np.meshgrid(abins, rbins)
        
        plt.close()
        
        totalmassinwindow = 0
        themassarray = []
        thefinalfinalbool = np.zeros(len(mass), dtype=bool)


        print('the len of mass is ', len(mass), 'num sig = ', len(mass[label==1]), 'num bkg = ', len(mass[label==0]))
        print('boring mass between 2900 and 3100 ', len(mass[(mass>2900) & (mass<3100)]))
        themax = hist.max()
        for b in range(1000):
            print('b=',b)
            sliced = massArraySlicer(hist, b, mass[loss2>bkgCut], loss1[loss2>bkgCut], loss2[loss2>bkgCut]-bkgCut, self.masspoint, themax, binstoadd, self.useWSBins, self.WSbin, self.windowsize)
            totalmassinwindow += sliced[0]
            print('len(thefinalfinalbool) =', len(thefinalfinalbool))
            print('len(sliced[3]) = ', len(sliced[3]))
            if(len(thefinalfinalbool) == len(sliced[3])):
                thefinalfinalbool = np.logical_or(thefinalfinalbool, sliced[3])
            #print(thefinalfinalbool.sum())
            themassarray = themassarray + list(sliced[1])
            if(self.useWSBins):
                #if(totalmassinwindow>minEventsInSR or sliced[2]):
                #if(totalmassinwindow>minEventsInSRv2 or sliced[2]):
                if(totalmassinwindow>max(350, 100.*np.power((self.masspoint/1000.)-3, 4)) or sliced[2]):
                    break
            if(self.doHighEta):
                #print(totalmassinwindow, sliced[2], 'num sig = ', len(mass[np.logical_and(thefinalfinalbool, label==1)]), 'num ttbar =', len(mass[np.logical_and(thefinalfinalbool, samp_label==-2)]), 'num single top =', len(mass[np.logical_and(thefinalfinalbool, samp_label==-1)]))
                #if(totalmassinwindow>max(300, 100.*np.power((self.masspoint/1000.)-3, 4)) or sliced[2]): #where I get the 3.3 sigma
                if(totalmassinwindow>max(self.reqEvs, 100.*np.power((self.masspoint/1000.)-3, 4)) or sliced[2]):
                    break
            else:
                print(totalmassinwindow, sliced[2], 'num sig = ', len(mass[np.logical_and(thefinalfinalbool, label==1)]))
                #if(totalmassinwindow>max(300, 100.*np.power((self.masspoint/1000.)-3, 4)) or sliced[2]):
                #if(totalmassinwindow>max(150, 100.*np.power((self.masspoint/1000.)-3, 4)) or sliced[2]):
                #if(totalmassinwindow>max(250, 100.*np.power((self.masspoint/1000.)-3, 4)) or sliced[2]):

                #this is what I used for March 13 scans
                #if(totalmassinwindow>max(200, 100.*np.power((self.masspoint/1000.)-3, 4)) or sliced[2]):

                #For first UNBLINDED_DATA run
                #if(totalmassinwindow>max(self.reqEvs, 100.*np.power((self.masspoint/1000.)-3, 4)) or sliced[2]):
                
                #This is what I used for data events...
                #if(totalmassinwindow>max(self.reqEvs, 100.*np.power((self.masspoint/1000.)-3, 3)) or sliced[2]):
                
                #This is what I used for data events August 15...
                if(totalmassinwindow>max(self.reqEvs, (self.reqEvs/1500.)*500.*np.power((self.masspoint/1000.)-3, 3)) or sliced[2]):
                

                #if(totalmassinwindow>500 or sliced[2]):

                #testing apr19 scans
                #if(totalmassinwindow>max(300, 100.*np.power((self.masspoint/1000.)-3, 4)) or sliced[2]):

                #if(totalmassinwindow>max(250, 700.*np.power((self.masspoint/1000.)-3, 1)) or sliced[2]):
                #if(totalmassinwindow>max(250, 1000.*np.power((self.masspoint/1000.)-3, 1)) or sliced[2]):
                    break


        themassarray = np.asarray(themassarray)
        massForWeighting = mass[np.logical_and(thefinalfinalbool, label==1)]
        #massForQCD = mass[np.logical_and(thefinalfinalbool, samp_label==0)]
        #massForTTBAR = mass[np.logical_and(thefinalfinalbool, samp_label==-2)]
        #massForST = mass[np.logical_and(thefinalfinalbool, samp_label==-1)]
        if(self.useWSBins):
            window_low2 = themassarray > binDict[self.WSbin][1][0]
            window_high2 = themassarray < binDict[self.WSbin][1][1]
            window_low3 = massForWeighting > binDict[self.WSbin][1][0]
            window_high3 = massForWeighting < binDict[self.WSbin][1][1]
        else:
            window_low2 = themassarray > (self.masspoint-(self.windowsize+windowOffset))
            window_high2 = themassarray < (self.masspoint+(self.windowsize-windowOffset))
            window_low3 = massForWeighting > (self.masspoint-(self.windowsize+windowOffset))
            window_high3 = massForWeighting < (self.masspoint+(self.windowsize-windowOffset))
        print('we will look at', len(themassarray), 'events.', len(massForWeighting), 'of those are signal events.', len(themassarray[window_low2 & window_high2]), 'events are in window, and', len(massForWeighting[window_low3 & window_high3]), 'of those are signal events')
        window2900 = themassarray > 2900
        window3100 = themassarray > 3100
        print('simple select 2900 to 3100: ', len(themassarray[window2900 & window2900]))
        print('in window, we have QCD:', len(mass[(thefinalfinalbool) & (samp_label==0) & (mass > (self.masspoint-(self.windowsize+windowOffset))) & (mass < (self.masspoint+(self.windowsize-windowOffset)))]), 'ttbar:', len(mass[(thefinalfinalbool) & (samp_label==-2) & (mass > (self.masspoint-(self.windowsize+windowOffset))) & (mass < (self.masspoint+(self.windowsize-windowOffset)))]), 'single top:', len(mass[(thefinalfinalbool) & (samp_label==-1) & (mass > (self.masspoint-(self.windowsize+windowOffset))) & (mass < (self.masspoint+(self.windowsize-windowOffset)))]), 'V+jets:', len(mass[(thefinalfinalbool) & (samp_label==-3) & (mass > (self.masspoint-(self.windowsize+windowOffset))) & (mass < (self.masspoint+(self.windowsize-windowOffset)))]))
        if(len(massForWeighting)>0):
            mainfit = stats.norm.fit(massForWeighting)
            print('mainfit = ', mainfit)

        if(self.saveEffs):
            effDict = {self.jsonKey: {"samSignalName": self.samSignalName, "selectionEfficiency": float(len(massForWeighting))/float(self.sigEvs)}}
            with open(self.samSignalName+"_selectionEfficiency.json", "w") as outfile:
                json.dump(effDict, outfile)
            
        if(self.doWeightingTest):

            systsvec = [0.15]
            
            #print('main crystalball = ', stats.crystalball.fit(massForWeighting))
            mass2 = branches['mjj']
            bkgloss_forCut = branches['loss1']
            JMEflag_forCut = branches['JME_flag']

            #loss2_bkg_forCut = branches2['loss1'] #bkg_loss
            #JME_flag_bkg_forCut = np.asarray([-1]*len(loss2_bkg_forCut))
            #bkgloss_forCut = np.concatenate((bkgloss_forCut, loss2_bkg_forCut))
            #JMEflag_forCut = np.concatenate((JMEflag_forCut, JME_flag_bkg_forCut))

            jmeselection2 = np.logical_and(mass2>1460, JMEflag_forCut==jmevar)
            bkgloss_forCut = bkgloss_forCut[jmeselection2]
            cutBool = bkgloss_forCut > bkgCut
            #selectionBool = np.logical_and(thefinalfinalbool, label==1)
            print('cutBool.sum() = ', cutBool.sum())
            #print('selectionBool size = ', len(selectionBool))
            lund_pt_var = branches['lund_pt_var']
            lund_pt_var = lund_pt_var[jmeselection2]
            lund_stat_var = branches['lund_stat_var']
            lund_stat_var = lund_stat_var[jmeselection2]
            lund_nom = branches['lund_nom']
            lund_nom = lund_nom[jmeselection2]
            #print('lund_nom.shape', lund_nom.shape, lund_nom)
            #print('lund_pt_var.shape', lund_pt_var.shape, lund_pt_var)

            lund_nom = branches['lund_nom']
            lund_nom = lund_nom[jmeselection2]
            lund_nom = lund_nom[cutBool]
            lund_pt_var = lund_pt_var[cutBool]
            lund_stat_var = lund_stat_var[cutBool]
            bkg_selection_bool = np.logical_and(loss2_bkg>bkgCut, mass_bkg>1460)
            print('len(lund_nom)', len(lund_nom), 'len(mass_bkg[bkg_selection_bool])', len(mass_bkg[bkg_selection_bool]), 'add =', len(lund_nom)+len(mass_bkg[bkg_selection_bool]))
            bkg_weights = np.asarray([1]*len(mass_bkg[bkg_selection_bool]))
            lund_nom = np.concatenate((lund_nom, bkg_weights))
            bkg_selection_bool = np.logical_and(loss2_bkg>bkgCut, mass_bkg>1460)
            bkg_weights_lund_pt_stat = np.asarray([[1]*100]*len(mass_bkg[bkg_selection_bool]))
            lund_pt_var = np.concatenate((lund_pt_var, bkg_weights_lund_pt_stat))
            lund_stat_var = np.concatenate((lund_stat_var, bkg_weights_lund_pt_stat))
            selectionBool = np.logical_and(thefinalfinalbool, label==1)
            lund_pt_var = lund_pt_var[selectionBool]
            lund_stat_var = lund_stat_var[selectionBool]
            lund_nom = lund_nom[selectionBool]

            ptvararray = []
            for i in range(100):
                ptvararray.append(sum(lund_pt_var[:,i]/lund_nom))
            ptvararray = np.asarray(ptvararray)
            statvararray = []
            for i in range(100):
                statvararray.append(sum(lund_stat_var[:,i]/lund_nom))
            statvararray = np.asarray(statvararray)

            print('ptvararray, avg:', np.average(ptvararray), ', std:', np.std(ptvararray))
            systsvec.append( np.average(ptvararray)/float(self.sigEvs) - float(len(massForWeighting))/float(self.sigEvs) + np.std(ptvararray)/float(self.sigEvs))
            print('statvararray, avg:', np.average(statvararray), ', std:', np.std(statvararray))
            systsvec.append( np.average(statvararray)/float(self.sigEvs) - float(len(massForWeighting))/float(self.sigEvs) + np.std(statvararray)/float(self.sigEvs))
            
            lund_matching = branches['lund_matching']
            print('lund_matching=', lund_matching[0])
            systsvec.append( lund_matching[0] * float(len(massForWeighting))/float(self.sigEvs) )

            updownsyst = []
            for w in theSysWeights:
                print(w)
                if('up' in w):
                    updownsyst = []
                lund_nom = branches['lund_nom']
                #lund_pt_var = branches['lund_pt_var']
                #lund_pt_var = lund_pt_var[JMEflag_forCut==jmevar]
                #lund_stat_var = branches['lund_stat_var']
                #lund_stat_var = lund_stat_var[JMEflag_forCut==jmevar]
                weights = branches[w]
                weights = weights[jmeselection2]
                lund_nom = lund_nom[jmeselection2]
                weights = weights[cutBool]
                lund_nom = lund_nom[cutBool]
                #lund_pt_var = lund_pt_var[cutBool]
                #lund_stat_var = lund_stat_var[cutBool]
                #print('weights', weights, np.average(weights))
                #print('lund_nom', lund_nom, np.average(lund_nom))
                ##weights = weights*lund_nom
                #print('weights again', weights, np.average(weights))
                bkg_selection_bool = np.logical_and(loss2_bkg>bkgCut, mass_bkg>1460)
                bkg_weights = np.asarray([1]*len(mass_bkg[bkg_selection_bool]))
                weights = np.concatenate((weights, bkg_weights))
                selectionBool = np.logical_and(thefinalfinalbool, label==1)
                weights = weights[selectionBool]

                print('weights final avg', np.average(weights), 'median', np.median(weights))
                print('len(weights) = ', len(weights))
                print('weights.sum() = ', weights.sum())
                #print(weights)
                #massplot = plt.hist(massForWeighting, bins=50, weights=weights)
                #normfit = stats.norm.fit(massplot)
                #print(normfit)
                sums, bins = np.histogram(massForWeighting, bins=50, weights=weights)
                hist_dist = stats.rv_histogram((sums, bins))
                weighted_data = hist_dist.rvs(size=100000)
                mean, stdev = stats.norm.fit(weighted_data)
                print(mean, stdev)
                updownsyst.append( float(abs(weights.sum() - len(weights))) / float(len(weights)) )
                
                if('down' in w):
                    """
                    maxval = 0.
                    for sys in updownsyst:
                        if(sys>maxval):
                            maxval = sys
                    systsvec.append(float(maxval))
                    """
                    systsvec.append(np.average(updownsyst))
                #print(stats.crystalball.fit(weighted_data))
                #stats.crystalball.fit(massForWeighting)

            if(self.saveSysts):
                totalsyst = 0.
                for sys in systsvec:
                    totalsyst += sys*sys
                totalsyst = np.sqrt(totalsyst)
                sysDict = {self.jsonKey: {"samSignalName": self.samSignalName, "signalSystematics": systsvec, "fullSystematic": totalsyst}}
                with open(self.samSignalName+"_signalSystematic.json", "w") as outfile:
                    json.dump(sysDict, outfile)
        #specialcut = (loss2 > bkgCut) & ((loss1)<1.8*(loss2-bkgCut)) & (np.sqrt((loss2-bkgCut)*(loss2-bkgCut) + (loss1)*(loss1)) > 5.) & (np.sqrt((loss2-bkgCut)*(loss2-bkgCut) + (loss1)*(loss1)) < 7.)
        #themassarray = mass[specialcut]
        #specialcut2 = (loss2 > bkgCut) & ((loss1)<1.8*(loss2-bkgCut)) & (np.sqrt((loss2-bkgCut)*(loss2-bkgCut) + (loss1)*(loss1)) > 5.) & (np.sqrt((loss2-bkgCut)*(loss2-bkgCut) + (loss1)*(loss1)) < 7.) & (label==1)
        

        """
        mj1plot = plt.hist(list(mj1[(thefinalfinalbool)]), bins=100)
        plt.title('mj1 spectrum for all selected events')
        plt.xlabel('mj1')
        plt.ylabel('count')
        plt.savefig("Mar27_all_mj1_spectrum.png")
        plt.close()
        mj2plot = plt.hist(list(mj2[(thefinalfinalbool)]), bins=100)
        plt.title('mj2 spectrum for all selected events')
        plt.xlabel('mj2')
        plt.ylabel('count')
        plt.savefig("Mar27_all_mj2_spectrum.png")
        plt.close()
        
        if('mj1' in branches):
            squarebool = np.logical_and((mj1>280) & (mj1<340), (mj2>260) & (mj2<320))
            #squarebool = np.logical_and(squarebool, (samp_label==-2))
            #squarebool = np.logical_and((mj1>180) & (mj1<240), (mj2>160) & (mj2<220))
            #squarebool = (mj1>160) & (mj1<200)

            #squarebool1 = np.logical_and((mj1>280) & (mj1<340), (mj2>260) & (mj2<320))
            #squarebool2 = (mj1>160) & (mj1<200)
            #squarebool = np.logical_or(squarebool1, squarebool2)
            #squarebool = np.logical_and((mj1>250) & (mj1<360), (mj2>220) & (mj2<340))
            thecutmassarray = mass[(thefinalfinalbool) & (squarebool) & (mass>2000)]
            binsx = [2037, 2132, 2231, 2332, 2438,
                 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854,
                 4010, 4171, 4337, 4509, 4700, 4900,  5100, 5300, 5500, 5800,
                 6100, 6400, 6800]
            realbins = []
            for b in binsx:
                realbins.append(b)
                if(b>max(thecutmassarray)):
                    break
            ttbarmassplot = plt.hist(list(thecutmassarray), bins=realbins)
            plt.title('all mass spectrum in high mj box')
            plt.xlabel('mjj')
            plt.ylabel('count')
            plt.savefig("Mar27_mjj_spectrum_inHighMJBox.png")
            plt.close()


            mj1plot = plt.hist(list(mj1[(thefinalfinalbool) & (mass>2000)]), bins=range(0,501,10))
            plt.title('mj1 spectrum for all selected events')
            plt.xlabel('mj1')
            plt.ylabel('count')
            plt.xlim([0,500])
            plt.savefig("Mar27_all_mj1_spectrum_m2000.png")
            plt.close()
            mj2plot = plt.hist(list(mj2[(thefinalfinalbool) & (mass>2000)]), bins=range(0,501,10))
            plt.title('mj2 spectrum for all selected events')
            plt.xlabel('mj2')
            plt.ylabel('count')
            plt.xlim([0,500])
            plt.savefig("Mar27_all_mj2_spectrum_m2000.png")
            plt.close()

            mj1plot = plt.hist2d(list(mj1[(thefinalfinalbool)]), list(mj2[(thefinalfinalbool)]), bins=100)
            plt.title('mj1 and mj2 spectrum for all selected events')
            plt.xlabel('mj1')
            plt.ylabel('mj2')
            plt.savefig("Mar27_all_mj1_and_mj2_spectrum.png")
            plt.close()

            mj1scatter = plt.scatter(list(mj1[(thefinalfinalbool) & (mass>2000)]), list(mj2[(thefinalfinalbool) & (mass>2000)]), c=list(mass[(thefinalfinalbool) & (mass>2000)]))
            plt.title('mj1 and mj2 spectrum for all selected events')
            plt.xlabel('mj1')
            plt.ylabel('mj2')
            plt.colorbar()
            plt.xlim([0,500])
            plt.ylim([0,500])
            plt.savefig("Mar27_all_mj1_and_mj2_scatter.png")
            plt.close()


            for i in range(10):
                bottomcut = 1900 + i*200
                topcut = 2100+i*200
                fig, ax = plt.subplots()
                ax.set_aspect("equal")
                #hist, xbins, ybins, im = ax.hist2d(x,y, bins=range(6))
                hist, xbins, ybins, im = ax.hist2d(list(mj1[(thefinalfinalbool) & (mass>bottomcut) & (mass<topcut)]), list(mj2[(thefinalfinalbool) & (mass>bottomcut) & (mass<topcut)]), bins=[10,10], range=[[0, 500], [0, 500]])
                for i in range(len(ybins)-1):
                    for j in range(len(xbins)-1):
                        ax.text(xbins[j]+25,ybins[i]+25, "{:.2f}".format(hist.T[i,j]/sum(sum(hist))), color="w", ha="center", va="center", fontweight="bold")
                #fillhist2d = plt.hist2d(list(mj1[(thefinalfinalbool) & (mass>bottomcut) & (mass<topcut)]), list(mj2[(thefinalfinalbool) & (mass>bottomcut) & (mass<topcut)]), bins=[10,10], range=[[0, 500], [0, 500]], normed=True)
                #fillhist2d[0] = np.divide(fillhist2d[0], sum(sum(fillhist2d[0])))
                plt.title('mj1 and mj2 spectrum for events with mjj between '+str(bottomcut)+' and '+str(topcut))
                plt.xlabel('mj1')
                plt.ylabel('mj2')
                plt.xlim([0,500])
                plt.ylim([0,500])
                plt.savefig("Mar27_mj1_and_mj2_hist"+str(bottomcut)+"_"+str(topcut)+".png")
                plt.close()


            lowlowcut = 2100
            lowtopcut = 2400
            mj1scatter = plt.scatter(list(mj1[(thefinalfinalbool) & (mass>lowlowcut) & (mass<lowtopcut)]), list(mj2[(thefinalfinalbool) & (mass>lowlowcut) & (mass<lowtopcut)]), c=list(mass[(thefinalfinalbool) & (mass>lowlowcut) & (mass<lowtopcut)]))
            numlow = len(mj1[(thefinalfinalbool) & (mass>lowlowcut) & (mass<lowtopcut)])
            numlowbox = len(mj1[(thefinalfinalbool) & (mass>lowlowcut) & (mass<lowtopcut) & (squarebool)])
            print('in the low mjj region, there are ', numlow, 'events. ', numlowbox, 'of these are in the box, which is ', float(numlowbox)/float(numlow))
            plt.title('mj1 and mj2 spectrum for low sideband events')
            plt.xlabel('mj1')
            plt.ylabel('mj2')
            plt.colorbar()
            plt.xlim([0,500])
            plt.ylim([0,500])
            plt.savefig("Mar27_low_mj1_and_mj2_scatter.png")
            plt.close()

            midlowcut = self.masspoint - 200
            midtopcut = self.masspoint + 200
            mj1scatter = plt.scatter(list(mj1[(thefinalfinalbool) & (mass>midlowcut) & (mass<midtopcut)]), list(mj2[(thefinalfinalbool) & (mass>midlowcut) & (mass<midtopcut)]), c=list(mass[(thefinalfinalbool) & (mass>midlowcut) & (mass<midtopcut)]))
            numlow = len(mj1[(thefinalfinalbool) & (mass>midlowcut) & (mass<midtopcut)])
            numlowbox = len(mj1[(thefinalfinalbool) & (mass>midlowcut) & (mass<midtopcut) & (squarebool)])
            print('in the SR mjj region, there are ', numlow, 'events. ', numlowbox, 'of these are in the box, which is ', float(numlowbox)/float(numlow))
            plt.title('mj1 and mj2 spectrum for SR events')
            plt.xlabel('mj1')
            plt.ylabel('mj2')
            plt.colorbar()
            plt.xlim([0,500])
            plt.ylim([0,500])
            plt.savefig("Mar27_SR_mj1_and_mj2_scatter.png")
            plt.close()

            highmin = self.masspoint + 300
            mj1scatter = plt.scatter(list(mj1[(thefinalfinalbool) & (mass>highmin)]), list(mj2[(thefinalfinalbool) & (mass>highmin)]), c=list(mass[(thefinalfinalbool) & (mass>highmin)]))
            numlow = len(mj1[(thefinalfinalbool) & (mass>highmin)])
            numlowbox = len(mj1[(thefinalfinalbool) & (mass>highmin) & (squarebool)])
            print('in the high mjj region, there are ', numlow, 'events. ', numlowbox, 'of these are in the box, which is ', float(numlowbox)/float(numlow))
            plt.title('mj1 and mj2 spectrum for high sideband events')
            plt.xlabel('mj1')
            plt.ylabel('mj2')
            plt.colorbar()
            plt.xlim([0,500])
            plt.ylim([0,500])
            plt.savefig("Mar27_high_mj1_and_mj2_scatter.png")
            plt.close()


            mj1plot = plt.hist(list(mj1[(thefinalfinalbool) & (mass>2500) & (mass<2900)]), bins=100)
            plt.title('mj1 spectrum for SR events')
            plt.xlabel('mj1')
            plt.ylabel('count')
            plt.savefig("Mar27_SR_mj1_spectrum.png")
            plt.close()
            mj2plot = plt.hist(list(mj2[(thefinalfinalbool) & (mass>2500) & (mass<2900)]), bins=100)
            plt.title('mj2 spectrum for SR events')
            plt.xlabel('mj2')
            plt.ylabel('count')
            plt.savefig("Mar27_SR_mj2_spectrum.png")
            plt.close()
            mj2plot = plt.hist(list(mj2[(thefinalfinalbool) & (mass>2500) & (mass<2900) & (mj1>290)]), bins=100)
            plt.title('mj2 spectrum for SR events with mj1>300')
            plt.xlabel('mj2')
            plt.ylabel('count')
            plt.savefig("Mar27_SR_mj2_spectrum_highmj1.png")
            plt.close()
            mj2plot = plt.hist(list(mj2[(thefinalfinalbool) & (mass>2500) & (mass<2900) & (mj1>160) & (mj1<200)]), bins=100)
            plt.title('mj2 spectrum for SR events with mj1 near top')
            plt.xlabel('mj2')
            plt.ylabel('count')
            plt.savefig("Mar27_SR_mj2_spectrum_neartop.png")
            plt.close()
            mj1plot = plt.hist2d(list(mj1[(thefinalfinalbool) & (mass>2500) & (mass<2900)]), list(mj2[(thefinalfinalbool) & (mass>2500) & (mass<2900)]), bins=100)
            plt.title('mj1 and mj2 spectrum for SR events')
            plt.xlabel('mj1')
            plt.ylabel('mj2')
            plt.savefig("Mar27_SR_mj1_and_mj2_spectrum.png")
            plt.close()


            mj1plot = plt.hist(list(mj1[(thefinalfinalbool) & (((mass>2000) & (mass<2300)) | ((mass>3000) & (mass<3500)))]), bins=100)
            plt.title('mj1 spectrum for CR events')
            plt.xlabel('mj1')
            plt.ylabel('count')
            plt.savefig("Mar27_CR_mj1_spectrum.png")
            plt.close()
            mj2plot = plt.hist(list(mj2[(thefinalfinalbool) & (((mass>2000) & (mass<2300)) | ((mass>3000) & (mass<3500)))]), bins=100)
            plt.title('mj2 spectrum for CR events')
            plt.xlabel('mj2')
            plt.ylabel('count')
            plt.savefig("Mar27_CR_mj2_spectrum.png")
            plt.close()

            thecutmassarray = mass[(thefinalfinalbool)]
            binsx = [2037, 2132, 2231, 2332, 2438,
                 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854,
                 4010, 4171, 4337, 4509, 4700, 4900,  5100, 5300, 5500, 5800,
                 6100, 6400, 6800]
            realbins = []
            for b in binsx:
                realbins.append(b)
                if(b>max(thecutmassarray)):
                    break
            ttbarmassplot = plt.hist(list(thecutmassarray), bins=realbins)
            plt.title('all mass spectrum')
            plt.xlabel('mjj')
            plt.ylabel('count')
            plt.savefig("Mar27_all_mjj_spectrum.png")
            plt.close()

            thecutmassarray = mass[(thefinalfinalbool) & ((mj1<160) | (mj1>200))]
            binsx = [2037, 2132, 2231, 2332, 2438,
                 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854,
                 4010, 4171, 4337, 4509, 4700, 4900,  5100, 5300, 5500, 5800,
                 6100, 6400, 6800]
            realbins = []
            for b in binsx:
                realbins.append(b)
                if(b>max(thecutmassarray)):
                    break
            ttbarmassplot = plt.hist(list(thecutmassarray), bins=realbins)
            plt.title('masked mass spectrum')
            plt.xlabel('mjj')
            plt.ylabel('count')
            plt.savefig("Mar27_masked_mjj_spectrum.png")
            plt.close()


            thecutmassarray = mass[(thefinalfinalbool) & ((mj1<290) | (mj1>360)) & ((mj2<250) | (mj2>360))]
            binsx = [2037, 2132, 2231, 2332, 2438,
                 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854,
                 4010, 4171, 4337, 4509, 4700, 4900,  5100, 5300, 5500, 5800,
                 6100, 6400, 6800]
            realbins = []
            for b in binsx:
                realbins.append(b)
                if(b>max(thecutmassarray)):
                    break
            ttbarmassplot = plt.hist(list(thecutmassarray), bins=realbins)
            plt.title('masked2 mass spectrum')
            plt.xlabel('mjj')
            plt.ylabel('count')
            plt.savefig("Mar27_masked2_mjj_spectrum.png")
            plt.close()

            if(len(mass[(thefinalfinalbool) & (samp_label==-2)]) > 0):
                thecutmassarray = mass[(thefinalfinalbool) & (samp_label==-2)]
                binsx = [2037, 2132, 2231, 2332, 2438,
                     2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854,
                     4010, 4171, 4337, 4509, 4700, 4900,  5100, 5300, 5500, 5800,
                     6100, 6400, 6800]
                realbins = []
                for b in binsx:
                    realbins.append(b)
                    if(b>max(thecutmassarray)):
                        break
                ttbarmassplot = plt.hist(list(thecutmassarray), bins=realbins)
                plt.title('ttbar mass spectrum')
                plt.xlabel('mjj')
                plt.ylabel('count')
                plt.savefig("Mar27_ttbar_spectrum.png")
                plt.close()


                mj1plot = plt.hist(list(mj1[(thefinalfinalbool) & (samp_label!=-2)]), bins=100)
                plt.title('mj1 spectrum for non ttbar events')
                plt.xlabel('mj1')
                plt.ylabel('count')
                plt.savefig("Mar27_other_mj1_spectrum.png")
                plt.close()
                mj2plot = plt.hist(list(mj2[(thefinalfinalbool) & (samp_label!=-2)]), bins=100)
                plt.title('mj2 spectrum for non ttbar events')
                plt.xlabel('mj2')
                plt.ylabel('count')
                plt.savefig("Mar27_other_mj2_spectrum.png")
                plt.close()
                mj1plot = plt.hist(list(mj1[(thefinalfinalbool) & (samp_label==-2) & (mass>2000)]), bins=100)
                plt.title('mj1 spectrum for ttbar events')
                plt.xlabel('mj1')
                plt.ylabel('count')
                plt.savefig("Mar27_ttbar_mj1_spectrum.png")
                plt.close()
                mj2plot = plt.hist(list(mj2[(thefinalfinalbool) & (samp_label==-2) & (mass>2000)]), bins=100)
                plt.title('mj2 spectrum for ttbar events')
                plt.xlabel('mj2')
                plt.ylabel('count')
                plt.savefig("Mar27_ttbar_mj2_spectrum.png")
                plt.close()
        """
        
        print('going to fit', len(themassarray), 'events')
        #print('there are', len(mass[specialcut2]), 'signal events fyi')
        #tagname = 'BKG_Template_singleBin'

        #squarebool1 = np.logical_and((mj1>270) & (mj1<350), (mj2>250) & (mj2<330))
        #squarebool2 = (mj1>150) & (mj1<200)
        #mjj_maskbox = mjj[(squarebool1) & (mjj>2037)]
        #mjj_masktop = mjj[(squarebool2) & (mjj>2037)]
        #mjj_others = mjj[(np.logical_not(np.logical_or(squarebool1,squarebool2))) & (mjj>2037)]


        hf = h5py.File(self.tagname+".h5", 'w')
        #hf.create_dataset('mjj', data=np.asarray(themassarray))
        hf.create_dataset('mjj', data=np.asarray(mass[thefinalfinalbool]))
        #hf.create_dataset('mjj', data=np.asarray(mass[(thefinalfinalbool) & (np.logical_not(np.logical_or(squarebool1,squarebool2)))]))
        hf.create_dataset('loss1', data=np.asarray(loss1[thefinalfinalbool]))
        hf.create_dataset('loss2', data=np.asarray(loss2[thefinalfinalbool]))
        hf.create_dataset('mj1', data=np.asarray(mj1[thefinalfinalbool]))
        hf.create_dataset('mj2', data=np.asarray(mj2[thefinalfinalbool]))
        hf.create_dataset('evt_idx', data=np.asarray(evt_idx[thefinalfinalbool]))
        hf.create_dataset('samp_label', data=np.asarray(samp_label[thefinalfinalbool]))
        #hf.create_dataset('fracSig', data=np.asarray([combinedFracSig]))
        hf.create_dataset('fracSig', data=np.asarray([0.]))
        hf.create_dataset('avgBkgLoss', data=np.asarray([0.]))
        hf.create_dataset('avgSigLoss', data=np.asarray([0.]))
        hf.create_dataset('numEv', data=np.asarray([len(themassarray)]))
        #hf.create_dataset('truth_label', data=np.asarray([0.]))
        hf.create_dataset('truth_label', data=np.asarray(label[thefinalfinalbool]))
        hf.close()        


    def runFit(self):
        run_fit = True
        last_change = 'start'

        #fit_start = 0.5*self.masspoint
        #fit_stop = -1
        #mult_min = 0.5

        #mult_min = 0.45
        #mult_min = 0.3
        mult_min = 1550

        mult_max = 7000
        
        #cmd = "python "+self.templateTop+"/dijetfit.py -i "+self.tagname+".h5 -M "+str(self.masspoint)+" --sig_shape "+self.templateTop+"/"+self.template+"/graviton_interpolation_M"+str(self.masspoint)+".0.root --dcb-model -p plots_M"+str(self.masspoint)+"_"+self.tagname+" --sig_norm_unc 0.15"
        #cmd = "python "+self.templateTop+"/dijetfit.py -i "+self.tagname+".h5 -M "+str(self.masspoint)+" --sig_shape "+self.templateTop+"/"+self.template+"/graviton_interpolation_M"+str(self.masspoint)+".0.root --dcb-model -p plots_M"+str(self.masspoint)+"_"+self.tagname+" --sig_norm_unc 0.17"
        #For 3000
        #For Wp: sig_norm_uncert 0.31
        #For Xyy 400-400: 0.15
        #For Qstar 400: 0.52
        #For Wkk 400: 0.27
        #For ZTT 400: 0.33
        #For YHH 400: 0.26
        

        #with open(self.templateTop+"/signalsEffJson.json","r") as sigjson2:
        if(self.samSignalName != ''):
            with open(self.templateTop+"/"+self.samSignalName+"_selectionEfficiency.json","r") as sigjson2:
                jsInfo2 = json.load(sigjson2)
        if not self.signalName == '':
            #sig_uncert = jsInfo2[self.signalName]["uncertainty"]
            sig_eff = jsInfo2[self.jsonKey]["selectionEfficiency"]
        else:
            #sig_uncert = 0.15
            sig_eff = 1.

        if(self.samSignalName != ''):
            with open(self.templateTop+"/"+self.samSignalName+"_signalSystematic.json","r") as sigjson2:
                jsInfo2 = json.load(sigjson2)
        if not self.signalName == '':
            #sig_uncert = jsInfo2[self.signalName]["uncertainty"]
            sig_uncert = jsInfo2[self.jsonKey]["fullSystematic"]
        else:
            sig_uncert = 0.15
            #sig_eff = 1.


        cmd = "timeout 3m python "+self.templateTop+"/dijetfit.py -i "+self.tagname+".h5 -M "+str(self.masspoint)+" --sig_shape "+self.templateTop+"/"+self.template+"/"+self.templateHeader+"_interpolation_M"+str(self.masspoint)+".0.root --dcb-model -p plots_M"+str(self.masspoint)+"_"+self.tagname+" --sig_norm_unc "+str(sig_uncert)

        iterator=0
        #while(run_fit and iterator < 10):
        while(run_fit and iterator < 30):
            iterator+=1
            
            cmdtorun = cmd
            """
            if(fit_start > 0):
                cmdtorun += " --mjj_min %.0f" % fit_start

            if(fit_stop > 0):
                cmdtorun += " --mjj_max %.0f" % fit_stop
            """

            cmdtorun += " --mult_min %.2f" % mult_min
            cmdtorun += " --mult_max %.2f" % mult_max
            
            print(cmdtorun)

            if(not self.verbose):
                cmdtorun = cmdtorun+" > /dev/null"
            os.system(cmdtorun)
            
            run_fit = False
            
            fit_file = "plots_M"+str(self.masspoint)+"_"+self.tagname+ "/" + 'fit_results_%.1f.json' % self.masspoint
            if(not os.path.exists(fit_file)):
                print("\nFit didn't converge")
                run_fit = True
            else:
                with open(fit_file, 'r') as f:
                    fit_params = json.load(f, encoding="latin-1")
                #print("fit_params['bkgfit_prob']", fit_params['bkgfit_prob'], "fit_params['sbfit_prob']", fit_params['sbfit_prob'], "fit_params['bkgfit_frac_err']", fit_params['bkgfit_frac_err'], 'range from ', fit_start, 'to', fit_stop)
                print("fit_params['bkgfit_prob']", fit_params['bkgfit_prob'], "fit_params['sbfit_prob']", fit_params['sbfit_prob'], "fit_params['bkgfit_frac_err']", fit_params['bkgfit_frac_err'], 'mult_min ', mult_min)

                if((fit_params['bkgfit_prob'] < 0.05 and fit_params['sbfit_prob'] < 0.05) or fit_params['bkgfit_frac_err'] > 0.15):
                #if((fit_params['bkgfit_prob'] < 0.05 and fit_params['sbfit_prob'] < 0.1) or fit_params['bkgfit_frac_err'] > 0.15):
                    run_fit = True

            if(run_fit):
                if(mult_min < 1550.):
                    mult_min = 1550.
                elif(mult_max < 0 or mult_max > 6500.):
                    mult_max = 6500.
                elif(last_change == "end" and (mult_min + 400. <= self.masspoint) ):
                    mult_min += 100.
                    last_change = 'start'
                elif(mult_max - 1000. > self.masspoint):
                    if(mult_max > 5000):
                        #mult_max -= 500.
                        mult_max -= min(500., mult_max-(self.masspoint+450.))
                    else:
                        #mult_max -= 200.
                        mult_max -= min(200., mult_max-(self.masspoint+450.)) 
                    last_change = 'end'
                elif((mult_min + 250. <= self.masspoint) ): 
                    mult_min += 100.
                else:
                    run_fit = False
                """
                #What I was running August 15
                if(self.masspoint - (mult_min+0.05)*self.masspoint > 600):
                    mult_min += 0.05
                else:
                    run_fit = False
                """
                """
                if(run_fit):
                    if(fit_start < 1550.):
                        fit_start = 1550.
                    elif(fit_stop < 0 or fit_stop > 6500.):
                        fit_stop = 6500.
                    elif(last_change == "end" and (fit_start + 400. <= self.masspoint) ): 
                        fit_start += 100.
                        last_change = 'start'
                    elif(fit_stop - 500. >= self.masspoint): 
                        if(fit_stop > 5000):
                            fit_stop -= 500
                        else: fit_stop -= 200.
                        last_change = 'end'
                    elif((fit_start + 250. <= self.masspoint) ): 
                        fit_start += 100.
                    else:
                        print("No boundary changes found!")
                        run_fit = False
                """


        if(self.allLimits):

            file1Sig = uproot.open('myOutput.root')
            tree1Sig = file1Sig['output']
            branchesSig = tree1Sig.arrays()
            massSig = branchesSig['mjj']
            JME_flagSig = branchesSig['JME_flag']
            massSig = massSig[JME_flagSig==0]
            massSig = massSig[massSig>1460.]

            lenMassPassMin = len(massSig)
            massCutRange = np.logical_and(massSig>mult_min, massSig<mult_max)
            massSig = massSig[massCutRange]
            lenMassPassActual = len(massSig)
            massRangeCutEfficiency = float(lenMassPassActual)/float(lenMassPassMin)

            sig_eff *= massRangeCutEfficiency

            #with open("/uscms_data/d1/wmccorma/SamQUAKSpaceMaker/final_oldTraining/masterJson.json","r") as sigjson:
            with open(self.templateTop+"/masterJson.json","r") as sigjson:
                jsInfo = json.load(sigjson)
            if not self.signalName == '':
                presel_eff = jsInfo[self.signalName]["preselection_eff"]
                deta_eff = jsInfo[self.signalName]["d_eta_eff"]
            else:
                presel_eff = 1.
                deta_eff = 1.

            sig_norm = 1680.
            lumi = 138.
            cmdtorun = "timeout 3m python "+self.templateTop+"/dijetfit.py -i "+self.tagname+".h5 -M "+str(self.masspoint)+" --sig_shape "+self.templateTop+"/"+self.template+"/"+self.templateHeader+"_interpolation_M"+str(self.masspoint)+".0.root --dcb-model -p plots_M"+str(self.masspoint)+"_"+self.tagname+" --sig_norm_unc 0."
            cmdtorun += " --mult_min %.2f" % mult_min
            cmdtorun += " --mult_max %.2f" % mult_max
            if(not self.verbose):
                cmdtorun = cmdtorun+" > /dev/null"
            os.system(cmdtorun)
            f_signif = r.TFile.Open('higgsCombinelim_test_raw.AsymptoticLimits.mH'+str(self.masspoint)+'.root', "READ")
            res1 = f_signif.Get("limit")
            eps = 0.01
            try:
                signif = [0,0,0,0,0,0]
                for i in range(6):
                    res1.GetEntry(i)
                    if(res1.quantileExpected == -1):  # obs limit
                        signif[5] = res1.limit
                    elif(abs(res1.quantileExpected - 0.5) < eps):  # exp limit
                        signif[2] = res1.limit
                    elif(abs(res1.quantileExpected - 0.025) < eps):  # 2sigma, low
                        signif[0] = res1.limit
                    elif(abs(res1.quantileExpected - 0.16) < eps):  # 1sigma, low
                        signif[1] = res1.limit
                    elif(abs(res1.quantileExpected - 0.84) < eps):  # 1sigma, high
                        signif[3] = res1.limit
                    elif(abs(res1.quantileExpected - 0.975) < eps):  # 2sigma, high
                        signif[4] = res1.limit
            except:
                signif = [0,0,0,0,0,0]

            print("no systematics, but including efficiencies and whatnot")
            print("Obs limit is %.3f (%.1f events)" % (signif[5], signif[5]*sig_norm))
            print("Expected was %.3f (%.1f events)" % (signif[2], signif[2]*sig_norm))
            print("Expected range %.1f-%.1f (one sigma), %.1f-%.1f (two sigma)" % (signif[1] * sig_norm, signif[3]*sig_norm, signif[0] * sig_norm, signif[4] * sig_norm))

            print("Full: Obs limit is %.3f (%.2f fb)" % (signif[5]*sig_norm/(presel_eff*deta_eff*sig_eff), signif[5]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi) ))
            print("Full: Expected was %.3f (%.2f fb)" % (signif[2]*sig_norm/(presel_eff*deta_eff*sig_eff), signif[2]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi) ))
            print("Full: Expected range %.2f-%.2f (one sigma), xs %.2f-%.2f (one sigma)" % ( signif[1] * sig_norm/(presel_eff*deta_eff*sig_eff), signif[3]*sig_norm/(presel_eff*deta_eff*sig_eff), signif[1] * sig_norm/(presel_eff*deta_eff*sig_eff*lumi), signif[3]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi) ))

            noSystDict = {"obs": signif[5]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi), "exp": signif[2]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi), "exp-1": signif[1] * sig_norm/(presel_eff*deta_eff*sig_eff*lumi), "exp+1": signif[3]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi)}
            with open(self.templateTop+"/sampleNamesEnhanced.json","r") as sigjson:
                sampleNameJson = json.load(sigjson)
            with open(sampleNameJson[self.jsonKey]["specificSignalName"]+"_noSyst.json", "w") as outfile:
                json.dump(noSystDict, outfile)

            cmdtorun = "timeout 3m python "+self.templateTop+"/dijetfit.py -i "+self.tagname+".h5 -M "+str(self.masspoint)+" --sig_shape "+self.templateTop+"/"+self.template+"/"+self.templateHeader+"_interpolation_M"+str(self.masspoint)+".0.root --dcb-model -p plots_M"+str(self.masspoint)+"_"+self.tagname+" --sig_norm_unc "+str(sig_uncert)
            cmdtorun += " --mult_min %.2f" % mult_min
            cmdtorun += " --mult_max %.2f" % mult_max
            if(not self.verbose):
                cmdtorun = cmdtorun+" > /dev/null"
            os.system(cmdtorun)
            f_signif = r.TFile.Open('higgsCombinelim_test_raw.AsymptoticLimits.mH'+str(self.masspoint)+'.root', "READ")
            res1 = f_signif.Get("limit")
            try:
                signif = [0,0,0,0,0,0]
                for i in range(6):
                    res1.GetEntry(i)
                    if(res1.quantileExpected == -1):  # obs limit
                        signif[5] = res1.limit
                    elif(abs(res1.quantileExpected - 0.5) < eps):  # exp limit
                        signif[2] = res1.limit
                    elif(abs(res1.quantileExpected - 0.025) < eps):  # 2sigma, low
                        signif[0] = res1.limit
                    elif(abs(res1.quantileExpected - 0.16) < eps):  # 1sigma, low
                        signif[1] = res1.limit
                    elif(abs(res1.quantileExpected - 0.84) < eps):  # 1sigma, high
                        signif[3] = res1.limit
                    elif(abs(res1.quantileExpected - 0.975) < eps):  # 2sigma, high
                        signif[4] = res1.limit
            except:
                signif = [0,0,0,0,0,0]

            print("WITH systematics, but including efficiencies and whatnot")
            print("Obs limit is %.3f (%.1f events)" % (signif[5], signif[5]*sig_norm))
            print("Expected was %.3f (%.1f events)" % (signif[2], signif[2]*sig_norm))
            print("Expected range %.1f-%.1f (one sigma), %.1f-%.1f (two sigma)" % (signif[1] * sig_norm, signif[3]*sig_norm, signif[0] * sig_norm, signif[4] * sig_norm))

            print("Full: Obs limit is %.3f (%.2f fb)" % (signif[5]*sig_norm/(presel_eff*deta_eff*sig_eff), signif[5]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi) ))
            print("Full: Expected was %.3f (%.2f fb)" % (signif[2]*sig_norm/(presel_eff*deta_eff*sig_eff), signif[2]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi) ))
            print("Full: Expected range %.2f-%.2f (one sigma), xs %.2f-%.2f (one sigma)" % ( signif[1] * sig_norm/(presel_eff*deta_eff*sig_eff), signif[3]*sig_norm/(presel_eff*deta_eff*sig_eff), signif[1] * sig_norm/(presel_eff*deta_eff*sig_eff*lumi), signif[3]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi) ))

            withSystDict = {"obs": signif[5]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi), "exp": signif[2]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi), "exp-1": signif[1] * sig_norm/(presel_eff*deta_eff*sig_eff*lumi), "exp+1": signif[3]*sig_norm/(presel_eff*deta_eff*sig_eff*lumi)}
            with open(sampleNameJson[self.jsonKey]["specificSignalName"]+".json", "w") as outfile:
                json.dump(withSystDict, outfile)



        if(self.evsInName):
            os.system("mv "+'higgsCombinesignificance_test_raw.Significance.mH'+str(self.masspoint)+'.root '+'higgsCombinesignificance_test_raw_'+str(self.reqEvs)+'evs.Significance.mH'+str(self.masspoint)+'.root')
            os.system("mv "+'higgsCombinepvalue_test_raw.Significance.mH'+str(self.masspoint)+'.root '+'higgsCombinepvalue_test_raw_'+str(self.reqEvs)+'evs.Significance.mH'+str(self.masspoint)+'.root')
            os.system("mv "+'higgsCombinelim_test_raw.AsymptoticLimits.mH'+str(self.masspoint)+'.root '+'higgsCombinelim_test_raw_'+str(self.reqEvs)+'evs.AsymptoticLimits.mH'+str(self.masspoint)+'.root')
            f_signif = r.TFile.Open('./higgsCombinesignificance_test_raw_'+str(self.reqEvs)+'evs.Significance.mH'+str(self.masspoint)+'.root', "READ")
        else:
            f_signif = r.TFile.Open('./higgsCombinesignificance_test_raw.Significance.mH'+str(self.masspoint)+'.root', "READ")    
        res1 = f_signif.Get("limit")
        try:
            res1.GetEntry(0)
            signif = res1.limit
        except:
            signif = 0.
        print('SIGNIFICANCE ', signif)
        if(signif>3.):
            print('ABOVE 3 SIGMA')
        if(self.evsInName):
            f_signif2 = r.TFile.Open('./higgsCombinepvalue_test_raw_'+str(self.reqEvs)+'evs.Significance.mH'+str(self.masspoint)+'.root', "READ")
        else:
            f_signif2 = r.TFile.Open('./higgsCombinepvalue_test_raw.Significance.mH'+str(self.masspoint)+'.root', "READ")    
        res2 = f_signif2.Get("limit")
        try:
            res2.GetEntry(0)
            signif2 = res2.limit
        except:
            signif2 = 0.
        print('PVAL ', signif2)


    def runBiasTest(self):
        fit_file =  "plots_M"+str(self.masspoint)+"_"+self.tagname + '/fit_results_%.1f.json' % self.masspoint
        print(fit_file)
        if(not os.path.exists(fit_file)):
            print("Fit results file not found. Run regular fit before bias test!")
            exit(1)
        outdir = "plots_M"+str(self.masspoint)+"_"+self.tagname
        with open(fit_file, 'r') as f:
            fit_params = json.load(f, encoding="latin-1")
            exp_lim = fit_params['exp_lim_events']

        alt_shape_ver = 2
        num_samples = 200
        dijet_cmd = "python "+self.templateTop+"/bias_test.py -i "+self.tagname+".h5 -p plots_M"+str(self.masspoint)+"_"+self.tagname+" --sig_norm_unc 0.15 --rebin --alt_shape_ver %i --num_samples %i " % (alt_shape_ver, num_samples)
        
        #roughly inject 0, 3sigma and 5sigma
        #sigs = [0., 2.5 * exp_lim, 4.0 * exp_lim]
        #labels = ["%isigma" % (sig) for sig in [0, 3, 5]]
        
        sigs = [0.0 * exp_lim, 1.5 * exp_lim, 4.0 * exp_lim]
        labels = ["%isigma" % (sig) for sig in [0, 2, 5]]        
        
        for i,num_sig in enumerate(sigs):
            dijet_cmd_iter =  dijet_cmd + " -M "+str(self.masspoint)+" --sig_shape "+self.templateTop+"/"+self.template+"/"+self.templateHeader+"_interpolation_M"+str(self.masspoint)+".0.root --dcb-model -l %s --num_sig %i " % (labels[i], num_sig)
            dijet_cmd_iter += " >& "+"plots_M"+str(self.masspoint)+"_"+self.tagname+"/bias_test_log%.0f_nsig%i.txt " % (self.masspoint, num_sig)
            full_fit_cmd = dijet_cmd_iter
            print(full_fit_cmd)
            os.system(full_fit_cmd)



    def run(self):
        os.chdir(self.directory)

        if(self.useWSBins):
            for b in binDict:
                if(self.masspoint in binDict[b][0]):
                    self.WSbin = b
                    self.WSbinLOW = b[0]+str(int(b[1])-1)
                    self.WSbinHIGH = b[0]+str(int(b[1])+1)

        if(self.doJMETest):
            for j in range(8):
                print('JME variation ', j+1)
                self.makeH5s(j+1)
            print('JME variation ', 0)
            self.makeH5s(0)
        else:
            self.makeH5s(0)

        """
        binsx = [2037, 2132, 2231, 2332, 2438,
                 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854,
                 4010, 4171, 4337, 4509, 4700, 4900,  5100, 5300, 5500, 5800,
                 6100, 6400, 6800]
        f = h5py.File("BKG_Template_singleBin_req300_2600.h5", "r")
        mjj = np.asarray(f['mjj'])
        mj1 = np.asarray(f['mj1'])
        mj2 = np.asarray(f['mj2'])
        squarebool1 = np.logical_and((mj1>280) & (mj1<340), (mj2>260) & (mj2<320))
        squarebool2 = (mj1>160) & (mj1<200)
        mjj_maskbox = mjj[(squarebool1) & (mjj>2037)]
        mjj_masktop = mjj[(squarebool2) & (mjj>2037)]
        mjj_others = mjj[(np.logical_not(np.logical_or(squarebool1,squarebool2))) & (mjj>2037)]
        print('testing len with mjj>2037', len(mjj[mjj>2037]), 'and len(mjj_maskbox)+len(mjj_masktop)', len(mjj_maskbox)+len(mjj_masktop), 'len(mjj_others)', len(mjj_others))
        """
        """
        boxhist = r.TH1F("mjj in box", "mjj in box", len(binsx) -1, array('d', binsx))
        tophist = r.TH1F("mjj in top", "mjj in top", len(binsx) -1, array('d', binsx))
        othershist = r.TH1F("mjj elsewhere", "mjj elsewhere", len(binsx) -1, array('d', binsx))
        for m in mjj_maskbox:
            boxhist.Fill(m)
        for m in mjj_masktop:
            tophist.Fill(m)
        for m in mjj_others:
            othershist.Fill(m)
        """

        if(not self.justPlots):
            self.runFit()

            if(self.doBiasTest):
                self.runBiasTest()

        os.chdir(self.templateTop)












parser = optparse.OptionParser()
parser.add_option("-d", "--dir", dest="directory", default=".")
parser.add_option("-t", "--template", dest="template", default="sigTemplateMakerWithInterpolation")
parser.add_option("-m","--masspoint", dest="masspoint", type=int, default=-1, help="starting mass for scan")
parser.add_option("-l","--binMin", dest="binMin", type=int, default=3000, help="lower boundary of mass scan")
parser.add_option("-u","--binMax", dest="binMax", type=int, default=0, help="upper boundary of mass scan")
#parser.add_option("-o","--sigOffset", dest="sigOffset", type=float, default=-7.5, help="signal loss offset")
#parser.add_option("-o","--sigOffset", dest="sigOffset", type=float, default=-8.5, help="signal loss offset")
parser.add_option("-o","--sigOffset", dest="sigOffset", type=float, default=-9, help="signal loss offset")
#parser.add_option("-c","--bkgCut", dest="bkgCut", type=float, default=-0.5, help="bkg loss cut")
parser.add_option("-b", "--doBiasTest", action="store_true", dest="doBiasTest", default=False, help="do bias test?")
parser.add_option("-e", "--doHighEta", action="store_true", dest="doHighEta", default=False, help="do high eta sideband?")
parser.add_option("-g", "--doBKGOnly", action="store_true", dest="doBKGOnly", default=False, help="do background only?")
parser.add_option("-s", "--useWSBins", action="store_true", dest="useWSBins", default=False, help="use bins from weak supervision?")
parser.add_option("-w", "--doWeightingTest", action="store_true", dest="doWeightingTest", default=False, help="do weight test?")
parser.add_option("-j", "--doJMETest", action="store_true", dest="doJMETest", default=False, help="do JME test?")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="don't print status messages to stdout")
#parser.add_option("-r","--reqEvs", dest="reqEvs", type=int, default=300, help="required number of events in SR")
parser.add_option("-r","--reqEvs", dest="reqEvs", type=int, default=1500, help="required number of events in SR")
parser.add_option("-n", "--evsInName", action="store_true", dest="evsInName", default=False, help="add reqEvs into tag name?")
parser.add_option("-p", "--justPlots", action="store_true", dest="justPlots", default=False, help="Only make event-related plots?")
parser.add_option("-a", "--allLimits", action="store_true", dest="allLimits", default=False, help="Run limits with no uncerts?")
parser.add_option("-c", "--signalName", dest="signalName", default="")
parser.add_option("-f", "--samSignalName", dest="samSignalName", default="")
parser.add_option("-k", "--jsonKey", dest="jsonKey", default="1000")
parser.add_option("-q","--sigEvs", dest="sigEvs", type=int, default=1000, help="")
parser.add_option("-i", "--saveEffs", action="store_true", dest="saveEffs", default=False, help="Save efficiency json?")
parser.add_option("-x", "--saveSysts", action="store_true", dest="saveSysts", default=False, help="Save systematics json?")


(options, args) = parser.parse_args()

binMin = options.binMin
if(options.binMax<options.binMin):
    binMax = options.binMin + 100
else:
    binMax = options.binMax

for mass in range(binMin,binMax,100):
    runner = QUAKRunner(options, mass)
    runner.run()


sys.exit("Don't need to go further")

