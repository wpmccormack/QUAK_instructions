import re
import os
import json
import sys
import numpy as np
import itertools
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from tqdm import tqdm
from scipy.stats import chi2, chi
from scipy.special import erf, gammainc, gammaincc

from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression

import torch
            
class LossTransform:
    def __init__(self,function,name):
        self.function = function
        self.name = name
    def transform(self,loss):
        return self.function(loss)
    
# Defining loss reductions (add yours here!)
def normLoss(order,name,offset=0):
    f = lambda x : -np.linalg.norm(x+offset,ord=order,axis=1)
    return LossTransform(f,name)

def signedNorm(order,name):
    def f(x):
        y = np.sum(np.sign(x)*np.abs(x)**order,axis=1)
        y = np.sign(y)*np.abs(y)**(1/order)
        return y
    return LossTransform(f,name)

def chi2zscore(x):
    df = x.shape[1]
    y = np.sum(x**2,axis=1)
    return -(y-chi2.mean(df))/chi2.std(df)

def chi2logsf(x):
    df = x.shape[1]
    y = np.sum(x**2,axis=1)
    return chi2.logsf(y,df)

def chizscore(x):
    df = x.shape[1]
    y = np.sqrt(np.sum(x**2,axis=1))
    return -(y-chi.mean(df))/chi.std(df)

def chilogsf(x):
    df = x.shape[1]
    y = np.sqrt(np.sum(x**2,axis=1))
    return chi.logsf(y,df)

def chicdf(x):
    df = x.shape[1]
    y = np.sqrt(np.sum(x**2,axis=1))
    return 1-chi.cdf(y,df)

def radialGauss(x):
    n = x.shape[1]
    r = np.sqrt(np.sum(x**2,axis=1))
    return 1-gammainc(n/2,r**2/2)
    

minLoss = LossTransform(lambda x : np.min(x,axis=1),"Min")
sumLoss = LossTransform(lambda x : np.sum(x,axis=1)/x.shape[1],"Sum")
sumProbLoss = LossTransform(lambda x : -np.log(np.sum(np.exp(-x),axis=1)), "SumProb")
logSumExpAbs = LossTransform(lambda x : -np.log(np.sum(np.exp(np.abs(x)),axis=1)), "LogSumExpAbs")
sumExpAbs = LossTransform(lambda x : -np.sum(np.exp(np.abs(x)),axis=1), "SumExpAbs")
Identity = LossTransform(lambda x : x.flatten(), "None")
sumProbOnlyNeg = LossTransform(lambda x : -np.log(np.sum(np.where(x<0,np.exp(-x),0),axis=1)), "SumProbNeg")

reductions = {
    "inv2norm":normLoss(-2,"inv2norm_o10"),
    "inv5norm":normLoss(-5,"inv5norm_o10"),
    "inv10norm":normLoss(-10,"inv10norm_o10"),
    "inv2norm_o10":normLoss(-2,"inv2norm_o10",offset=10),
    "inv5norm_o10":normLoss(-5,"inv5norm_o10",offset=10),
    "inv10norm_o10":normLoss(-10,"inv10norm_o10",offset=10),
    "l5":normLoss(5,"l5"),
    "l4":normLoss(4,'l4'),
    "l3":normLoss(3,"l3"),
    "l2":normLoss(2,"l2"),
    "l1":normLoss(1,"l1"),
    "Min":minLoss,
    "Sum":sumLoss,
    "SumProb":sumProbLoss,
    "SumProbNeg":sumProbOnlyNeg,
    "SumExpAbs":sumExpAbs,
    "LogSumExpAbs":logSumExpAbs,
    "chi2z":LossTransform(chi2zscore,'chi2z'),
    "chi2sf":LossTransform(chi2logsf,'chi2sf'),
    "chiz":LossTransform(chizscore,'chiz'),
    "chisf":LossTransform(chilogsf,'chisf'),
    "chicdf":LossTransform(chicdf,'chicdf'),
    "radialGauss":LossTransform(radialGauss,"radialGauss"),
    "None":Identity
}

for n in range(0,11):
    reductions[f"signedL{n}"] = signedNorm(n,f"signedL{n}")

# Defining loss transformations (add yours here!)
transformations = {
    "None":Identity
}

# convenience map between systematic string and array index
syst_map = {
    "nom_weight":0,
    "pdf_up":1,
    "pdf_down":2,
    "prefire_up":3,
    "prefire_down":4,
    "pileup_up":5,
    "pileup_down":6,
    "btag_up":7,
    "btag_down":8,
    "PS_ISR_up":9,
    "PS_ISR_down":10,
    "PS_FSR_up":11,
    "PS_FSR_down":12,
    "F_up":13,
    "F_down":14,
    "R_up":15,
    "R_down":16,
    "RF_up":17,
    "RF_down":18,
    "top_ptrw_up":19,
    "top_ptrw_down":20
}

lund_syst_map = {
    "lund_nom":0,
    "lund_plane_up":1,
    "lund_plane_down":2,
    "lund_bjet_up":3,
    "lund_bjet_down":4
}

inv_syst_map = {idx:name for name,idx in syst_map.items()}
inv_lund_syst_map = {idx:name for name,idx in lund_syst_map.items()}

JME_syst_map = {
    "nominal":0,
    "JES_up":1,
    "JES_down":2,
    "JER_up":3,
    "JER_down":4,
    "JMS_up":5,
    "JMS_down":6,
    "JMR_up":7,
    "JMR_down":8
}

class QuakSpaceMaker:
    def __init__(self,bkg_axis,sig_axes,toInject,nInject=None,xs_inject=None,loss_reduction="None",loss_transform="None",decorrelation=False,loss_dir="lossEvals",mode=0):
        """
        bkg_axis = background training to use for the background axis
        sig_axes = signal training(s) to use for signal axis/axes
        toInject = signal to inject
        nInject = Number of signal events to inject into the SR (dEta < 1.3) -- NOTE all events are pre-selected to be in SR
        xs_inject = Signal cross section to inject in fb (number of events derived from presel_eff and dEta_cut_eff)
        loss_transform = transformation to apply to all losses before making QUAK space
        loss_reduction = reduction to apply to signal losses (e.g. sum, norm, etc.)
        """
        assert type(bkg_axis) == str
        assert type(sig_axes) == list
        self.bkg_axis = bkg_axis
        self.sig_axes = sig_axes
        self.axes = [bkg_axis] + sig_axes
        self.toInject = toInject
        self.decorrelation = decorrelation
        self.loss_dir = loss_dir
        self.JME_systs = ["JES_up","JES_down","JER_up","JER_down","JMS_up","JMS_down","JMR_up","JMR_down"]
        
        if not (nInject or xs_inject):
            print("You need to specify how much signal to inject")
        if nInject and xs_inject:
            print("Can't specify nInject and xs_inject simultaneously")
            
        self.loss_transform = transformations[loss_transform]
        self.loss_reduction = reductions[loss_reduction]

        #self.outBase = f"injectedQuakSpaces_noPCA/transform-{self.loss_transform.name}_reduce-{self.loss_reduction.name}/"
        self.outBase = f"injectedQuakSpaces/transform-{self.loss_transform.name}_reduce-{self.loss_reduction.name}/"
        if not os.path.isdir(self.outBase):
            os.makedirs(self.outBase)
        
        self.outDir = self.outBase + f"{self.bkg_axis}-" + "-".join(self.sig_axes) + "/"
        if not os.path.isdir(self.outDir):
            os.makedirs(self.outDir)
            
        with open("masterJson.json","r") as f:
            self.jsInfo = json.load(f)
            
        # calculate nEvents to inject for each JME syst variation
        self.inject_systs = ["nominal"] + self.JME_systs
        self.nInject = {}
        self.xs_inject = {}
        self.lumi = 26.81 # combined 2016-2018 background luminosity
        for syst_var in self.inject_systs:
            if "JES" in syst_var or "JER" in syst_var:
                presel_eff = self.jsInfo[self.toInject][f"preselection_eff_{syst_var}"]
            else:
                presel_eff = self.jsInfo[self.toInject]["preselection_eff"]
            deta_eff = self.jsInfo[self.toInject]["d_eta_eff"]
            if nInject:
                self.nInject[syst_var] = nInject
                self.xs_inject[syst_var] = nInject/(self.lumi * presel_eff * deta_eff)
                print('self.lumi = ', self.lumi, 'presel_eff = ', presel_eff, 'deta_eff = ', deta_eff)
            else:
                self.xs_inject[syst_var] = xs_inject
                self.nInject[syst_var] = int(xs_inject * self.lumi * presel_eff * deta_eff)
        
        if nInject:
            print(f"Injecting {nInject} events / {list(self.xs_inject.values())} fb")
        else:
            print(f"Injecting {xs_inject} fb / {list(self.nInject.values())} events")
        print(f"Bkg axis : {self.bkg_axis}")
        print(f"Sig axes : {self.sig_axes}")
        print(f"Sig Reduction : {loss_reduction}")
        print(f"bkg_loss / sig_loss decorrelation : {decorrelation}")
        
        if nInject:
            self.outFile = self.outDir+f"inject-{toInject}_N{nInject}.root"
        else:
            self.outFile = self.outDir+f"inject-{toInject}_xs{xs_inject:.2f}.root"
        
        if mode==0:
            self.bkgLossFile = f"{self.loss_dir}/BKG.h5"
            self.bkgOutFile = self.outDir+f"BKG.root"
        elif mode==1:
            self.bkgLossFile = f"{self.loss_dir}/DATA_SB.h5"
            self.bkgOutFile = self.outDir+f"DATA_SB.root"
        elif mode==2:
            self.bkgLossFile = f"{self.loss_dir}/DATA.h5"
            self.bkgOutFile = self.outDir+f"DATA.root"
        elif mode==3:
            self.bkgLossFile = f"{self.loss_dir}/BKG_SB.h5"
            self.bkgOutFile = self.outDir+f"BKG_SB.root"
        else:
            print("Unrecognized input for background-type file (0=MC, 1=data sideband, 2=data)")
            exit()
        self.sigLossFile = f"{self.loss_dir}/{toInject}.h5"
        
    def run(self):
        model = torch.load("./xyy_400_400_massCut_fullmodel.pt")
        model.eval()

        # Load bkg losses
        bkg_losses = {}
        bkg_nnscores = {}
        bkg_mjj = None
        with h5py.File(self.bkgLossFile,"r") as f:
            bkg_mjj = f["mjj"][()]
            #if "BKG" in self.bkgLossFile:
            #    bkg_label = f["truth_label"][()]
            #bkg_mj1 = f["jet1_mass"][()]
            #bkg_mj2 = f["jet2_mass"][()]
            #bkg_btag1 = f["jet1_btagscore"][()]
            #bkg_btag2 = f["jet2_btagscore"][()]
            nn_branches = []
            for k in f.keys():
                if "NNScore" in k:
                    nn_branches.append(k)
                    print("    "+k)
            for nnb in nn_branches:
                bkg_nnscores[nnb] = f[nnb][()]
            for ax in self.axes:
                bkg_losses[ax] = self.loss_transform.transform(f[ax][()])
        bkg_JME_label = (-1*np.ones(bkg_mjj.shape[0])).astype('int')
        #bkg_evtIdx = np.arange(bkg_mj1.shape[0])
        
        # load injection losses
        inj_losses = {ax:[] for ax in self.axes}
        inj_nnscores = {}
        inj_mjj = []
        inj_sys = []
        inj_JME_label = []
        lund_matching = None
        inj_lund_sys = []
        inj_lund_pt_var = []
        inj_lund_stat_var = []
        if self.nInject['nominal'] != 0:
            with h5py.File(f"{self.loss_dir}/{self.toInject}.h5","r") as f:
                full_mjj = f["mjj"][()]
                full_sys = f["sys_weights"][()]
                lund_matching = f['lund_weights_matching_unc'][()]
                full_lund_sys = np.concatenate((f['lund_weights'][()].reshape(-1,1),f['lund_weights_sys_var'][()]),axis=1)
                full_lund_pt_var = f['lund_weights_pt_var'][()]
                full_lund_stat_var = f['lund_weights_stat_var'][()]
                # quantities for sampling
                #wgt = full_sys[:,syst_map["nom_weight"]]
                #output_lund_systematics = {s:inj_lund_sys[:,lund_syst_map[s]] for s in lund_syst_map.keys()}
                print('full_sys', full_sys)
                print('full_sys[:,syst_map["nom_weight"]]', full_sys[:,syst_map["nom_weight"]])
                print('full_lund_sys', full_lund_sys)
                print('full_lund_sys[:,lund_syst_map["lund_nom"]]', full_lund_sys[:,lund_syst_map["lund_nom"]])
                #wgt = full_sys[:,syst_map["nom_weight"]]*full_lund_sys[:,lund_syst_map["lund_nom"]]
                wgt = full_sys[:,syst_map["nom_weight"]]
                print('first wgt', wgt, np.average(wgt))
                p = wgt / np.sum(wgt)
                print('first p', p, np.average(p))
                wgt = full_sys[:,syst_map["nom_weight"]]*full_lund_sys[:,lund_syst_map["lund_nom"]]
                print('second wgt', wgt, np.average(wgt))
                p = wgt / np.sum(wgt)
                print('second p', p, np.average(p))
                idx_array = np.arange(len(p))
                # loop over all JME systematics & sample
                nn_branches = []
                for k in f.keys():
                    if "NNScore" in k:
                        nn_branches.append(k)
                        print("    "+k)
                for nnb in nn_branches:
                    inj_nnscores[nnb] = []
                for syst in self.inject_systs:
                    # sampling indices
                    nSample = self.nInject[syst]
                    sample_idx = np.random.choice(idx_array,size=nSample,replace=False,p=p)
                    # taking quantities
                    inj_mjj.append(full_mjj[sample_idx])
                    inj_sys.append(full_sys[sample_idx])
                    inj_lund_sys.append(full_lund_sys[sample_idx])
                    inj_lund_pt_var.append(full_lund_pt_var[sample_idx])
                    inj_lund_stat_var.append(full_lund_stat_var[sample_idx])
                    inj_JME_label.append((JME_syst_map[syst]*np.ones(nSample)).astype('int'))
                    for nnb in nn_branches:
                        inj_nnscores[nnb].append(f[nnb][()][sample_idx])
                    for ax in self.axes:
                        if syst == "nominal":
                            axName = ax
                        else:
                            axName = f"{ax}_{syst}"
                        inj_losses[ax].append(self.loss_transform.transform(f[axName][()][sample_idx]))

        # combine all the JME variations for the injections
        inj_mjj = np.concatenate(inj_mjj,axis=0)
        inj_sys = np.concatenate(inj_sys,axis=0)
        inj_lund_sys = np.concatenate(inj_lund_sys,axis=0)
        inj_lund_pt_var = np.concatenate(inj_lund_pt_var,axis=0)
        inj_lund_stat_var = np.concatenate(inj_lund_stat_var,axis=0)
        inj_JME_label = np.concatenate(inj_JME_label,axis=0)
        for k in inj_nnscores.keys():
            inj_nnscores[k] = np.concatenate(inj_nnscores[k],axis=0)
        for ax in self.axes:
            inj_losses[ax] = np.concatenate(inj_losses[ax],axis=0)
            
        
        # bkg_loss / sig_loss decorrelation
        if self.decorrelation:
            v1_bkg = bkg_losses[self.bkg_axis].reshape(-1,1)
            v2_bkg = np.concatenate([bkg_losses[a].reshape(-1,1) for a in self.sig_axes],axis=1)
            fit_bkg = LinearRegression().fit(v1_bkg,v2_bkg)
            v2corr_bkg = v2_bkg - fit_bkg.predict(v1_bkg)
            
            for ia,a in enumerate(self.sig_axes):
                bkg_losses[a] = v2corr_bkg[:,ia]
            
            if self.nInject['nominal'] != 0:
                v1_inj = inj_losses[self.bkg_axis].reshape(-1,1)
                v2_inj = np.concatenate([inj_losses[a].reshape(-1,1) for a in self.sig_axes],axis=1)
                v2corr_inj = v2_inj - fit_bkg.predict(v1_inj)
                for ia,a in enumerate(self.sig_axes):
                    inj_losses[a] = v2corr_inj[:,ia]
                
        # perform loss reduction
        bkg_loss1 = bkg_losses[self.bkg_axis]
        bkg_loss2 = self.loss_reduction.transform(np.concatenate([bkg_losses[ax].reshape(-1,1) for ax in self.sig_axes],axis=1))
        inj_loss1 = inj_losses[self.bkg_axis]
        if self.nInject['nominal'] != 0:
            inj_loss2 = self.loss_reduction.transform(np.concatenate([inj_losses[ax].reshape(-1,1) for ax in self.sig_axes],axis=1))
        else:
            inj_loss2 = inj_losses[self.bkg_axis]
        
        #loss_sig = np.stack((inj_loss1, inj_loss2), axis=1)
        #loss_bkg = np.stack((bkg_loss1, bkg_loss2), axis=1)
        loss_sig = np.stack((inj_loss2, inj_loss1), axis=1)
        loss_bkg = np.stack((bkg_loss2, bkg_loss1), axis=1)

        inj_ptmodel = model(torch.tensor(loss_sig).float())
        bkg_ptmodel = model(torch.tensor(loss_bkg).float())

        print('inj mean', inj_ptmodel.mean())
        print(inj_ptmodel.shape)
        print('bkg mean', bkg_ptmodel.mean())
        print(bkg_ptmodel.shape)


        # prep nn scores
        nn_avail_inj = list(inj_nnscores.keys())
        nn_avail_bkg = list(bkg_nnscores.keys())
        nn_toSave = [k for k in nn_avail_bkg if k in nn_avail_inj]
        
        # create output file with injected signal events
        import uproot
        n_sig = 0
        if self.nInject['nominal'] != 0:
            print("Making signal injection file")
            n_sig = inj_mjj.shape[0]
            systematics = [s for s in syst_map.keys() if s != "nom_weight"]
            output_systematics = {s:inj_sys[:,syst_map[s]] for s in systematics}
            output_lund_systematics = {s:inj_lund_sys[:,lund_syst_map[s]] for s in lund_syst_map.keys()}
            
            outTree = {}
            outTree["label"] = np.ones(n_sig)
            outTree["JME_flag"] = inj_JME_label
            outTree["mjj"] = inj_mjj
            outTree["loss1"] = inj_loss1
            outTree["loss2"] = inj_loss2
            outTree["ptmodel"] = inj_ptmodel.detach().numpy()
            outTree["samp_label"] = np.ones(n_sig)
            outTree["mj1"] = np.zeros(n_sig)
            outTree["mj2"] = np.zeros(n_sig)
            outTree["btag1"] = np.zeros(n_sig)
            outTree["btag2"] = np.zeros(n_sig)
            outTree["evt_idx"] = -1*np.ones(n_sig)
            for s,sf in output_systematics.items():
                outTree[s] = sf.astype('float')
            for s,sf in output_lund_systematics.items():
                outTree[s] = sf.astype('float')
            outTree['lund_matching'] = np.repeat(lund_matching,n_sig)
            outTree['lund_pt_var'] = inj_lund_pt_var
            outTree['lund_stat_var'] = inj_lund_stat_var
            nnScores_save = {}
            for k in nn_toSave:
                nnScores_save[k] = inj_nnscores[k]
            for k,scores in nnScores_save.items():
                newk = ''
                for l in k:
                    newk += l if l != '-' else '_'
                outTree[newk] = scores
            with uproot.recreate(self.outFile) as fout:
                fout["output"] = outTree

            print(f"Injection output at {self.outFile}")
        
        if not os.path.exists(self.bkgOutFile):
            print("Making large data or bkg file")
            n_bkg = bkg_mjj.shape[0]
            #if "BKG" in self.bkgLossFile:
            #    samp_labels = bkg_label
            #else:
            #    samp_labels = np.ones(bkg_mjj.shape[0])
                
            outTree = {}
            outTree["label"] = np.zeros(n_bkg)
            outTree["mjj"] = bkg_mjj
            outTree["loss1"] = bkg_loss1
            outTree["loss2"] = bkg_loss2
            outTree["ptmodel"] = bkg_ptmodel.detach().numpy()
            #outTree["samp_label"] = samp_labels
            #outTree["mj1"] = bkg_mj1
            #outTree["mj2"] = bkg_mj2
            #outTree["btag1"] = bkg_btag1
            #outTree["btag2"] = bkg_btag2
            #outTree["evt_idx"] = bkg_evtIdx
            nnScores_save = {}
            for k in nn_toSave:
                nnScores_save[k] = bkg_nnscores[k]
            for k,scores in nnScores_save.items():
                newk = ''
                for l in k:
                    newk += l if l != '-' else '_'
                outTree[newk] = scores
            with uproot.recreate(self.bkgOutFile) as fout:
                fout["output"] = outTree

            print(f"Evaluation on data/bkg saved at {self.bkgOutFile}")

        
                
