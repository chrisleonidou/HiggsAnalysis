#!/usr/bin/env python

import re

from HiggsAnalysis.HeavyChHiggsToTauNu.tools.multicrab import *

#step = "skim"
step = "embedding"
#step = "analysis"
#step = "analysisTau"
#step = "signalAnalysis"
#step = "muonAnalysis"
#step = "caloMetEfficiency"

dirPrefix = ""
#dirPrefix += "_caloMet45"
#dirPrefix += "_caloMet60"
#dirPrefix += "_taueff"
#dirPrefix += "_noTauMatching"
#dirPrefix += "_noTauPtCut"
#dirPrefix += "_tauPt50"
#dirPrefix += "_nJet40"
#dirPrefix += "_noEmuVeto"
#dirPrefix += "_noEmuVetoEnd"
#dirPrefix += "_MCGT"
#dirPrefix += "_forClosureTest"
#dirPrefix = "_TauIdScan"
#dirPrefix = "_iso05"
#dirPrefix = "_test"

#pt = "_pt30"
#pt = "_pt40"
#if step in ["generation", "embedding", "analysis", "signalAnalysis"]:
#    dirPrefix += pt

if step == "signalAnalysis":
    #dirPrefix += "_triggerVertex2010"
    #dirPrefix += "_triggerVertex2011"
    #dirPrefix += "_trigger2010"
    #dirPrefix += "_trigger2011"
    pass

config = {"skim":           {"input": "AOD",                           "config": "muonSkim_cfg.py", "output": "skim.root"},
          "embedding":      {"input": "tauembedding_skim_v13", "config": "embed.py",   "output": "embedded.root"},
#          "analysis":       {"input": "tauembedding_embedding_v13"+pt,  "config": "embeddingAnalysis_cfg.py"},
          "analysis":       {"input": "tauembedding_embedding_v13",  "config": "embeddingAnalysis_cfg.py"},
#          "analysisTau":    {"input": "pattuple_v17",                  "config": "tauAnalysis_cfg.py"},
#          "signalAnalysis": {"input": "tauembedding_embedding_v13"+pt,  "config": "../signalAnalysis_cfg.py"},
          "signalAnalysis": {"input": "tauembedding_embedding_v13",  "config": "../signalAnalysis_cfg.py"},
          "muonAnalysis":   {"input": "tauembedding_skim_v13",          "config": "muonAnalysisFromSkim_cfg.py"},
          "caloMetEfficiency": {"input": "tauembedding_skim_v13",         "config": "caloMetEfficiency_cfg.py"},
          }

crabcfg = "crab.cfg"
if step in ["analysis", "analysisTau", "signalAnalysis", "muonAnalysis", "caloMetEfficiency"]:
    crabcfg = "../crab_analysis.cfg"


multicrab = Multicrab(crabcfg, config[step]["config"], lumiMaskDir="..")

datasetsData2010 = [
    "Mu_136035-144114_Apr21", # HLT_Mu9
    "Mu_146428-147116_Apr21", # HLT_Mu9
    "Mu_147196-149294_Apr21", # HLT_Mu15_v1
]
datasetsData2011 = [
    "SingleMu_Mu_160431-163261_May10",  # HLT_Mu20_v1
    "SingleMu_Mu_163270-163869_May10",  # HLT_Mu24_v2
    "SingleMu_Mu_165088-166150_Prompt", # HLT_Mu30_v3
    "SingleMu_Mu_166161-166164_Prompt", # HLT_Mu40_v1
    "SingleMu_Mu_166346-166346_Prompt", # HLT_Mu40_v2
    "SingleMu_Mu_166374-167043_Prompt", # HLT_Mu40_v1
    "SingleMu_Mu_167078-167913_Prompt", # HLT_Mu40_v3
    "SingleMu_Mu_170722-172619_Aug05",  # HLT_Mu40_v5
    "SingleMu_Mu_172620-173198_Prompt", # HLT_Mu40_v5
    "SingleMu_Mu_173236-173692_Prompt", # HLT_Mu40_eta2p1_v1
]
datasetsMCnoQCD = [
    "TTJets_TuneZ2_Summer11",
    "WJets_TuneZ2_Summer11",
    "DYJetsToLL_M50_TuneZ2_Summer11",
    "T_t-channel_TuneZ2_Summer11",
    "Tbar_t-channel_TuneZ2_Summer11",
    "T_tW-channel_TuneZ2_Summer11",
    "Tbar_tW-channel_TuneZ2_Summer11",
    "T_s-channel_TuneZ2_Summer11",
    "Tbar_s-channel_TuneZ2_Summer11",
    "WW_TuneZ2_Summer11",
    "WZ_TuneZ2_Summer11",
    "ZZ_TuneZ2_Summer11",
]
datasetsMCQCD = [
    "QCD_Pt20_MuEnriched_TuneZ2_Summer11",
]
datasetsTest = [
    "TTToHplusBWB_M120_Summer11"
]
    
datasets = []
if step in ["analysis", "analysisTau"]:
    datasets.extend(datasetsMCnoQCD)
else:
#    datasets.extend(datasetsData2010)
    datasets.extend(datasetsData2011)
    datasets.extend(datasetsMCnoQCD)
    datasets.extend(datasetsMCQCD)

#    if step in ["skim", "generation", "embedding", "caloMetEfficiency"]:
#        datasets.extend(datasetsTest)

multicrab.extendDatasets(config[step]["input"], datasets)

multicrab.appendLineAll("GRID.maxtarballsize = 15")
#if step != "skim":
#    multicrab.extendBlackWhiteListAll("ce_white_list", ["jade-cms.hip.fi"])


path_re = re.compile("_tauembedding_.*")
tauname = "_tauembedding_%s_v13_1" % step
#if step in ["generation", "embedding"]:
#    tauname += pt

reco_re = re.compile("^Run[^_]+_(?P<reco>[^_]+_v\d+_[^_]+_)")

# Goal: ~5 hour jobs
skimNjobs = {
    "WJets_TuneZ2_Summer11": 1000, # ~10 hours
    "TTJets_TuneZ2_Summer11": 500,
    "QCD_Pt20_MuEnriched_TuneZ2_Summer11": 490,
    "DYJetsToLL_M50_TuneZ2_Summer11": 1000,
    "T_t-channel_TuneZ2_Summer11": 490,
    "Tbar_t-channel_TuneZ2_Summer11": 160,
    "T_tW-channel_TuneZ2_Summer11": 90,
    "Tbar_tW-channel_TuneZ2_Summer11": 90,
    "T_s-channel_TuneZ2_Summer11": 50,
    "Tbar_s-channel_TuneZ2_Summer11": 10,
    "WW_TuneZ2_Summer11": 200,
    "WZ_TuneZ2_Summer11": 200,
    "ZZ_TuneZ2_Summer11": 350,
    }

muonAnalysisNjobs = { # goal: 30k events/job
    "SingleMu_Mu_160431-163261_May10": 2,
    "SingleMu_Mu_163270-163869_May10": 5,
    "SingleMu_Mu_165088-166150_Prompt": 6,
    "SingleMu_Mu_166161-166164_Prompt": 1,
    "SingleMu_Mu_166346-166346_Prompt": 1,
    "SingleMu_Mu_166374-167043_Prompt": 4,
    "SingleMu_Mu_167078-167913_Prompt": 3,
    "SingleMu_Mu_170722-172619_Aug05": 5,
    "SingleMu_Mu_172620-173198_Prompt": 8,
    "SingleMu_Mu_173236-173692_Prompt": 4,
    
    "WJets_TuneZ2_Summer11": 60,
    "TTJets_TuneZ2_Summer11": 17,
    "QCD_Pt20_MuEnriched_TuneZ2_Summer11": 5,
    "DYJetsToLL_M50_TuneZ2_Summer11": 15,
    "T_t-channel_TuneZ2_Summer11": 3,
    "Tbar_t-channel_TuneZ2_Summer11": 2,
    "T_tW-channel_TuneZ2_Summer11": 4,
    "Tbar_tW-channel_TuneZ2_Summer11": 4,
    "T_s-channel_TuneZ2_Summer11": 1,
    "Tbar_s-channel_TuneZ2_Summer11": 1,
    "WW_TuneZ2_Summer11": 8,
    "WZ_TuneZ2_Summer11": 8,
    "ZZ_TuneZ2_Summer11": 8,
    }
   

def modify(dataset):
    name = ""

    if dataset.isData():
        dataset.appendLine("CMSSW.total_number_of_lumis = -1")
    else:
        dataset.appendLine("CMSSW.total_number_of_events = -1")

    path = dataset.getDatasetPath().split("/")
    if step == "skim":
        name = path[2].replace("-", "_")
        name += "_"+path[3]
        name += tauname

        if dataset.isData():
            frun = dataset.getName().split("_")[1].split("-")[0]
            m = reco_re.search(name)
            name = reco_re.sub(m.group("reco")+frun+"_", name)

        dataset.useServer(False)

        try:
            njobs = skimNjobs[dataset.getName()]
            dataset.setNumberOfJobs(njobs)
#            if njobs > 490:
#                dataset.useServer(True)
        except KeyError:
            pass


    else:
        name = path_re.sub(tauname, path[2])
        name = name.replace("local-", "")

    if dataset.isData() and step in ["generation", "embedding"]:
        dataset.appendArg("overrideBeamSpot=1")

    dataset.appendLine("USER.publish_data_name = "+name)
    dataset.appendLine("CMSSW.output_file = "+config[step]["output"])

def modifyAnalysis(dataset):
    if step == "signalAnalysis":
        dataset.appendArg("tauEmbeddingInput=1")
        dataset.appendArg("doPat=1")
#    if step == "analysisTau":
#        if dataset.getName() == "WJets":
#            dataset.setNumberOfJobs(100)

def modifyMuonAnalysis(dataset):
    try:
        dataset.setNumberOfJobs(muonAnalysisNjobs[dataset.getName()])
    except KeyError:
        pass
    

if step in ["analysis", "analysisTau","signalAnalysis"]:
    multicrab.appendLineAll("CMSSW.output_file = histograms.root")
    multicrab.forEachDataset(modifyAnalysis)
elif step in ["muonAnalysis", "caloMetEfficiency"]:
    multicrab.appendLineAll("CMSSW.output_file = histograms.root")
    multicrab.forEachDataset(modifyMuonAnalysis)
else:
    multicrab.forEachDataset(modify)

multicrab.extendBlackWhiteListAll("se_black_list", defaultSeBlacklist)

prefix = "multicrab_"+step+dirPrefix

multicrab.createTasks(prefix=prefix)
#multicrab.createTasks(configOnly=True,prefix=prefix)
