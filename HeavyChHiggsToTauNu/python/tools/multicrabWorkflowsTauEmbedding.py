## \package multicrabDatasetsTauEmbedding
# Functions for embedding workflow definitions

import re

from multicrabWorkflowsTools import Dataset, Workflow, Data, Source, updatePublishName, TaskDef, updateTaskDefinitions
import multicrabDatasetsCommon as common
from multicrabWorkflowsPattuple import constructProcessingWorkflow_44X, constructProcessingWorkflow_53X

def addEmbeddingGenTauSkim_44X(version, datasets, updateDefinitions):
    # Tau+MET trigger has 5 % efficiency, GenTauSkim has 10 %, so 2x jobs
    defaultDefinitions = {
        "TTJets_TuneZ2_Fall11":              TaskDef(njobsIn=4000, njobsOut=50),
        }
    workflowName = "tauembedding_gentauskim_"+version
    updateTaskDefinitions(defaultDefinitions, updateDefinitions, workflowName)
    for datasetName, taskDef in defaultDefinitions.iteritems():
        dataset = datasets.getDataset(datasetName)
        wf = constructProcessingWorkflow_44X(dataset, taskDef, sourceWorkflow="AOD", workflowName=workflowName)
        if dataset.isData():
            raise Exception("GenTauSkim workflow is not supported for data")
        wf.addCrabLine("CMSSW.total_number_of_events = -1")
        name = updatePublishName(dataset, wf.source.getDataForDataset(dataset).getDatasetPath(), workflowName)
        wf.addCrabLine("USER.publish_data_name = "+name)

        wf.addArg("customizeConfig", "CustomGenTauSkim")

        dataset.addWorkflow(wf)
        if wf.output is not None:
            dataset.addWorkflow(Workflow("tauembedding_gentauanalysis_"+version, source=Source(workflowName),
                                         args=wf.args, output_file="histograms.root"))

def addEmbeddingGenTauSkim_44X(version, datasets, updateDefinitions):
    # Tau+MET trigger has 5 % efficiency, GenTauSkim has 10 %, so 2x jobs
    defaultDefinitions = {
        "TTJets_TuneZ2_Fall11":              TaskDef(njobsIn=4000, njobsOut=50),
        }
    workflowName = "tauembedding_gentauskim_"+version
    updateTaskDefinitions(defaultDefinitions, updateDefinitions, workflowName)
    for datasetName, taskDef in defaultDefinitions.iteritems():
        dataset = datasets.getDataset(datasetName)
        wf = constructProcessingWorkflow_44X(dataset, taskDef, sourceWorkflow="AOD", workflowName=workflowName)
        if dataset.isData():
            raise Exception("GenTauSkim workflow is not supported for data")
        wf.addCrabLine("CMSSW.total_number_of_events = -1")
        name = updatePublishName(dataset, wf.source.getDataForDataset(dataset).getDatasetPath(), workflowName)
        wf.addCrabLine("USER.publish_data_name = "+name)

        wf.addArg("customizeConfig", "CustomGenTauSkim")

        dataset.addWorkflow(wf)
        if wf.output is not None:
            dataset.addWorkflow(Workflow("tauembedding_gentauanalysis_"+version, source=Source(workflowName),
                                         args=wf.args, output_file="histograms.root", crabLines=["CMSSW.total_number_of_lumis = -1"]))

def addEmbeddingAodAnalysis_44X(datasets):
    njobs = {
        "WJets_TuneZ2_Fall11":               TaskDef(njobsIn=490),
#        "W2Jets_TuneZ2_Fall11":              TaskDef(njobsIn=300),
#        "W3Jets_TuneZ2_Fall11":              TaskDef(njobsIn=120),
#        "W4Jets_TuneZ2_Fall11":              TaskDef(njobsIn=200),
        "TTJets_TuneZ2_Fall11":              TaskDef(njobsIn=1000),
        "DYJetsToLL_M50_TuneZ2_Fall11":      TaskDef(njobsIn=350),
        "T_t-channel_TuneZ2_Fall11":         TaskDef(njobsIn=50),
        "Tbar_t-channel_TuneZ2_Fall11":      TaskDef(njobsIn=50),
        "T_tW-channel_TuneZ2_Fall11":        TaskDef(njobsIn=20),
        "Tbar_tW-channel_TuneZ2_Fall11":     TaskDef(njobsIn=20),
        "T_s-channel_TuneZ2_Fall11":         TaskDef(njobsIn=10),
        "Tbar_s-channel_TuneZ2_Fall11":      TaskDef(njobsIn=10),
        "WW_TuneZ2_Fall11":                  TaskDef(njobsIn=50),
        "WZ_TuneZ2_Fall11":                  TaskDef(njobsIn=50),
        "ZZ_TuneZ2_Fall11":                  TaskDef(njobsIn=50),
#        "QCD_Pt20_MuEnriched_TuneZ2_Fall11": TaskDef(njobsIn=20),
        }

    for datasetName, taskDef in njobs.iteritems():
        dataset = datasets.getDataset(datasetName)
        source = Source("AOD", number_of_jobs=taskDef.njobsIn)
        wf = Workflow("embeddingAodAnalysis_44X", source=source, output_file="histograms.root")
        wf.addCrabLine("CMSSW.total_number_of_lumis = -1")
        wf.addArg("doPat", 1)
        dataset.addWorkflow(wf)

def getDefaultDefinitions_44X():
    mcTrigger = "HLT_Mu40_eta2p1_v1"
    def TaskDefMC(**kwargs):
        return TaskDef(triggerOR=[mcTrigger], **kwargs)

    return {
        "SingleMu_160431-163261_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu20_v1"]),
        "SingleMu_163270-163869_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu24_v2"]),
        "SingleMu_165088-166150_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu30_v3"]),

        "SingleMu_166161-166164_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu40_v1"]),
        "SingleMu_166346-166346_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu40_v2"]),
        "SingleMu_166374-167043_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu40_v1"]),
        "SingleMu_167078-167913_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu40_v3"]),
        "SingleMu_170722-172619_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu40_v5"]),
        "SingleMu_172620-173198_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu40_v5"]),
        "SingleMu_166161-173198_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu40_v1", "HLT_Mu40_v2", "HLT_Mu40_v3", "HLT_Mu40_v5"], triggerThrow=False),

        "SingleMu_173236-173692_2011A_Nov08": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v1"]),

        "SingleMu_173693-177452_2011B_Nov19": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v1"]),
        "SingleMu_177453-178380_2011B_Nov19": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v1"]),
        "SingleMu_178411-179889_2011B_Nov19": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v4"]),
        "SingleMu_179942-180371_2011B_Nov19": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v5"]),
        "SingleMu_175832-180252_2011B_Nov19": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v1", "HLT_Mu40_eta2p1_v4", "HLT_Mu40_eta2p1_v5"], triggerThrow=False),

        # MC, triggered with mcTrigger
        "WJets_TuneZ2_Fall11":               TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": -1}),
        "W1Jets_TuneZ2_Fall11":              TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": 1}),
        "W2Jets_TuneZ2_Fall11":              TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": 2}),
        "W3Jets_TuneZ2_Fall11":           TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": 3}),
        "W3Jets_TuneZ2_v2_Fall11":           TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": 3}),
        "W4Jets_TuneZ2_Fall11":              TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": 4}),
        "TTJets_TuneZ2_Fall11":              TaskDefMC(),
        "DYJetsToLL_M50_TuneZ2_Fall11":      TaskDefMC(),
        "T_t-channel_TuneZ2_Fall11":         TaskDefMC(),
        "Tbar_t-channel_TuneZ2_Fall11":      TaskDefMC(),
        "T_tW-channel_TuneZ2_Fall11":        TaskDefMC(),
        "Tbar_tW-channel_TuneZ2_Fall11":     TaskDefMC(),
        "T_s-channel_TuneZ2_Fall11":         TaskDefMC(),
        "Tbar_s-channel_TuneZ2_Fall11":      TaskDefMC(),
        "WW_TuneZ2_Fall11":                  TaskDefMC(),
        "WZ_TuneZ2_Fall11":                  TaskDefMC(),
        "ZZ_TuneZ2_Fall11":                  TaskDefMC(),
        "QCD_Pt20_MuEnriched_TuneZ2_Fall11": TaskDefMC(),
        "TTToHplusBWB_M80_Fall11":           TaskDefMC(),
        "TTToHplusBWB_M90_Fall11":           TaskDefMC(),
        "TTToHplusBWB_M100_Fall11":          TaskDefMC(),
        "TTToHplusBWB_M120_Fall11":          TaskDefMC(),
        "TTToHplusBWB_M140_Fall11":          TaskDefMC(),
        "TTToHplusBWB_M150_Fall11":          TaskDefMC(),
        "TTToHplusBWB_M155_Fall11":          TaskDefMC(),
        "TTToHplusBWB_M160_Fall11":          TaskDefMC(),
        "TTToHplusBHminusB_M80_Fall11":      TaskDefMC(),
        "TTToHplusBHminusB_M90_Fall11":      TaskDefMC(),
        "TTToHplusBHminusB_M100_Fall11":     TaskDefMC(),
        "TTToHplusBHminusB_M120_Fall11":     TaskDefMC(),
        "TTToHplusBHminusB_M140_Fall11":     TaskDefMC(),
        "TTToHplusBHminusB_M150_Fall11":     TaskDefMC(),
        "TTToHplusBHminusB_M155_Fall11":     TaskDefMC(),
        "TTToHplusBHminusB_M160_Fall11":     TaskDefMC(),
        }

def getDefaultDefinitions_53X():
    mcTrigger = "HLT_Mu40_eta2p1_v9"
    def TaskDefMC(**kwargs):
        return TaskDef(triggerOR=[mcTrigger], **kwargs)

    return {
        "SingleMu_190456-193621_2012A_Jan22": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v9"]),
        "SingleMu_193834-196531_2012B_Jan22": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v9"]),
        "SingleMu_198022-200381_2012C_Jan22": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v10", "HLT_Mu40_eta2p1_v11"], triggerThrow=False),
        "SingleMu_200466-203742_2012C_Jan22": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v11"]),
        "SingleMu_203777-205834_2012D_Jan22": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v11"]),
        "SingleMu_205908-207100_2012D_Jan22": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v11"]),
        "SingleMu_207214-208686_2012D_Jan22": TaskDef(triggerOR=["HLT_Mu40_eta2p1_v11"]),

        # MC, triggered with mcTrigger
        "WJets_TuneZ2star_v1_Summer12":            TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": -1}),
        "WJets_TuneZ2star_v2_Summer12":            TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": -1}),
        "W1Jets_TuneZ2star_Summer12":              TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": 1}),
        "W2Jets_TuneZ2star_Summer12":              TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": 2}),
        "W3Jets_TuneZ2star_Summer12":              TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": 3}),
        "W4Jets_TuneZ2star_Summer12":              TaskDefMC(args={"wjetsWeighting": 1, "wjetBin": 4}),
        "TTJets_TuneZ2star_Summer12":              TaskDefMC(),
        "DYJetsToLL_M50_TuneZ2star_Summer12":      TaskDefMC(),
        "T_t-channel_TuneZ2star_Summer12":         TaskDefMC(),
        "Tbar_t-channel_TuneZ2star_Summer12":      TaskDefMC(),
        "T_tW-channel_TuneZ2star_Summer12":        TaskDefMC(),
        "Tbar_tW-channel_TuneZ2star_Summer12":     TaskDefMC(),
        "T_s-channel_TuneZ2star_Summer12":         TaskDefMC(),
        "Tbar_s-channel_TuneZ2star_Summer12":      TaskDefMC(),
        "WW_TuneZ2star_Summer12":                  TaskDefMC(),
        "WZ_TuneZ2star_Summer12":                  TaskDefMC(),
        "ZZ_TuneZ2star_Summer12":                  TaskDefMC(),
        "QCD_Pt20_MuEnriched_TuneZ2star_Summer12": TaskDefMC(),
        "TTToHplusBWB_M80_ext_Summer12":           TaskDefMC(),
        "TTToHplusBWB_M90_ext_Summer12":           TaskDefMC(),
        "TTToHplusBWB_M100_ext_Summer12":          TaskDefMC(),
        "TTToHplusBWB_M120_ext_Summer12":          TaskDefMC(),
        "TTToHplusBWB_M140_ext_Summer12":          TaskDefMC(),
        "TTToHplusBWB_M150_ext_Summer12":          TaskDefMC(),
        "TTToHplusBWB_M155_ext_Summer12":          TaskDefMC(),
        "TTToHplusBWB_M160_ext_Summer12":          TaskDefMC(),
        "TTToHplusBHminusB_M80_ext_Summer12":      TaskDefMC(),
        "TTToHplusBHminusB_M90_Summer12":          TaskDefMC(),
        "TTToHplusBHminusB_M100_ext_Summer12":     TaskDefMC(),
        "TTToHplusBHminusB_M120_ext_Summer12":     TaskDefMC(),
        "TTToHplusBHminusB_M140_ext_Summer12":     TaskDefMC(),
        "TTToHplusBHminusB_M150_ext_Summer12":     TaskDefMC(),
        "TTToHplusBHminusB_M155_ext_Summer12":     TaskDefMC(),
        "TTToHplusBHminusB_M160_ext_Summer12":     TaskDefMC(),
        }


def addEmbeddingSkim_44X(version, datasets, updateDefinitions):
    defaultDefinitions = getDefaultDefinitions_44X()
    # Specifies the default
    # - number of jobs in skim
    # - number of jobs for those who read the skim (default value for skimAnalysis)
    # - triggers
    # for all relevant Datasets
    # Goal for skim jobs is ~5 hour jobs
    # Goal for skimAnalysis jobs is 30k events/job
    njobs = {
        "SingleMu_160431-163261_2011A_Nov08": TaskDef(njobsIn=100, njobsOut= 1),
        "SingleMu_163270-163869_2011A_Nov08": TaskDef(njobsIn=250, njobsOut= 2),
        "SingleMu_165088-166150_2011A_Nov08": TaskDef(njobsIn=490, njobsOut= 4),

        "SingleMu_166161-166164_2011A_Nov08": TaskDef(njobsIn=  2, njobsOut= 1),
        "SingleMu_166346-166346_2011A_Nov08": TaskDef(njobsIn=  2, njobsOut= 1),
        "SingleMu_166374-167043_2011A_Nov08": TaskDef(njobsIn=300, njobsOut= 6),
        "SingleMu_167078-167913_2011A_Nov08": TaskDef(njobsIn=230, njobsOut= 3),
        "SingleMu_170722-172619_2011A_Nov08": TaskDef(njobsIn=200, njobsOut= 6),
        "SingleMu_172620-173198_2011A_Nov08": TaskDef(njobsIn=230, njobsOut= 6),
        "SingleMu_166161-173198_2011A_Nov08": TaskDef(njobsIn=1700, njobsOut= 25),

        "SingleMu_173236-173692_2011A_Nov08": TaskDef(njobsIn=200, njobsOut= 4),

        "SingleMu_173693-177452_2011B_Nov19": TaskDef(njobsIn=480, njobsOut=16),
        "SingleMu_177453-178380_2011B_Nov19": TaskDef(njobsIn=300, njobsOut=11),
        "SingleMu_178411-179889_2011B_Nov19": TaskDef(njobsIn=300, njobsOut=11),
        "SingleMu_179942-180371_2011B_Nov19": TaskDef(njobsIn= 60, njobsOut= 2),
        "SingleMu_175832-180252_2011B_Nov19": TaskDef(njobsIn=4000, njobsOut=40),

        # MC, triggered with mcTrigger
        "TTJets_TuneZ2_Fall11":              TaskDef(njobsIn=4990, njobsOut=100),
        "WJets_TuneZ2_Fall11":               TaskDef(njobsIn= 700, njobsOut=12),
        "W1Jets_TuneZ2_Fall11":              TaskDef(njobsIn= 700, njobsOut=20),
        "W2Jets_TuneZ2_Fall11":              TaskDef(njobsIn= 750, njobsOut=20),
        "W3Jets_TuneZ2_Fall11":              TaskDef(njobsIn= 750, njobsOut=20),
        "W3Jets_TuneZ2_v2_Fall11":          TaskDef(njobsIn= 750, njobsOut=20),
        "W4Jets_TuneZ2_Fall11":              TaskDef(njobsIn= 750, njobsOut=20),
        "DYJetsToLL_M50_TuneZ2_Fall11":      TaskDef(njobsIn=1700, njobsOut=40),
        "T_t-channel_TuneZ2_Fall11":         TaskDef(njobsIn= 150, njobsOut= 2),
        "Tbar_t-channel_TuneZ2_Fall11":      TaskDef(njobsIn=  70, njobsOut= 1),
        "T_tW-channel_TuneZ2_Fall11":        TaskDef(njobsIn= 100, njobsOut= 2),
        "Tbar_tW-channel_TuneZ2_Fall11":     TaskDef(njobsIn= 100, njobsOut= 2),
        "T_s-channel_TuneZ2_Fall11":         TaskDef(njobsIn=  15, njobsOut= 1),
        "Tbar_s-channel_TuneZ2_Fall11":      TaskDef(njobsIn=  10, njobsOut= 1),
        "WW_TuneZ2_Fall11":                  TaskDef(njobsIn= 200, njobsOut= 4),
        "WZ_TuneZ2_Fall11":                  TaskDef(njobsIn= 200, njobsOut= 4),
        "ZZ_TuneZ2_Fall11":                  TaskDef(njobsIn= 200, njobsOut= 4),
        "QCD_Pt20_MuEnriched_TuneZ2_Fall11": TaskDef(njobsIn= 200, njobsOut= 3),
        "TTToHplusBWB_M80_Fall11":           TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBWB_M90_Fall11":           TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBWB_M100_Fall11":          TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBWB_M120_Fall11":          TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBWB_M140_Fall11":          TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBWB_M150_Fall11":          TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBWB_M155_Fall11":          TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBWB_M160_Fall11":          TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBHminusB_M80_Fall11":      TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBHminusB_M90_Fall11":      TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBHminusB_M100_Fall11":     TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBHminusB_M120_Fall11":     TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBHminusB_M140_Fall11":     TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBHminusB_M150_Fall11":     TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBHminusB_M155_Fall11":     TaskDef(njobsIn=  20, njobsOut= 1),
        "TTToHplusBHminusB_M160_Fall11":     TaskDef(njobsIn=  20, njobsOut= 1),
        }

    workflowName = "tauembedding_skim_"+version

    # Update the default definitions from the argument
    updateTaskDefinitions(defaultDefinitions, njobs, workflowName)
    updateTaskDefinitions(defaultDefinitions, updateDefinitions, workflowName)

    # Add skim Workflow for each dataset
    for datasetName, taskDef in defaultDefinitions.iteritems():
        dataset = datasets.getDataset(datasetName)

        # Construct processing workflow
        wf = constructProcessingWorkflow_44X(dataset, taskDef, sourceWorkflow="AOD", workflowName=workflowName)

        # CRAB configuration lines
        if dataset.isData():
            wf.addCrabLine("CMSSW.total_number_of_lumis = -1")
        else:
            # split by events can only be used for MC and in skim step
            # embedding step is impossible, because the counters are saved
            # in the lumi sections, and will get doubly counted in split
            # by events mode
            wf.addCrabLine("CMSSW.total_number_of_events = -1")

        # Setup the publish name
        if "v44_4" in version:
            name = updatePublishName(dataset, wf.source.getDataForDataset(dataset).getDatasetPath(), workflowName.replace("tau", ""))
        else:
            name = updatePublishName(dataset, wf.source.getDataForDataset(dataset).getDatasetPath(), workflowName)
        wf.addCrabLine("USER.publish_data_name = "+name)

        # Add the skim Workflow to Dataset
        dataset.addWorkflow(wf)

        # If have skim output, define the workflows which depend on it
        if wf.output != None:
            dataset.addWorkflow(Workflow("tauembedding_skimAnalysis_"+version, source=Source("tauembedding_skim_"+version),
                                         triggerOR=taskDef.triggerOR, args=wf.args, output_file="histograms.root", crabLines=["CMSSW.total_number_of_lumis = -1"]))

def addEmbeddingSkim_53X(version, datasets, updateDefinitions):
    defaultDefinitions = getDefaultDefinitions_53X()
    njobs = {
        # Goal for data: 9000 events / job (selection efficiency 11 %)
        "SingleMu_190456-193621_2012A_Jan22": TaskDef(nlumisPerJobIn=1, njobsOut=20),
        "SingleMu_193834-196531_2012B_Jan22": TaskDef(nlumisPerJobIn=1, njobsOut=100),
        "SingleMu_198022-200381_2012C_Jan22": TaskDef(nlumisPerJobIn=1, njobsOut=70),
        "SingleMu_200466-203742_2012C_Jan22": TaskDef(nlumisPerJobIn=1, njobsOut=90),
        "SingleMu_203777-205834_2012D_Jan22": TaskDef(nlumisPerJobIn=1, njobsOut=60),
        "SingleMu_205908-207100_2012D_Jan22": TaskDef(nlumisPerJobIn=1, njobsOut=60),
        "SingleMu_207214-208686_2012D_Jan22": TaskDef(nlumisPerJobIn=1, njobsOut=60),

        # MC, triggered with mcTrigger
        "WJets_TuneZ2star_v1_Summer12":            TaskDef(njobsIn= 250, njobsOut=5),
        "WJets_TuneZ2star_v2_Summer12":            TaskDef(njobsIn= 750, njobsOut=15),
        "W1Jets_TuneZ2star_Summer12":              TaskDef(njobsIn= 400, njobsOut=10),
        "W2Jets_TuneZ2star_Summer12":              TaskDef(njobsIn=2500, njobsOut=40),
        "W3Jets_TuneZ2star_Summer12":              TaskDef(njobsIn=2300, njobsOut=35),
        "W4Jets_TuneZ2star_Summer12":              TaskDef(njobsIn=2200, njobsOut=40),
        "TTJets_TuneZ2star_Summer12":              TaskDef(njobsIn=1700, njobsOut=25),
        "DYJetsToLL_M50_TuneZ2star_Summer12":      TaskDef(njobsIn=2700, njobsOut=40),
        "T_t-channel_TuneZ2star_Summer12":         TaskDef(njobsIn= 200, njobsOut=3),
        "Tbar_t-channel_TuneZ2star_Summer12":      TaskDef(njobsIn= 100, njobsOut=1),
        "T_tW-channel_TuneZ2star_Summer12":        TaskDef(njobsIn= 250, njobsOut=2),
        "Tbar_tW-channel_TuneZ2star_Summer12":     TaskDef(njobsIn= 250, njobsOut=2),
        "T_s-channel_TuneZ2star_Summer12":         TaskDef(njobsIn=  15, njobsOut=1),
        "Tbar_s-channel_TuneZ2star_Summer12":      TaskDef(njobsIn=  10, njobsOut=1),
        "WW_TuneZ2star_Summer12":                  TaskDef(njobsIn= 120, njobsOut=7),
        "WZ_TuneZ2star_Summer12":                  TaskDef(njobsIn= 120, njobsOut=7),
        "ZZ_TuneZ2star_Summer12":                  TaskDef(njobsIn= 120, njobsOut=7),
        "QCD_Pt20_MuEnriched_TuneZ2star_Summer12": TaskDef(njobsIn= 300, njobsOut=4),
        "TTToHplusBWB_M80_ext_Summer12":           TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBWB_M90_ext_Summer12":           TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBWB_M100_ext_Summer12":          TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBWB_M120_ext_Summer12":          TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBWB_M140_ext_Summer12":          TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBWB_M150_ext_Summer12":          TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBWB_M155_ext_Summer12":          TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBWB_M160_ext_Summer12":          TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBHminusB_M80_ext_Summer12":      TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBHminusB_M90_Summer12":          TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBHminusB_M100_ext_Summer12":     TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBHminusB_M120_ext_Summer12":     TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBHminusB_M140_ext_Summer12":     TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBHminusB_M150_ext_Summer12":     TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBHminusB_M155_ext_Summer12":     TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        "TTToHplusBHminusB_M160_ext_Summer12":     TaskDef(njobsIn= 150, njobsOut=1), # FIXME: set njobsOut
        }

    workflowName = "tauembedding_skim_"+version

    # Update the default definitions from the argument
    updateTaskDefinitions(defaultDefinitions, njobs, workflowName)
    updateTaskDefinitions(defaultDefinitions, updateDefinitions, workflowName)

    # Add skim Workflow for each dataset
    for datasetName, taskDef in defaultDefinitions.iteritems():
        dataset = datasets.getDataset(datasetName)

        # Construct processing workflow
        wf = constructProcessingWorkflow_53X(dataset, taskDef, sourceWorkflow="AOD", workflowName=workflowName)

        # CRAB configuration lines
        if dataset.isData():
            wf.addCrabLine("CMSSW.total_number_of_lumis = -1")
        else:
            # split by events can only be used for MC and in skim step
            # embedding step is impossible, because the counters are saved
            # in the lumi sections, and will get doubly counted in split
            # by events mode
            wf.addCrabLine("CMSSW.total_number_of_events = -1")

        # Setup the publish name
        name = updatePublishName(dataset, wf.source.getDataForDataset(dataset).getDatasetPath(), workflowName)
        wf.addCrabLine("USER.publish_data_name = "+name)

        # Add the skim Workflow to Dataset
        dataset.addWorkflow(wf)

        # If have skim output, define the workflows which depend on it
        if wf.output != None:
            dataset.addWorkflow(Workflow("tauembedding_skimAnalysis_"+version, source=Source("tauembedding_skim_"+version),
                                         triggerOR=taskDef.triggerOR, args=wf.args, output_file="histograms.root", crabLines=["CMSSW.total_number_of_lumis = -1"]))


def addEmbeddingEmbedding_44X(sourceWorkflow, version, datasets, updateDefinitions):
    defaultDefinitions = getDefaultDefinitions_44X()
    # njobsIn default value is for embedding
    njobs = {
        "SingleMu_160431-163261_2011A_Nov08": TaskDef(njobsIn=100, njobsOut=1),
        "SingleMu_163270-163869_2011A_Nov08": TaskDef(njobsIn= 70, njobsOut=1),
        "SingleMu_165088-166150_2011A_Nov08": TaskDef(njobsIn=100, njobsOut=1),
        "SingleMu_166161-166164_2011A_Nov08": TaskDef(njobsIn=  1, njobsOut=1),
        "SingleMu_166346-166346_2011A_Nov08": TaskDef(njobsIn=  1, njobsOut=1),
        "SingleMu_166374-167043_2011A_Nov08": TaskDef(njobsIn=110, njobsOut=2),
        "SingleMu_167078-167913_2011A_Nov08": TaskDef(njobsIn= 60, njobsOut=1),
        "SingleMu_170722-172619_2011A_Nov08": TaskDef(njobsIn=120, njobsOut=2),
        "SingleMu_172620-173198_2011A_Nov08": TaskDef(njobsIn=120, njobsOut=2),
        "SingleMu_173236-173692_2011A_Nov08": TaskDef(njobsIn= 70, njobsOut=1),
        "SingleMu_173693-177452_2011B_Nov19": TaskDef(njobsIn=300, njobsOut=4),
        "SingleMu_177453-178380_2011B_Nov19": TaskDef(njobsIn=210, njobsOut=3),
        "SingleMu_178411-179889_2011B_Nov19": TaskDef(njobsIn=200, njobsOut=3),
        "SingleMu_179942-180371_2011B_Nov19": TaskDef(njobsIn= 35, njobsOut=1),

        # MC, triggered with mcTrigger
        "TTJets_TuneZ2_Fall11":              TaskDef(njobsIn=2490, njobsOut=60), # embedding: 5 s per event
        "WJets_TuneZ2_Fall11":               TaskDef(njobsIn= 200, njobsOut= 5),
#        "W3Jets_TuneZ2_Fall11":              TaskDef(njobsIn= 140, njobsOut=42),
        "DYJetsToLL_M50_TuneZ2_Fall11":      TaskDef(njobsIn= 500, njobsOut=12),
        "T_t-channel_TuneZ2_Fall11":         TaskDef(njobsIn=  40, njobsOut= 1),
        "Tbar_t-channel_TuneZ2_Fall11":      TaskDef(njobsIn=  20, njobsOut= 1),
        "T_tW-channel_TuneZ2_Fall11":        TaskDef(njobsIn=  30, njobsOut= 1),
        "Tbar_tW-channel_TuneZ2_Fall11":     TaskDef(njobsIn=  30, njobsOut= 1),
        "T_s-channel_TuneZ2_Fall11":         TaskDef(njobsIn=   3, njobsOut= 1),
        "Tbar_s-channel_TuneZ2_Fall11":      TaskDef(njobsIn=   2, njobsOut= 1),
        "WW_TuneZ2_Fall11":                  TaskDef(njobsIn=  55, njobsOut= 2),
        "WZ_TuneZ2_Fall11":                  TaskDef(njobsIn=  60, njobsOut= 2),
        "ZZ_TuneZ2_Fall11":                  TaskDef(njobsIn=  50, njobsOut= 1),
        "QCD_Pt20_MuEnriched_TuneZ2_Fall11": TaskDef(njobsIn=  45, njobsOut= 1),
        }
    # Update the default definitions from the argument
    updateTaskDefinitions(defaultDefinitions, njobs)
    updateTaskDefinitions(defaultDefinitions, updateDefinitions)

    path_re = re.compile("_tauembedding_.*")

    # Add embedding Workflow for each dataset
    for datasetName, taskDef in defaultDefinitions.iteritems():
        dataset = datasets.getDataset(datasetName)

        wf = constructProcessingWorkflow_44X(dataset, taskDef, sourceWorkflow=sourceWorkflow, workflowName="tauembedding_embedding_"+version)
        wf.source.lumiMask = None # be agnostic for lumi mask
        wf.setOutputFile("embedded.root")

        # CRAB configuration lines
        wf.addCrabLine("CMSSW.total_number_of_lumis = -1")

        # Override the beamspot only for data
        if dataset.isData():
            wf.args["overrideBeamSpot"] = "1"

        # Setup the publish name
        path = wf.source.getDataForDataset(dataset).getDatasetPath().split("/")
        postfix = ""
        if taskDef.publishPostfix is not None:
            postfix = taskDef.publishPostfix
        name = path_re.sub("_tauembedding_embedding_"+version+postfix, path[2])
        name = name.replace("local-", "")
        wf.addCrabLine("USER.publish_data_name = "+name)

        # Add the skim Workflow to Dataset
        dataset.addWorkflow(wf)

        # If have embedding output, define the workflows which depend on it
        if wf.output != None:
            args = {}
            args.update(wf.args)
            args["tauEmbeddingInput"] = "1"
            # these is needed only for embedding jobs
            for a in ["overrideBeamSpot", "tauDecayMode", "tauMinVisPt"]:
                try:
                    del args[a] 
                except KeyError:
                    pass
            wf_analysis = Workflow("tauembedding_analysis_"+version, source=Source("tauembedding_embedding_"+version),
                                   triggerOR=taskDef.triggerOR, args=args, output_file="histograms.root")
            wf_analysis.addCrabLine("CMSSW.total_number_of_lumis = -1")
            dataset.addWorkflow(wf_analysis)
 
def addEmbeddingSkim_v44_4_2(datasets):
    definitions = {
        # 113 files, min 86 MB, max, 111 MB
        "SingleMu_160431-163261_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_160431_tauembedding_skim_v44_4_2-d9eec32ec495d473673489f496724114/USER"),
        # 125 files, min 69 MB, max 235 MB
        "SingleMu_163270-163869_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_163270_tauembedding_skim_v44_4_2-f2f022e2f1824ce67ef0c386cf58329f/USER", njobsIn=140),
        # 373 files, min 3.5 MB, max 134 MB
        "SingleMu_165088-166150_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_165088_tauembedding_skim_v44_4_2-70b45374fbabe2248133e36a2cbe01e2/USER"),
        # 2 files, min 9.8 MB, max 138 MB
        "SingleMu_166161-166164_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_166161_tauembedding_skim_v44_4_2-0a9e46f53bac3a3199fc5d08e63772d3/USER"),
        # 2 files, min 116 MB, max 142 MB
        "SingleMu_166346-166346_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_166346_tauembedding_skim_v44_4_2-91fdceb3c5af341c67f22e9c64363c60/USER"),
        # 260 files, min 32 MB, max 142 MB
        "SingleMu_166374-167043_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_166374_tauembedding_skim_v44_4_2-0a9e46f53bac3a3199fc5d08e63772d3/USER"),
        # 191 files, min 13 MB, max 113 MB
        "SingleMu_167078-167913_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_167078_tauembedding_skim_v44_4_2-c4a8eba430793887b650ad5b083c7ed7/USER"),
        # 177 files, min 60 MB, max 272 MB
        "SingleMu_170722-172619_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_170722_tauembedding_skim_v44_4_2-859a149425fa4c077a5b33666a7b993a/USER"),
        # 203 files, min 44 MB, max 217 MB
        "SingleMu_172620-173198_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_172620_tauembedding_skim_v44_4_2-859a149425fa4c077a5b33666a7b993a/USER"),
        # 108 files, min 25 MB, max 246 MB
        "SingleMu_173236-173692_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_173236_tauembedding_skim_v44_4_2-e4bb75f1eb1d67eb9cca9dc53de3dd14/USER", njobsIn=120),
        # 417 files, min 18 MB, max 405 MB
        "SingleMu_173693-177452_2011B_Nov19": TaskDef("/SingleMu/local-Run2011B_19Nov2011_v1_AOD_173693_tauembedding_skim_v44_4_2-e4bb75f1eb1d67eb9cca9dc53de3dd14/USER"),
        # 263 files, min 21 MB, max 449 MB
        "SingleMu_177453-178380_2011B_Nov19": TaskDef("/SingleMu/local-Run2011B_19Nov2011_v1_AOD_177453_tauembedding_skim_v44_4_2-e4bb75f1eb1d67eb9cca9dc53de3dd14/USER"),
        # 259 files, min 104 MB, max 421 MB
        "SingleMu_178411-179889_2011B_Nov19": TaskDef("/SingleMu/local-Run2011B_19Nov2011_v1_AOD_178411_tauembedding_skim_v44_4_2-ce29afbda8ffd97a0d906ada6ad5e907/USER"),
        # 53 files, min 134 MB, max 291 MB
        "SingleMu_179942-180371_2011B_Nov19": TaskDef("/SingleMu/local-Run2011B_19Nov2011_v1_AOD_179942_tauembedding_skim_v44_4_2-2fa526d919657b36eff79f06bc501a9c/USER"),
        # 4961 files, min 3.3 MB, max 232 MB
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER"),
        # 998 files, min 19 MB, max 104 MB
        "WJets_TuneZ2_Fall11":                TaskDef("/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=990),
        # 989 files, min 74 MB, max 252 MB
        "DYJetsToLL_M50_TuneZ2_Fall11":       TaskDef("/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=990),
        # 492 files, min 9.8 MB, max 45 MB
        "T_t-channel_TuneZ2_Fall11":          TaskDef("/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=300),
        # 162 files, min 6.1 MB, max 63 MB
        "Tbar_t-channel_TuneZ2_Fall11":       TaskDef("/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=70),
        # 91 files, min 5.1 MB, max 156 MB
        "T_tW-channel_TuneZ2_Fall11":         TaskDef("/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=90),
        # 92 files, min 11 MB, max 149 MB
        "Tbar_tW-channel_TuneZ2_Fall11":      TaskDef("/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=90),
        # 51 files, min 4.9 MB, max 40 MB
        "T_s-channel_TuneZ2_Fall11":          TaskDef("/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=25),
        # 12 files, min 20 MB, max 83 MB
        "Tbar_s-channel_TuneZ2_Fall11":       TaskDef("/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=5),
        # 201 files, min 3.9 MB, max 136 MB
        "WW_TuneZ2_Fall11":                   TaskDef("/WW_TuneZ2_7TeV_pythia6_tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=150),
        # 201 files, min 26 MB, max 122 MB
        "WZ_TuneZ2_Fall11":                   TaskDef("/WZ_TuneZ2_7TeV_pythia6_tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=150),
        # 352 files, min 19 MB, max 77 MB
        "ZZ_TuneZ2_Fall11":                   TaskDef("/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER"),
        # 493 files, min 14 MB, max 54 MB
        "QCD_Pt20_MuEnriched_TuneZ2_Fall11":  TaskDef("/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_4_2-80448358a193f69c52fb3eaa57e02bff/USER", njobsIn=490),
        }
    addEmbeddingSkim_44X("v44_4_2", datasets, definitions)

def addEmbeddingEmbedding_v44_4_2(datasets):
    skimVersion = "tauembedding_skim_v44_4_2"

    def addEmbedding(version, definitions):
        addEmbeddingEmbedding_44X(skimVersion, version, datasets, definitions)

    addEmbedding("v44_4_2_muiso0", {
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_embedding_v44_4_2_muiso0-50da2d6a5b0c9c8a2f96f633ada0c1c6/USER"),
        })
    addEmbedding("v44_4_2_muiso1", {
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_embedding_v44_4_2_muiso1-50da2d6a5b0c9c8a2f96f633ada0c1c6/USER"),
        })

    addEmbedding("v44_4_2_seed0", {
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_embedding_v44_4_2_seed0-2dedf078d8faded30b2dddce6fe8cdec/USER"),
        })
    
    addEmbedding("v44_4_2_seed1", {
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_embedding_v44_4_2_seed1-2dedf078d8faded30b2dddce6fe8cdec/USER"),
        })

def addEmbeddingGenTauSkim_v44_5(datasets):
    definitions = {
        # 6662721 events, 4002 jobs
        # User mean 2775.2, min 916.6, max 4447.0
        # Mean 130.1 MB, min 44.1 MB, max 141.2 MB
        # The files in this dataset are actually inside pattuple directory
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_gentauskim_v44_5-9ecb3a23e436fc2ffd8a803eac2a3a15/USER"),
        }
    addEmbeddingGenTauSkim_44X("v44_5", datasets, definitions)

def addEmbeddingSkim_v44_5(datasets):
    # Expecting 33 % file size increase w.r.t. v44_2 (237.7/186.5=27%
    # for event size, 107/102=5% for number of events)
    definitions = {
        # 8055624 events, 4989 jobs
        # User mean 3192.2, min 16.7, max 6071.3
        # Mean 339.4 MB, min 4.1 MB, max 376.1 MB
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5-a9adb1d2c9d25e1e9802345c8c130cf6/USER", args={"triggerMC": 0}), # disable trigger in skim jobs for TTJets
        }
    addEmbeddingSkim_44X("v44_5", datasets, definitions)

def addEmbeddingEmbedding_v44_5(datasets):
    skimVersion = "tauembedding_skim_v44_5"

    def addEmbedding(version, definitions):
        addEmbeddingEmbedding_44X(skimVersion, version, datasets, definitions)

    addEmbedding("v44_5_notrg2", {
        # 5458414 events, 2066 jobs
        # User mean 12701.9, min 7665.2, max 17593.6
        # Mean 414.2 MB, min 321.6 MB, max 543.0 MB
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_embedding_v44_5_notrg2-398a1e0c60b666df745f371020ad555a/USER", args={"triggerMC": 0}),
        })


def addEmbeddingSkim_v44_5_1(datasets):
    # Not doing chi2<10 ID cut in the skim job (for possibility of TuneP)
#        # 89420 events, 113 jobs
#        # User mean 2432.4, min 470.6, max 7155.4
#        # Mean 67.3 MB, min 9.1 MB, max 151.4 MB
#        "SingleMu_160431-163261_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_160431_163261_tauembedding_skim_v44_5_1-e6e5b12b8c9c249a4985635b30210544/USER"),
#        # 258958 events, 206 jobs
#        # User mean 3153.1, min 839.5, max 10264.4
#        # Mean 101.6 MB, min 18.1 MB, max 217.9 MB
#        "SingleMu_163270-163869_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_163270_163869_tauembedding_skim_v44_5_1-b86f83b25a65adafc510540db2b385d1/USER"),
#        # 386851 events, 373 jobs
#        # User mean 3321.1, min 60.2, max 8545.0
#        # Mean 87.6 MB, min 3.7 MB, max 181.6 MB
#        "SingleMu_165088-166150_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_165088_166150_tauembedding_skim_v44_5_1-856f19da8dd50633d34c9707f7e229d2/USER"),
#        # 1539955 events, 953 jobs
#        # User mean 4592.0, min 1223.0, max 11895.6
#        # Mean 134.1 MB, min 39.8 MB, max 287.3 MB
#        "SingleMu_166161-173198_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_166161_173198_tauembedding_skim_v44_5_1-8a2a027b46e0bf70954ec35b9e70b0b2/USER"),
#        # 245801 events, 162 jobs
#        # User mean 4544.6, min 1874.6, max 10131.6
#        # Mean 140.6 MB, min 55.2 MB, max 226.0 MB
#        "SingleMu_173236-173692_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_173236_173692_tauembedding_skim_v44_5_1-27084e7b9fd52190a9502d45e7513543/USER"),
#        # 2842141 events, 1924 jobs
#        # User mean 4694.9, min 1334.5, max 9042.2
#        # Mean 163.3 MB, min 47.4 MB, max 318.7 MB
#        "SingleMu_175832-180252_2011B_Nov19": TaskDef("/SingleMu/local-Run2011B_19Nov2011_v1_AOD_175832_180252_tauembedding_skim_v44_5_1-39b55c7eb0c557d9c273b43a0be4a2ff/USER"),
    # Bugfix for calculating the "standard" MET filters for data (v44_5_2)
    definitions = {
        # 89420 events, 95 jobs
        # User mean 2858.8, min 633.8, max 6299.9
        # Mean 78.8 MB, min 23.5 MB, max 180.3 MB
        "SingleMu_160431-163261_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_160431_163261_tauembedding_skim_v44_5_2-5da894984d43d7affb331f2626592c2b/USER"),
        # 258958 events, 206 jobs
        # User mean 3522.0, min 595.6, max 7687.9
        # Mean 101.6 MB, min 18.1 MB, max 217.9 MB
        "SingleMu_163270-163869_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_163270_163869_tauembedding_skim_v44_5_2-5dcc6e23548084ce0af9542cb166052d/USER"),
        # 386851 events, 373 jobs
        # User mean 2765.5, min 45.5, max 7547.7
        # Mean 87.6 MB, min 3.7 MB, max 181.6 MB
        "SingleMu_165088-166150_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_165088_166150_tauembedding_skim_v44_5_2-62bc199df7c6de3914ddfd4e7a4dedf1/USER"),
        # 1539955 events, 1295 jobs
        # User mean 3107.2, min 374.8, max 9280.2
        # Mean 100.6 MB, min 15.2 MB, max 219.0 MB
        "SingleMu_166161-173198_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_166161_173198_tauembedding_skim_v44_5_2-9279e439e9124cd35585c1246e432cfd/USER"),
        # 245801 events, 162 jobs
        # User mean 4546.0, min 1653.9, max 10561.0
        # Mean 140.6 MB, min 55.2 MB, max 226.0 MB
        "SingleMu_173236-173692_2011A_Nov08": TaskDef("/SingleMu/local-Run2011A_08Nov2011_v1_AOD_173236_173692_tauembedding_skim_v44_5_2-883c99ed2824db47683a70a86f21fc2e/USER"),
        # 2842141 events, 2638 jobs
        # User mean 3423.2, min 1095.3, max 7534.9
        # Mean 120.9 MB, min 47.4 MB, max 253.1 MB
        "SingleMu_175832-180252_2011B_Nov19": TaskDef("/SingleMu/local-Run2011B_19Nov2011_v1_AOD_175832_180252_tauembedding_skim_v44_5_2-ff6b0c8176db587721dc02d5ab403d42/USER"),

        # 8183108 events, 4994 jobs
        # User mean 3079.8, min 14.8, max 6305.4
        # Mean 346.4 MB, min 4.1 MB, max 384.5 MB
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-29c7b20153e31b2c86faa8316de20ff4/USER", args={"triggerMC": 0}), # disable trigger in skim jobs for TTJets
        # 624557 events, 707 jobs
        # User mean 6162.0, min 487.4, max 16838.2
        # Mean 161.1 MB, min 18.8 MB, max 181.9 MB
        "WJets_TuneZ2_Fall11":                TaskDef("/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 902288 events, 705 jobs
        # User mean 7644.3, min 104.5, max 20925.5
        # Mean 236.6 MB, min 6.3 MB, max 261.4 MB
        "W1Jets_TuneZ2_Fall11":               TaskDef("/W1Jet_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 1295031 events, 752 jobs
        # User mean 5582.2, min 981.8, max 6633.3
        # Mean 290.2 MB, min 46.2 MB, max 318.4 MB
        "W2Jets_TuneZ2_Fall11":               TaskDef("/W2Jets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 607948 events, 751 jobs
        # User mean 1551.8, min 193.7, max 2137.6
        # Mean 148.9 MB, min 9.5 MB, max 165.7 MB
        "W3Jets_TuneZ2_v2_Fall11":            TaskDef("/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v2_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 1374644 events, 753 jobs
        # User mean 5269.6, min 355.8, max 24077.4
        # Mean 351.7 MB, min 18.6 MB, max 378.0 MB
        "W4Jets_TuneZ2_Fall11":               TaskDef("/W4Jets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 1601361 events, 1703 jobs
        # User mean 4622.9, min 495.2, max 28479.3
        # Mean 167.9 MB, min 14.6 MB, max 188.5 MB
        "DYJetsToLL_M50_TuneZ2_Fall11":       TaskDef("/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 106575 events, 151 jobs
        # User mean 1820.1, min 1053.7, max 2366.1
        # Mean 132.4 MB, min 76.3 MB, max 147.6 MB
        "T_t-channel_TuneZ2_Fall11":          TaskDef("/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 57734 events, 73 jobs
        # User mean 1879.4, min 235.2, max 3751.7
        # Mean 146.8 MB, min 14.5 MB, max 165.3 MB
        "Tbar_t-channel_TuneZ2_Fall11":       TaskDef("/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 81230 events, 102 jobs
        # User mean 1796.3, min 332.5, max 5872.3
        # Mean 157.1 MB, min 31.7 MB, max 177.9 MB
        "T_tW-channel_TuneZ2_Fall11":         TaskDef("/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 80555 events, 101 jobs
        # User mean 1981.5, min 601.9, max 6133.1
        # Mean 157.2 MB, min 65.3 MB, max 175.5 MB
        "Tbar_tW-channel_TuneZ2_Fall11":      TaskDef("/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 9585 events, 16 jobs
        # User mean 1320.4, min 664.1, max 1720.0
        # Mean 116.2 MB, min 60.0 MB, max 133.6 MB
        "T_s-channel_TuneZ2_Fall11":          TaskDef("/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 5157 events, 12 jobs
        # User mean 916.4, min 317.6, max 1261.1
        # Mean 84.0 MB, min 24.8 MB, max 103.6 MB
        "Tbar_s-channel_TuneZ2_Fall11":       TaskDef("/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 171457 events, 201 jobs
        # User mean 2300.6, min 147.8, max 3415.4
        # Mean 151.7 MB, min 4.9 MB, max 172.0 MB
        "WW_TuneZ2_Fall11":                   TaskDef("/WW_TuneZ2_7TeV_pythia6_tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 155345 events, 201 jobs
        # User mean 1877.3, min 639.6, max 2842.6
        # Mean 140.1 MB, min 30.9 MB, max 154.0 MB
        "WZ_TuneZ2_Fall11":                   TaskDef("/WZ_TuneZ2_7TeV_pythia6_tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 148668 events, 202 jobs
        # User mean 1551.9, min 545.3, max 2315.1
        # Mean 136.3 MB, min 43.4 MB, max 152.0 MB
        "ZZ_TuneZ2_Fall11":                   TaskDef("/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),
        # 163078 events, 203 jobs
        # User mean 2669.4, min 283.0, max 4216.4
        # Mean 159.9 MB, min 18.6 MB, max 174.4 MB
        "QCD_Pt20_MuEnriched_TuneZ2_Fall11":  TaskDef("/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_skim_v44_5_1-ce9fe259c534bb39b2e660365c1c0aaf/USER"),

        "TTToHplusBWB_M80_Fall11":          TaskDef(""),
        "TTToHplusBWB_M90_Fall11":          TaskDef(""),
        "TTToHplusBWB_M100_Fall11":         TaskDef(""),
        "TTToHplusBWB_M120_Fall11":         TaskDef(""),
        "TTToHplusBWB_M140_Fall11":         TaskDef(""),
        "TTToHplusBWB_M150_Fall11":         TaskDef(""),
        "TTToHplusBWB_M155_Fall11":         TaskDef(""),
        "TTToHplusBWB_M160_Fall11":         TaskDef(""),

        "TTToHplusBHminusB_M80_Fall11":     TaskDef(""),
        "TTToHplusBHminusB_M90_Fall11":     TaskDef(""),
        "TTToHplusBHminusB_M100_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M120_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M140_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M150_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M155_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M160_Fall11":    TaskDef(""),
        }
    addEmbeddingSkim_44X("v44_5_1", datasets, definitions)


def addEmbeddingEmbedding_v44_5_1(datasets):
    skimVersion = "tauembedding_skim_v44_5_1"

    def addEmbedding(version, definitions):
        addEmbeddingEmbedding_44X(skimVersion, version, datasets, definitions)

    addEmbedding("v44_5_1_notrg2", {
        # 5296176 events, 2254 jobs
        # User mean 11404.4, min 9613.9, max 16618.0
        # Mean 369.0 MB, min 314.1 MB, max 530.3 MB
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_embedding_v44_5_1_notrg2-b4444849cbd68cba8058d20690fa09f4/USER", args={"triggerMC": 0}),
        })

    addEmbedding("v44_5_1", {
        # 4741432 events, 2254 jobs
        # User mean 10503.6, min 8873.8, max 15635.7
        # Mean 341.6 MB, min 290.8 MB, max 498.7 MB
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_embedding_v44_5_1-7e75763c5b9284e37debaf454cd27985/USER", args={"triggerMC": 0, "triggerMCInAnalysis": 1}),
        })

    addEmbedding("v44_5_1_tauhad", {
        # 4883090 events, 2254 jobs
        # User mean 11790.0, min 9915.0, max 17016.8
        # Mean 344.0 MB, min 291.0 MB, max 499.5 MB
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_embedding_v44_5_1_tauhad-d57ea742826c3abce18a6ceed0c3bca3/USER", args={"triggerMC": 0, "triggerMCInAnalysis": 1, "tauDecayMode": 230}),
        })

    # addEmbedding("v44_5_1_tauhad_vispt10", {
    #     "TTJets_TuneZ2_Fall11":               TaskDef("", args={"triggerMC": 1, "tauDecayMode": 230, "tauMinVisPt": 10}),
    #     })

    # addEmbedding("v44_5_1_tauhad_vispt20", {
    #     "TTJets_TuneZ2_Fall11":               TaskDef("", args={"triggerMC": 1, "tauDecayMode": 230, "tauMinVisPt": 20}),
    #     })

    addEmbedding("v44_5_1_tauhad_vispt30", {
        # 4933858 events, 2254 jobs
        # User mean 13113.8, min 10899.2, max 19112.7
        # Mean 345.6 MB, min 291.8 MB, max 500.5 MB
        "TTJets_TuneZ2_Fall11":               TaskDef("/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_embedding_v44_5_1_tauhad_vispt30_b-d57ea742826c3abce18a6ceed0c3bca3/USER", args={"triggerMC": 0, "triggerMCInAnalysis": 1, "tauDecayMode": 230, "tauMinVisPt": 30}, publishPostfix="_b"),
        })


def addEmbeddingSkim_v53_3(datasets):
    definitions = {
        # 1324642 events, 817 jobs
        # User mean 4515.8, min 57.0, max 14153.3
        # Mean 225.4 MB, min 2.1 MB, max 341.5 MB
        "SingleMu_190456-193621_2012A_Jan22":      TaskDef("/SingleMu/local-Run2012A_22Jan2013_v1_AOD_190456_193621_tauembedding_skim_v53_3-b1df9b924b2706d50674792311c27d18/USER"),
        # 6409898 events, 4204 jobs
        # User mean 4189.4, min 380.1, max 13122.2
        # Mean 227.0 MB, min 30.8 MB, max 286.6 MB
        "SingleMu_193834-196531_2012B_Jan22":      TaskDef("/SingleMu/local-Run2012B_22Jan2013_v1_AOD_193834_196531_tauembedding_skim_v53_3-b1df9b924b2706d50674792311c27d18/USER"),
        # 4485614 events, 3478 jobs
        # User mean 3986.9, min 28.2, max 11804.8
        # Mean 187.6 MB, min 2.0 MB, max 276.1 MB
        "SingleMu_198022-200381_2012C_Jan22":      TaskDef("/SingleMu/local-Run2012C_22Jan2013_v1_AOD_198022_200381_tauembedding_skim_v53_3-544fbabef383a6c007d0a78670accc6a/USER"),
        # 5372985 events, 4277 jobs
        # User mean 3117.4, min 25.9, max 11646.1
        # Mean 197.3 MB, min 1.9 MB, max 290.7 MB
        "SingleMu_200466-203742_2012C_Jan22":      TaskDef("/SingleMu/local-Run2012C_22Jan2013_v1_AOD_200466_203742_tauembedding_skim_v53_3-bf6e26bfda4583e5a02e30bcb8e788ff/USER"),
        # 3643631 events, 3593 jobs
        # User mean 2161.4, min 58.3, max 6779.1
        # Mean 162.8 MB, min 2.2 MB, max 359.8 MB
        "SingleMu_203777-205834_2012D_Jan22":      TaskDef("/SingleMu/local-Run2012D_22Jan2013_v1_AOD_203777_205834_tauembedding_skim_v53_3-bf6e26bfda4583e5a02e30bcb8e788ff/USER"),
        # 3361174 events, 3613 jobs
        # User mean 2089.3, min 25.7, max 8453.0
        # Mean 156.4 MB, min 1.8 MB, max 281.4 MB
        "SingleMu_205908-207100_2012D_Jan22":      TaskDef("/SingleMu/local-Run2012D_22Jan2013_v1_AOD_205908_207100_tauembedding_skim_v53_3-bf6e26bfda4583e5a02e30bcb8e788ff/USER"),
        # 3697195 events, 3888 jobs
        # User mean 2305.3, min 26.9, max 6659.2
        # Mean 156.9 MB, min 2.0 MB, max 278.9 MB
        "SingleMu_207214-208686_2012D_Jan22":      TaskDef("/SingleMu/local-Run2012D_22Jan2013_v1_AOD_207214_208686_tauembedding_skim_v53_3-bf6e26bfda4583e5a02e30bcb8e788ff/USER"),

        # 204429 events, 254 jobs
        # User mean 3740.7, min 310.5, max 5605.4
        # Mean 202.1 MB, min 15.6 MB, max 227.9 MB
        "WJets_TuneZ2star_v1_Summer12":            TaskDef("/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 641634 events, 758 jobs
        # User mean 3969.3, min 130.1, max 7533.6
        # Mean 212.5 MB, min 7.9 MB, max 237.0 MB
        "WJets_TuneZ2star_v2_Summer12":            TaskDef("/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/local-Summer12_DR53X_PU_S10_START53_V7A_v2_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 481424 events, 405 jobs
        # User mean 4812.2, min 117.9, max 7553.7
        # Mean 293.5 MB, min 7.8 MB, max 322.1 MB
        "W1Jets_TuneZ2star_Summer12":              TaskDef("/W1JetsToLNu_TuneZ2Star_8TeV-madgraph/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 1892218 events, 2506 jobs
        # User mean 1763.5, min 317.8, max 2734.5
        # Mean 184.1 MB, min 24.1 MB, max 205.1 MB
        "W2Jets_TuneZ2star_Summer12":              TaskDef("/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 1253103 events, 2305 jobs
        # User mean 1178.5, min 403.8, max 1744.8
        # Mean 140.7 MB, min 38.5 MB, max 159.0 MB
        "W3Jets_TuneZ2star_Summer12":              TaskDef("/W3JetsToLNu_TuneZ2Star_8TeV-madgraph/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 1346368 events, 2201 jobs
        # User mean 1311.4, min 448.9, max 1942.2
        # Mean 168.2 MB, min 69.8 MB, max 189.8 MB
        "W4Jets_TuneZ2star_Summer12":              TaskDef("/W4JetsToLNu_TuneZ2Star_8TeV-madgraph/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 959991 events, 1702 jobs
        # User mean 1246.4, min 257.4, max 2356.3
        # Mean 165.9 MB, min 32.3 MB, max 188.3 MB
        "TTJets_TuneZ2star_Summer12":              TaskDef("/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-9a24e6fe0421ec76a55ad5183bef176f/USER", args={"triggerMC": 0}), # disable trigger in skim jobs for TTJets
        # 1739021 events, 2704 jobs
        # User mean 1997.3, min 473.7, max 3827.5
        # Mean 162.0 MB, min 36.4 MB, max 186.1 MB
        "DYJetsToLL_M50_TuneZ2star_Summer12":      TaskDef("/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 110024 events, 202 jobs
        # User mean 1394.4, min 321.1, max 2378.6
        # Mean 144.1 MB, min 33.3 MB, max 160.3 MB
        "T_t-channel_TuneZ2star_Summer12":         TaskDef("/T_t-channel_TuneZ2star_8TeV-powheg-tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 60692 events, 101 jobs
        # User mean 1436.4, min 97.8, max 2590.2
        # Mean 156.5 MB, min 2.0 MB, max 172.6 MB
        "Tbar_t-channel_TuneZ2star_Summer12":      TaskDef("/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 49646 events, 251 jobs
        # User mean 447.5, min 118.8, max 713.1
        # Mean 59.6 MB, min 10.1 MB, max 68.6 MB
        "T_tW-channel_TuneZ2star_Summer12":        TaskDef("/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 49385 events, 252 jobs
        # User mean 443.2, min 130.3, max 731.1
        # Mean 59.3 MB, min 11.4 MB, max 68.0 MB
        "Tbar_tW-channel_TuneZ2star_Summer12":     TaskDef("/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 5172 events, 11 jobs
        # User mean 1042.7, min 26.5, max 1417.7
        # Mean 125.1 MB, min 2.0 MB, max 145.8 MB
        "T_s-channel_TuneZ2star_Summer12":         TaskDef("/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 5172 events, 11 jobs
        # User mean 1042.7, min 26.5, max 1417.7
        # Mean 125.1 MB, min 2.0 MB, max 145.8 MB
        "Tbar_s-channel_TuneZ2star_Summer12":      TaskDef("/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 444533 events, 122 jobs
        # User mean 13294.8, min 1379.8, max 16795.1
        # Mean 868.2 MB, min 79.1 MB, max 925.0 MB
        "WW_TuneZ2star_Summer12":                  TaskDef("/WW_TuneZ2star_8TeV_pythia6_tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 381239 events, 123 jobs
        # User mean 9276.8, min 826.0, max 13801.0
        # Mean 762.4 MB, min 92.0 MB, max 811.4 MB
        "WZ_TuneZ2star_Summer12":                  TaskDef("/WZ_TuneZ2star_8TeV_pythia6_tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 357677 events, 125 jobs
        # User mean 9069.8, min 157.8, max 11537.0
        # Mean 716.3 MB, min 12.1 MB, max 789.2 MB
        "ZZ_TuneZ2star_Summer12":                  TaskDef("/ZZ_TuneZ2star_8TeV_pythia6_tauola/local-Summer12_DR53X_PU_S10_START53_V7A_v1_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        # 155596 events, 303 jobs
        # User mean 2023.1, min 308.4, max 3119.9
        # Mean 140.8 MB, min 13.9 MB, max 159.0 MB
        "QCD_Pt20_MuEnriched_TuneZ2star_Summer12": TaskDef("/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/local-Summer12_DR53X_PU_S10_START53_V7A_v3_AODSIM_tauembedding_skim_v53_3-8832fd6cee6cb431f37d162280985223/USER"),
        "TTToHplusBWB_M80_ext_Summer12":           TaskDef(""),
        "TTToHplusBWB_M90_ext_Summer12":           TaskDef(""),
        "TTToHplusBWB_M100_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M120_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M140_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M150_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M155_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M160_ext_Summer12":          TaskDef(""),
        "TTToHplusBHminusB_M80_ext_Summer12":      TaskDef(""),
        "TTToHplusBHminusB_M90_Summer12":          TaskDef(""),
        "TTToHplusBHminusB_M100_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M120_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M140_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M150_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M155_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M160_ext_Summer12":     TaskDef(""),
        }
    addEmbeddingSkim_53X("v53_3", datasets, definitions)


def addEmbedding_SKELETON(datasets):
    definitions = {
        "SingleMu_160431-163261_2011A_Nov08": TaskDef(""),
        "SingleMu_163270-163869_2011A_Nov08": TaskDef(""),
        "SingleMu_165088-166150_2011A_Nov08": TaskDef(""),
        "SingleMu_166161-173198_2011A_Nov08": TaskDef(""),
        "SingleMu_173236-173692_2011A_Nov08": TaskDef(""),
        "SingleMu_175832-180252_2011B_Nov19": TaskDef(""),

        "TTJets_TuneZ2_Fall11":               TaskDef(""),
        "WJets_TuneZ2_Fall11":                TaskDef(""),
        "W1Jets_TuneZ2_Fall11":               TaskDef(""),
        "W2Jets_TuneZ2_Fall11":               TaskDef(""),
        "W3Jets_TuneZ2_v2_Fall11":            TaskDef(""),
        "W4Jets_TuneZ2_Fall11":               TaskDef(""),
        "DYJetsToLL_M50_TuneZ2_Fall11":       TaskDef(""),
        "T_t-channel_TuneZ2_Fall11":          TaskDef(""),
        "Tbar_t-channel_TuneZ2_Fall11":       TaskDef(""),
        "T_tW-channel_TuneZ2_Fall11":         TaskDef(""),
        "Tbar_tW-channel_TuneZ2_Fall11":      TaskDef(""),
        "T_s-channel_TuneZ2_Fall11":          TaskDef(""),
        "Tbar_s-channel_TuneZ2_Fall11":       TaskDef(""),
        "WW_TuneZ2_Fall11":                   TaskDef(""),
        "WZ_TuneZ2_Fall11":                   TaskDef(""),
        "ZZ_TuneZ2_Fall11":                   TaskDef(""),
        "QCD_Pt20_MuEnriched_TuneZ2_Fall11":  TaskDef(""),

        "TTToHplusBWB_M80_Fall11":          TaskDef(""),
        "TTToHplusBWB_M90_Fall11":          TaskDef(""),
        "TTToHplusBWB_M100_Fall11":         TaskDef(""),
        "TTToHplusBWB_M120_Fall11":         TaskDef(""),
        "TTToHplusBWB_M140_Fall11":         TaskDef(""),
        "TTToHplusBWB_M150_Fall11":         TaskDef(""),
        "TTToHplusBWB_M155_Fall11":         TaskDef(""),
        "TTToHplusBWB_M160_Fall11":         TaskDef(""),

        "TTToHplusBHminusB_M80_Fall11":     TaskDef(""),
        "TTToHplusBHminusB_M90_Fall11":     TaskDef(""),
        "TTToHplusBHminusB_M100_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M120_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M140_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M150_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M155_Fall11":    TaskDef(""),
        "TTToHplusBHminusB_M160_Fall11":    TaskDef(""),
        }


def addEmbedding_SKELETON_53X(datasets):
    definitions = {
        "SingleMu_190456-193621_2012A_Jan22":      TaskDef(""),
        "SingleMu_193834-196531_2012B_Jan22":      TaskDef(""),
        "SingleMu_198022-200381_2012C_Jan22":      TaskDef(""),
        "SingleMu_200466-203742_2012C_Jan22":      TaskDef(""),
        "SingleMu_203777-205834_2012D_Jan22":      TaskDef(""),
        "SingleMu_205908-207100_2012D_Jan22":      TaskDef(""),
        "SingleMu_207214-208686_2012D_Jan22":      TaskDef(""),

        "WJets_TuneZ2star_v1_Summer12":            TaskDef(""),
        "WJets_TuneZ2star_v2_Summer12":            TaskDef(""),
        "W1Jets_TuneZ2star_Summer12":              TaskDef(""),
        "W2Jets_TuneZ2star_Summer12":              TaskDef(""),
        "W3Jets_TuneZ2star_Summer12":              TaskDef(""),
        "W4Jets_TuneZ2star_Summer12":              TaskDef(""),
        "TTJets_TuneZ2star_Summer12":              TaskDef(""),
        "DYJetsToLL_M50_TuneZ2star_Summer12":      TaskDef(""),
        "T_t-channel_TuneZ2star_Summer12":         TaskDef(""),
        "Tbar_t-channel_TuneZ2star_Summer12":      TaskDef(""),
        "T_tW-channel_TuneZ2star_Summer12":        TaskDef(""),
        "Tbar_tW-channel_TuneZ2star_Summer12":     TaskDef(""),
        "T_s-channel_TuneZ2star_Summer12":         TaskDef(""),
        "Tbar_s-channel_TuneZ2star_Summer12":      TaskDef(""),
        "WW_TuneZ2star_Summer12":                  TaskDef(""),
        "WZ_TuneZ2star_Summer12":                  TaskDef(""),
        "ZZ_TuneZ2star_Summer12":                  TaskDef(""),
        "QCD_Pt20_MuEnriched_TuneZ2star_Summer12": TaskDef(""),
        "TTToHplusBWB_M80_ext_Summer12":           TaskDef(""),
        "TTToHplusBWB_M90_ext_Summer12":           TaskDef(""),
        "TTToHplusBWB_M100_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M120_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M140_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M150_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M155_ext_Summer12":          TaskDef(""),
        "TTToHplusBWB_M160_ext_Summer12":          TaskDef(""),
        "TTToHplusBHminusB_M80_ext_Summer12":      TaskDef(""),
        "TTToHplusBHminusB_M90_Summer12":          TaskDef(""),
        "TTToHplusBHminusB_M100_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M120_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M140_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M150_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M155_ext_Summer12":     TaskDef(""),
        "TTToHplusBHminusB_M160_ext_Summer12":     TaskDef(""),
        }

############################# below are old definitions, which still exist in DBS and disk (waiting for Matti's thesis defence)


    # Deleted to save space, but not invalidated in DBS due to parentage issues (v13 is a child dataset)
    #add("skim", "v13_original", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
            #})

    #add("skim", "v13", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_skim_v13-9cdf0eb8900aa637b61ccc82f152c6ed/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_skim_v13-c61b4e9bdf1ffeec75e33c6424b25cdb/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_skim_v13-a425e273cc15b943339437b345a2e98d/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_skim_v13_2-46d2cb331c30c7a82a30d66aa5cd1785/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_skim_v13_2-5838507f6e8a0c16a243ee2686c03016/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_skim_v13-46d2cb331c30c7a82a30d66aa5cd1785/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_skim_v13-0cb7eab6610b4f67800452f2685ea477/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_skim_v13-23d3791c980dcc1a59a485ca5d8ad22e/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_skim_v13-23d3791c980dcc1a59a485ca5d8ad22e/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_skim_v13-19b79d7dc2d65f948070055429ce37d3/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-4bab7d1d8563adf9003d830a78cf7377/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13_2-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_copy_v13-4711bccee3dc7cbc7d82816b57065063/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13_2-84f25f495c5d27ed442aae803bfc12e6/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13_2-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_skim_v13-6ce8de2c5b6c0c9ed414998577b7e28d/USER",
        #})

    #add("embedding", "v13_1", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_1-947a4a88c33687e763c591af079fc279/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    ## Deleted to save space, but not invalidated in DBS due to parentage issues (v13_3 is a child dataset)
    #add("embedding", "v13_2", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_2-947a4a88c33687e763c591af079fc279/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    ## Deleted to save space, but not invalidated in DBS due to parentage issues (v13_3_seedTest1 is a child dataset)
    #add("embedding", "v13_2_seedTest1", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_2_seedTest1-947a4a88c33687e763c591af079fc279/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    ## Deleted to save space, but not invalidated in DBS due to parentage issues (v13_3_seedTest2 is a child dataset)
    #add("embedding", "v13_2_seedTest2", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_2_seedTest2-947a4a88c33687e763c591af079fc279/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    ## Deleted to save space, but not invalidated in DBS due to parentage issues (v13_3_seedTest3 is a child dataset)
    #add("embedding", "v13_2_seedTest3", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_2_seedTest3-947a4a88c33687e763c591af079fc279/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_2_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})


    #add("embedding", "v13_3", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3-4711bccee3dc7cbc7d82816b57065063/USER",
        #})

    #add("embedding", "v13_3_seedTest1", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest1-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest1-4711bccee3dc7cbc7d82816b57065063/USER",
        #})

    #add("embedding", "v13_3_seedTest2", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest2-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest2-4711bccee3dc7cbc7d82816b57065063/USER",
        #})

    #add("embedding", "v13_3_seedTest3", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest3-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_copy_v13_3_seedTest3-4711bccee3dc7cbc7d82816b57065063/USER",
        #})

    #add("embedding", "v13_3_seedTest4", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_3_seedTest4-947a4a88c33687e763c591af079fc279/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest4-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    #add("embedding", "v13_3_seedTest5", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_3_seedTest5-947a4a88c33687e763c591af079fc279/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest5-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    #add("embedding", "v13_3_seedTest6", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_3_seedTest6-947a4a88c33687e763c591af079fc279/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest6-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    #add("embedding", "v13_3_seedTest7", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_3_seedTest7-947a4a88c33687e763c591af079fc279/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest7-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    #add("embedding", "v13_3_seedTest8", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_3_seedTest8-947a4a88c33687e763c591af079fc279/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest8-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    #add("embedding", "v13_3_seedTest9", {
        #"SingleMu_Mu_160431-163261_May10":     "/SingleMu/local-May10ReReco_v1_AOD_160431_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_163270-163869_May10":     "/SingleMu/local-May10ReReco_v1_AOD_163270_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_165088-166150_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_165088_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166161-166164_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166161_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166346-166346_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166346_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_166374-167043_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_166374_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_167078-167913_Prompt":    "/SingleMu/local-PromptReco_v4_AOD_167078_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_170722-172619_Aug05":     "/SingleMu/local-05Aug2011_v1_AOD_170722_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_172620-173198_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_172620_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"SingleMu_Mu_173236-173692_Prompt":    "/SingleMu/local-PromptReco_v6_AOD_173236_tauembedding_embedding_v13_3_seedTest9-947a4a88c33687e763c591af079fc279/USER",
        #"TTToHplusBWB_M80_Summer11":           "/TTToHplusBWB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M90_Summer11":           "/TTToHplusBWB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M100_Summer11":          "/TTToHplusBWB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M120_Summer11":          "/TTToHplusBWB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M140_Summer11":          "/TTToHplusBWB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M150_Summer11":          "/TTToHplusBWB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M155_Summer11":          "/TTToHplusBWB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBWB_M160_Summer11":          "/TTToHplusBWB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M80_Summer11":      "/TTToHplusBHminusB_M-80_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M90_Summer11":      "/TTToHplusBHminusB_M-90_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M100_Summer11":     "/TTToHplusBHminusB_M-100_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M120_Summer11":     "/TTToHplusBHminusB_M-120_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M140_Summer11":     "/TTToHplusBHminusB_M-140_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M150_Summer11":     "/TTToHplusBHminusB_M-150_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M155_Summer11":     "/TTToHplusBHminusB_M-155_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTToHplusBHminusB_M160_Summer11":     "/TTToHplusBHminusB_M-160_7TeV-pythia6-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"TTJets_TuneZ2_Summer11":              "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"W3Jets_TuneZ2_Summer11":              "/W3Jets_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"DYJetsToLL_M50_TuneZ2_Summer11":      "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_t-channel_TuneZ2_Summer11":         "/T_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_t-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_tW-channel_TuneZ2_Summer11":        "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_tW-channel_TuneZ2_Summer11":     "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"T_s-channel_TuneZ2_Summer11":         "/T_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"Tbar_s-channel_TuneZ2_Summer11":      "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WW_TuneZ2_Summer11":                  "/WW_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"WZ_TuneZ2_Summer11":                  "/WZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"ZZ_TuneZ2_Summer11":                  "/ZZ_TuneZ2_7TeV_pythia6_tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #"QCD_Pt20_MuEnriched_TuneZ2_Summer11": "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest9-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})

    #add("embedding", "v13_3_seedTest10", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest10-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest11", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest11-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest12", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest12-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest13", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest13-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest14", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest14-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest15", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest15-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest16", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest16-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest17", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest17-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest18", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest18-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest19", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest19-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest20", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest20-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest21", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest21-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest22", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest22-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest23", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest23-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest24", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest24-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest25", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest25-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest26", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest26-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest27", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest27-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest28", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest28-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})
    #add("embedding", "v13_3_seedTest29", {
        #"WJets_TuneZ2_Summer11":               "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/local-Summer11_PU_S4_START42_V11_v1_AODSIM_tauembedding_embedding_v13_3_seedTest29-22559ec2c5e66c0c33625ecb67add84e/USER",
        #})


    #add("embedding", "v44_1_isoTest", {
        #"TTJets_TuneZ2_Fall11":               "/TTJets_TuneZ2_7TeV-madgraph-tauola/local-Tauembedding_skim_v44_1_TTJets_TuneZ2_Fall11-6b597f00b6cbae0a4858f059d847b478/USER",
        #})
