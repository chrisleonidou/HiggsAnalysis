#!/usr/bin/env python

from HiggsAnalysis.NtupleAnalysis.main import Process, PSet, Analyzer

process = Process()

# Example of adding a dataset which has its files defined in data/<dataset_name>.txt file
#process.addDatasets(["TTbar_HBWB_HToTauNu_M_160_13TeV_pythia6"])

# Example of adding datasets from a multicrab directory
import sys
if len(sys.argv) != 2:
    print "Usage: ./exampleAnalysis.py <path-to-multicrab-directory>"
    sys.exit(0)
process.addDatasetsFromMulticrab(sys.argv[1])

# Add config
from HiggsAnalysis.NtupleAnalysis.parameters.signalAnalysisParameters import allSelections
process.addAnalyzer("SignalAnalysis", Analyzer("SignalAnalysis", config=allSelections))

# Example of adding an analyzer whose configuration depends on dataVersion
#def createAnalyzer(dataVersion):
    #a = Analyzer("ExampleAnalysis")
    #if dataVersion.isMC():
        #a.tauPtCut = 10
    #else:
        #a.tauPtCut = 20
    #return a
#process.addAnalyzer("test2", createAnalyzer)

# Run the analysis
process.run()

# Run the analysis with PROOF
# By default it uses all cores, but you can give proofWorkers=<N> as a parameter
#process.run(proof=True)
