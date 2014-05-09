#!/usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import sys
import os
import time
from optparse import OptionParser
import cProfile

import HiggsAnalysis.HeavyChHiggsToTauNu.tools.dataset as dataset
import HiggsAnalysis.HeavyChHiggsToTauNu.tools.plots as plots
from HiggsAnalysis.HeavyChHiggsToTauNu.tools.analysisModuleSelector import *
from HiggsAnalysis.HeavyChHiggsToTauNu.qcdInverted.qcdInvertedResult import *
from HiggsAnalysis.HeavyChHiggsToTauNu.tools.ShellStyles import *
from HiggsAnalysis.HeavyChHiggsToTauNu.tools.pseudoMultiCrabCreator import *

myNormalizationFactorSource = "QCDInvertedNormalizationFactors.py"

def doNominalModule(myMulticrabDir,era,searchMode,optimizationMode,myOutputCreator,myShapeString,myNormFactors,myDisplayStatus):
    # Construct info string of module
    myModuleInfoString = "%s_%s_%s"%(era, searchMode, optimizationMode)
    # Obtain dataset manager
    dsetMgrCreator = dataset.readFromMulticrabCfg(directory=myMulticrabDir)
    dsetMgr = dsetMgrCreator.createDatasetManager(dataEra=era,searchMode=searchMode,optimizationMode=optimizationMode)
    # Do the usual normalisation
    dsetMgr.updateNAllEventsToPUWeighted()
    dsetMgr.loadLuminosities()
    plots.mergeRenameReorderForDataMC(dsetMgr)
    dsetMgr.merge("EWK", ["TTJets","WJets","DYJetsToLL","SingleTop","Diboson"])
    # Obtain luminosity
    myLuminosity = dsetMgr.getDataset("Data").getLuminosity()
    # Print info so that user can check that merge went correct
    if myDisplayStatus:
        dsetMgr.printDatasetTree()
        print "Luminosity = %f 1/fb"%(myLuminosity / 1000.0)
        print
        myDisplayStatus = False
    # Gather results
    # Create containers for results
    myModuleResults = PseudoMultiCrabModule(dsetMgr, era, searchMode, optimizationMode)
    myQCDNormalizationSystUpResults = PseudoMultiCrabModule(dsetMgr, era, searchMode, optimizationMode, "SystVarQCDNormPlus")
    myQCDNormalizationSystDownResults = PseudoMultiCrabModule(dsetMgr, era, searchMode, optimizationMode, "SystVarQCDNormMinus")
    # Obtain results
    myResult = QCDInvertedResultManager(myShapeString, "AfterCollinearCuts", dsetMgr, myLuminosity, myModuleInfoString, myNormFactors, shapeOnly=False, displayPurityBreakdown=True)
    # Store results
    myModuleResults.addShape(myResult.getShape(), myShapeString)
    myModuleResults.addShape(myResult.getShapeMCEWK(), myShapeString+"_MCEWK")
    myModuleResults.addShape(myResult.getShapePurity(), myShapeString+"_Purity")
    myModuleResults.addDataDrivenControlPlots(myResult.getControlPlots(),myResult.getControlPlotLabels())
    myOutputCreator.addModule(myModuleResults)
    # Up variation of QCD normalization (i.e. ctrl->signal region transition)
    myQCDNormalizationSystUpResults.addShape(myResult.getRegionSystUp(), myShapeString)
    myQCDNormalizationSystUpResults.addDataDrivenControlPlots(myResult.getRegionSystUpCtrlPlots(),myResult.getControlPlotLabelsForQCDSyst())
    myOutputCreator.addModule(myQCDNormalizationSystUpResults)
    # Down variation of QCD normalization (i.e. ctrl->signal region transition)
    myQCDNormalizationSystDownResults.addShape(myResult.getRegionSystDown(), myShapeString)
    myQCDNormalizationSystDownResults.addDataDrivenControlPlots(myResult.getRegionSystDownCtrlPlots(),myResult.getControlPlotLabelsForQCDSyst())
    myOutputCreator.addModule(myQCDNormalizationSystDownResults)
    myResult.delete()
    dsetMgr.close()
    dsetMgrCreator.close()
    ROOT.gROOT.CloseFiles()
    ROOT.gROOT.GetListOfCanvases().Delete()

def doSystematicsVariation(myMulticrabDir,era,searchMode,optimizationMode,syst,myOutputCreator,myShapeString,myNormFactors):
    myModuleInfoString = "%s_%s_%s_%s"%(era, searchMode, optimizationMode,syst)
    dsetMgrCreator = dataset.readFromMulticrabCfg(directory=myMulticrabDir)
    systDsetMgr = dsetMgrCreator.createDatasetManager(dataEra=era,searchMode=searchMode,optimizationMode=optimizationMode,systematicVariation=syst)
    # Do the usual normalisation
    systDsetMgr.updateNAllEventsToPUWeighted()
    systDsetMgr.loadLuminosities()
    plots.mergeRenameReorderForDataMC(systDsetMgr)
    systDsetMgr.merge("EWK", ["TTJets","WJets","DYJetsToLL","SingleTop","Diboson"])
    myLuminosity = systDsetMgr.getDataset("Data").getLuminosity()
    # Obtain results
    mySystModuleResults = PseudoMultiCrabModule(systDsetMgr, era, searchMode, optimizationMode, syst)
    mySystResult = QCDInvertedResultManager(myShapeString, "AfterCollinearCuts", systDsetMgr, myLuminosity, myModuleInfoString, myNormFactors, shapeOnly=False, displayPurityBreakdown=False)
    mySystModuleResults.addShape(mySystResult.getShape(), myShapeString)
    mySystModuleResults.addShape(mySystResult.getShapeMCEWK(), myShapeString+"_MCEWK")
    mySystModuleResults.addShape(mySystResult.getShapePurity(), myShapeString+"_Purity")
    mySystModuleResults.addDataDrivenControlPlots(mySystResult.getControlPlots(),mySystResult.getControlPlotLabels())
    mySystResult.delete()
    ## Save module info
    myOutputCreator.addModule(mySystModuleResults)
    systDsetMgr.close()
    dsetMgrCreator.close()
    ROOT.gROOT.CloseFiles()
    ROOT.gROOT.GetListOfCanvases().Delete()

def printTimeEstimate(globalStart, localStart, nCurrent, nAll):
    myLocalDelta = time.time() - localStart
    myGlobalDelta = time.time() - globalStart
    myEstimate = myGlobalDelta / float(nCurrent) * float(nAll-nCurrent)
    s = "%02d:"%(myEstimate/60)
    myEstimate -= int(myEstimate/60)*60
    s += "%02d"%(myEstimate)
    print "Module finished in %.1f s, estimated time to complete: %s"%(myLocalDelta, s)

if __name__ == "__main__":
    # Obtain normalization factors
    myNormFactors = None
    myNormFactorsSafetyCheck = None
    if os.path.exists(myNormalizationFactorSource):
        myNormFactorsImport = getattr(__import__(myNormalizationFactorSource.replace(".py","")), "QCDInvertedNormalization")
        myNormFactorsSafetyCheck = getattr(__import__(myNormalizationFactorSource.replace(".py","")), "QCDInvertedNormalizationSafetyCheck")
        #QCDInvertedNormalizationSafetyCheck(era)
        myNormFactors = myNormFactorsImport.copy()
    else:
        raise Exception(ErrorLabel()+"Normalisation factors ('%s.py') not found!\nRun script InvertedTauID_Normalization.py to generate the normalization factors."%myNormalizationFactorSource)
    # Object for selecting data eras, search modes, and optimization modes
    myModuleSelector = AnalysisModuleSelector() 
    # Parse command line options
    parser = OptionParser(usage="Usage: %prog [options]",add_help_option=True,conflict_handler="resolve")
    myModuleSelector.addParserOptions(parser)
    parser.add_option("--mdir", dest="multicrabDir", action="store", help="Multicrab directory")
    parser.add_option("--mtonly", dest="transverseMassOnly", action="store_true", default=False, help="Create only transverse mass plots")
    parser.add_option("--invmassonly", dest="invariantMassOnly", action="store_true", default=False, help="Create only invariant mass plots")
    parser.add_option("--test", dest="test", action="store_true", default=False, help="Make short test by limiting number of syst. variations")
    # Add here parser options, if necessary, following line is an example
    #parser.add_option("--showcard", dest="showDatacard", action="store_true", default=False, help="Print datacards also to screen")

    # Parse options
    (opts, args) = parser.parse_args()

    # Obtain multicrab directory
    myMulticrabDir = "."
    if opts.multicrabDir != None:
        myMulticrabDir = opts.multicrabDir
    if not os.path.exists("%s/multicrab.cfg"%myMulticrabDir):
        raise Exception(ErrorLabel()+"No multicrab directory found at path '%s'! Please check path or specify it with --mdir!"%(myMulticrabDir)+NormalStyle())

    myShapes = ["mt","invmass"]
    if opts.transverseMassOnly != None and opts.transverseMassOnly:
        myShapes = ["mt"]
    elif opts.invariantMassOnly != None and opts.invariantMassOnly:
        myShapes = ["invmass"]

    # Obtain dsetMgrCreator and register it to module selector
    dsetMgrCreator = dataset.readFromMulticrabCfg(directory=myMulticrabDir)
    # Obtain systematics names
    mySystematicsNamesRaw = dsetMgrCreator.getSystematicVariationSources()
    mySystematicsNames = []
    for item in mySystematicsNamesRaw:
        mySystematicsNames.append("%sPlus"%item)
        mySystematicsNames.append("%sMinus"%item)
    if opts.test:
        mySystematicsNames = [mySystematicsNames[0]]

    myModuleSelector.setPrimarySource("analysis", dsetMgrCreator)
    # Select modules
    myModuleSelector.doSelect(opts)
    #dsetMgrCreator.close()
    # Loop over era/searchMode/optimizationMode combos
    myDisplayStatus = True
    myTotalModules = myModuleSelector.getSelectedCombinationCount() * (len(mySystematicsNames)+1) * len(myShapes)
    myModuleSelector.printSelectedCombinationCount()
    # Loop over era/searchMode/optimizationMode options

    # Create pseudo-multicrab creator
    myOutputCreator = PseudoMultiCrabCreator("QCDinverted", myMulticrabDir)
    n = 0
    myGlobalStartTime = time.time()
    for massType in myShapes:
        myShapeString = ""
        if massType == "mt":
            myShapeString = "shapeTransverseMass"
        elif massType == "invmass":
            myShapeString = "shapeInvariantMass"
        myOutputCreator.initialize(massType)
        print HighlightStyle()+"Creating dataset for shape: %s%s"%(massType,NormalStyle())
        for era in myModuleSelector.getSelectedEras():
            # Check if normalization coefficients are suitable for era
            myNormFactorsSafetyCheck(era)
            for searchMode in myModuleSelector.getSelectedSearchModes():
                for optimizationMode in myModuleSelector.getSelectedOptimizationModes():
                    myModuleInfoString = "%s_%s_%s"%(era, searchMode, optimizationMode)
                    n += 1
                    print CaptionStyle()+"Module %d/%d: %s/%s%s"%(n,myTotalModules,myModuleInfoString,massType,NormalStyle())
                    myStartTime = time.time()
                    doNominalModule(myMulticrabDir,era,searchMode,optimizationMode,myOutputCreator,myShapeString,myNormFactors,myDisplayStatus)
                    printTimeEstimate(myGlobalStartTime, myStartTime, n, myTotalModules)
                    # Now do the rest of systematics variations
                    for syst in mySystematicsNames:
                        n += 1
                        print CaptionStyle()+"Analyzing systematics variations %d/%d: %s/%s/%s%s"%(n,myTotalModules,myModuleInfoString,syst,massType,NormalStyle())
                        myStartTime = time.time()
                        doSystematicsVariation(myMulticrabDir,era,searchMode,optimizationMode,syst,myOutputCreator,myShapeString,myNormFactors)
                        printTimeEstimate(myGlobalStartTime, myStartTime, n, myTotalModules)
        print "\nPseudo-multicrab ready for mass %s...\n"%massType
    # Create rest of pseudo multicrab directory
    myOutputCreator.finalize()
    print "Average processing time of one module: %.1f s, total elapsed time: %.1f s"%((time.time()-myGlobalStartTime)/float(myTotalModules), (time.time()-myGlobalStartTime))
