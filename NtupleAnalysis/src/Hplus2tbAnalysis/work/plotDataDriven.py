#!/usr/bin/env python
'''
DESCRIPTION:


USAGE:
./plotDataDriven.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
./plotDataDriven.py -m Hplus2tbAnalysis_3bjets40_MVA0p88_MVA0p88_TopMassCutOff600GeV_180113_050540 -n FakeBMeasurement_PreSel_3bjets40_SigSel_MVA0p85_InvSel_EE2CSVM_MVA0p60to085_180120_092605/FakeB/ --gridX --gridY --logY --signal 500
/plotDataDriven.py -m <mcrab1> -n <mcrab2>/FakeBMeasurement/ --gridX --gridY --logY --useMC --unblind


LAST USED:
./plotDataDriven.py -m Hplus2tbAnalysis_NewLeptonVeto_PreSel_3bjets40_SigSel_MVA0p85_180125_123633 -n FakeBMeasurement_NewLeptonVeto_PreSel_3bjets40_SigSel_MVA0p85_InvSel_EE2CSVM_MVA0p60to085_180125_123834/FakeBMeasurement/ --gridX --gridY --logY --unblind

'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
from optparse import OptionParser

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

import HiggsAnalysis.NtupleAnalysis.tools.dataset as dataset
import HiggsAnalysis.NtupleAnalysis.tools.histograms as histograms
import HiggsAnalysis.NtupleAnalysis.tools.counter as counter
import HiggsAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HiggsAnalysis.NtupleAnalysis.tools.styles as styles
import HiggsAnalysis.NtupleAnalysis.tools.plots as plots
import HiggsAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HiggsAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles

#================================================================================================ 
# Function Definition
#================================================================================================ 
def Print(msg, printHeader=False):
    fName = __file__.split("/")[-1]
    if printHeader==True:
        print "=== ", fName
        print "\t", msg
    else:
        print "\t", msg
    return

def Verbose(msg, printHeader=True, verbose=False):
    if not opts.verbose:
        return
    Print(msg, printHeader)
    return

def rchop(myString, endString):
  if myString.endswith(endString):
    return myString[:-len(endString)]
  return myString

def GetLumi(datasetsMgr):
    lumi = 0.0
    for d in datasetsMgr.getAllDatasets():
        if d.isMC():
            continue
        else:
            lumi += d.getLuminosity()
    Verbose("Luminosity = %s (pb)" % (lumi), True)
    return lumi

def GetListOfEwkDatasets():
    return  ["TT", "WJetsToQQ_HT_600ToInf", "SingleTop", "DYJetsToQQHT", "TTZToQQ",  "TTWJetsToQQ", "Diboson", "TTTT"]

def GetDatasetsFromDir(opts, otherDir=False):
    
    myDir    = opts.mcrab1
    analysis = opts.analysisName
    if otherDir:
        myDir    = opts.mcrab2        
    datasets = dataset.getDatasetsFromMulticrabDirs([myDir],
                                                    dataEra=opts.dataEra,
                                                    searchMode=opts.searchMode, 
                                                    analysisName=analysis,
                                                    optimizationMode=opts.optMode)
    return datasets
    
def GetBinWidthMinMax(binList):
    if not isinstance(binList, list):
        raise Exception("Argument is not a list instance!")

    minWidth = +1e6
    maxWidth = -1e6
    # For-loop: All bin values (centre)
    for i in range(0, len(binList)-1):
        j = i + 1
        iBin = binList[i]
        jBin = binList[j]
        wBin = jBin-iBin
        if wBin < minWidth:
            minWidth = wBin

        if wBin > maxWidth:
            maxWidth = wBin
    return minWidth, maxWidth

def main(opts):

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    style.setLogX(opts.logX)
    # If you want BOTH pads (main and ratio) in log scale 
    if 0:
        style.setLogY(opts.logY) 

    # Overwrite default legends
    plots._legendLabels["MCStatError"] = "Bkg. stat."
    plots._legendLabels["MCStatSystError"] = "Bkg. stat.#oplussyst."
    plots._legendLabels["BackgroundStatError"] = "Bkg. stat. unc"
    plots._legendLabels["BackgroundStatSystError"] = "Bkg. stat.#oplussyst. unc."
    
    # Define optimisatio modes to run on
    optModes = [""] #["", "OptChiSqrCutValue50", "OptChiSqrCutValue100"]

    if opts.optMode != None:
        optModes = [opts.optMode]
        
    # For-loop: All opt Mode
    for opt in optModes:
        opts.optMode = opt

        # Setup & configure the dataset manager 
        dsetMgr1 = GetDatasetsFromDir(opts, False) 
        dsetMgr2 = GetDatasetsFromDir(opts, True)

        # Setup the dataset managers
        dsetMgr1.updateNAllEventsToPUWeighted()
        dsetMgr2.updateNAllEventsToPUWeighted()

        # Load luminosities
        dsetMgr1.loadLuminosities() # from lumi.json
        # dsetMgr2.loadLuminosities()

        # Print PSets. Perhaps i can use this to ensure parameters are matching!
        if 0:
            dsetMgr1.printSelections()
            dsetMgr2.printSelections()
            PrintPSet("FakeBMeasurement", dsetMgr1)
            PrintPSet("TopSelectionBDT" , dsetMgr2)

        # Remove datasets with overlap?
        removeList = ["QCD-b"]
        dsetDY     = "DYJetsToQQ_HT180"
        dsetZJ     = "ZJetsToQQ_HT600toInf"
        dsetRM     = dsetZJ # datasets with overlap
        removeList.append(dsetRM)

        # Set/Overwrite cross-sections. Remove all but 1 signal mass 
        for d in dsetMgr1.getAllDatasets():
            if "ChargedHiggs" in d.getName():
                dsetMgr1.getDataset(d.getName()).setCrossSection(1.0) # ATLAS 13 TeV H->tb exclusion limits
                if d.getName() != opts.signal:
                    removeList.append(d.getName())

        # Print useful information?
        if opts.verbose:
            dsetMgr1.PrintCrossSections()
            dsetMgr1.PrintLuminosities()
            dsetMgr2.PrintCrossSections()
            dsetMgr2.PrintLuminosities()

        # Merge histograms
        plots.mergeRenameReorderForDataMC(dsetMgr1) 
        # plots.mergeRenameReorderForDataMC(dsetMgr2) 
   
        # Get the luminosity
        if opts.intLumi < 0:
            opts.intLumi = dsetMgr1.getDataset("Data").getLuminosity()

        # Custom Filtering of datasets 
        for i, d in enumerate(removeList, 1):
            msg = "Removing datasets %s from dataset manager" % (ShellStyles.NoteStyle() + d + ShellStyles.NormalStyle())
            Verbose(msg, i==1)
            dsetMgr1.remove(filter(lambda name: d == name, dsetMgr1.getAllDatasetNames()))

        # Replace MC datasets with data-driven
        if not opts.useMC: #fixme
            replaceQCD(dsetMgr1, dsetMgr2, "FakeBMeasurementTrijetMass", "FakeB") #dsetMgr1 now contains "FakeB" pseudo-dataset

        # Print dataset information
        dsetMgr1.PrintInfo()
        dsetMgr2.PrintInfo()

        # Definitions
        histoNames  = []    
        #allHistos   = dsetMgr2.getAllDatasets()[0].getDirectoryContent(opts.folder)
        #histoList1  = [h for h in allHistos if "MCEWK" not in h]
        #histoList2  = [h for h in histoList1 if "Purity" not in h]
        #histoList3  = [h for h in histoList2 if "Bjet" not in h]
        #histoList4  = [h for h in histoList3 if "Jet" not in h]
        #histoList5  = [h for h in histoList4 if "Njet" not in h]
        #histoList6  = [h for h in histoList5 if "Njet" not in h]
        #histoList7  = [h for h in histoList6 if "MVA" not in h]
        #histoList8  = [h for h in histoList7 if "Sub" not in h]
        #histoList9  = [h for h in histoList8 if "Dijet" not in h]
        #histoPaths  = [opts.folder + "/" + h for h in histoList9]
        
        allHistos   = dsetMgr1.getAllDatasets()[0].getDirectoryContent(opts.folder)
        histoPaths  = []
        for h in allHistos:
            if "StandardSelections" in h or "Subldg" in h or "BJetEta" in h or "BtagDiscriminator" in h or "JetPt" in h or "JetEta" in h:
                continue
            if "Njets" in h or "NBjets" in h or "Dijet" in h or "WMassRatio" in h or h == "MET" or "METPhi" in h or "BJetPt" in h:
                continue
            # Please add me
            if "LdgTrijetBjetPt" in h or "LdgTrijetBjetEta" in h or "TetrajetBjetPt" in h or "TetrajetBjetEta" in h or "NVertices" in h or "HT" in h:
                continue
            histoPaths.append(os.path.join(opts.folder, h))


        # For-loop: All histograms in list
        for i, hName in enumerate(histoPaths, 1):
            
            msg = "{:<9} {:>3} {:<1} {:<3} {:<50}".format("Histogram", "%i" % i, "/", "%s:" % (len(histoPaths)), hName)
            Print(ShellStyles.SuccessStyle() + msg + ShellStyles.NormalStyle(), i==1)

            PlotHistogram(dsetMgr1, hName, opts)

    Print("All plots saved under directory %s" % (opts.saveDir), True)
    return

def GetHistoKwargs(hName, opts):
    '''
    Dictionary with 
    key   = histogramName
    value = kwargs
    '''
    ymaxF  = 1.2
    ymin   = 1e-1
    kwargs = {
        "ratioCreateLegend": True,
        "ratioType"        : "errorScale",
        "ratioErrorOptions": {"numeratorStatSyst": False},
        "ratioMoveLegend"  : {"dx": -0.51, "dy": 0.03, "dh": -0.05},
        "ylabel"           : "Events / %.0f",
        "rebinX"           : 1,
        "rebinY"           : None,
        "ratioYlabel"      : "Data/Bkg. ",
        "ratio"            : True, 
        "stackMCHistograms": True,
        "ratioInvert"      : False,
        "addMCUncertainty" : True, 
        "addLuminosityText": True,
        "addCmsText"       : True,
        "cmsExtraText"     : "Preliminary",
        "opts"             : {"ymin": ymin, "ymaxfactor": ymaxF},
        "opts2"            : {"ymin": 0.2, "ymax": 2.0-0.2},
        "log"              : opts.logY,
        "moveLegend"       : {"dx": -0.1, "dy": -0.01, "dh": 0.1},
        "cutBoxY"          : {"cutValue": 1.2, "fillColor": 16, "box": False, "line": True, "greaterThan": True, "mainCanvas": False, "ratioCanvas": False}
        }

    if "MET" in hName:
        myBins = []
        for j in range(0, 100, 10):
            myBins.append(j)
        for k in range(100, 200, 20):
            myBins.append(k)
        for k in range(200, 300, 50):
            myBins.append(k)
        for k in range(300, 400+100, 100):
            myBins.append(k)
        units            = "GeV/c"
        binWmin, binWmax = GetBinWidthMinMax(myBins)
        kwargs["ylabel"] = "Events / %.0f-%.0f %s" % (binWmin, binWmax, units)
        kwargs["xlabel"] = "E_{T}^{miss} (%s)" % units
        kwargs["rebinX"] = myBins
        kwargs["opts"]   = {"ymin": ymin, "ymaxfactor": ymaxF}
        kwargs["log"]    = True
    if "NVertices" in hName:
        kwargs["ylabel"] = "Events / %.0f"
        kwargs["xlabel"] = "Vertices"
    if "Njets" in hName:                
        kwargs["ylabel"] = "Events / %.0f"
        kwargs["xlabel"] = "Jets Multiplicity"
        kwargs["cutBox"] = {"cutValue": 7.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 7.0, "xmax": +16.0, "ymin": ymin, "ymaxfactor": ymaxF}
        ROOT.gStyle.SetNdivisions(10, "X")
    if "BtagDiscriminator" in hName:                
        kwargs["ylabel"]     = "Events / %.2f"
        kwargs["xlabel"]     = "b-tag discriminator"
        kwargs["cutBox"]     = {"cutValue": 0.8484, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["moveLegend"] = {"dx": -0.5, "dy": -0.01, "dh": 0.0}
        kwargs["opts"]       = {"xmin": 0.0, "xmax": 1.05, "ymin": ymin, "ymaxfactor": ymaxF}
    if "HT" in hName:
        units            = "GeV"
        kwargs["ylabel"] = "Events / %.0f " + units
        kwargs["xlabel"] = "H_{T} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 500.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 5
        kwargs["opts"]   = {"xmin": 500.0, "xmax": 3000, "ymin": ymin, "ymaxfactor": ymaxF}
    if "LdgTrijetPt" in hName:
        units = "GeV/c"
        kwargs["ylabel"] = "Events / %.0f " + units
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmax": +800.0, "ymin": ymin, "ymaxfactor": ymaxF}
        kwargs["rebinX"] = 2
        kwargs["log"]    = True
    if "LdgTrijetMass" in hName:
        startBlind       = 140
        endBlind         = 200
        kwargs["rebinX"] = 1
        units            = "GeV/c^{2}"
        kwargs["ylabel"] = "Events / %.0f " + units
        kwargs["xlabel"] = "m_{jjb} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 173.21, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 50.0, "xmax": 350.0, "ymin": ymin, "ymaxfactor": ymaxF}
        kwargs["log"]    = True
        # Blind data
        kwargs["blindingRangeString"] = "%s-%s" % (startBlind, endBlind)
        kwargs["moveBlindedText"]     = {"dx": -0.22, "dy": +0.08, "dh": -0.1}
    if "LdgTrijetBjetPt" in hName:
        units            = "GeV/c"
        kwargs["ylabel"] = "Events / %.0f " + units
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["log"]    = True
    if "LdgTrijetBjetEta" in hName:
        kwargs["ylabel"] = "Events / %.2f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -2.5, "xmax": +2.5, "ymin": ymin, "ymaxfactor": ymaxF}
    if "LdgTetrajetPt" in hName:
        units            = "GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)"  % units
        kwargs["opts"]   = {"xmax": +900.0, "ymin": ymin, "ymaxfactor": ymaxF}
        kwargs["log"]    = True
        myBins = []
        for j in range(0, 500, 20):
            myBins.append(j)
        for k in range(500, 950, 50):
            myBins.append(k)
        kwargs["rebinX"] = myBins #2
        binWmin, binWmax = GetBinWidthMinMax(myBins)
        kwargs["ylabel"] = "Events / %.0f-%.0f %s" % (binWmin, binWmax, units)
        # kwargs["ylabel"] = "Events / %.0f " + units
    if "LdgTetrajetMass" in hName:
        myBins = []
        for j in range(0, 1000, 50):
            myBins.append(j)
        for k in range(1000, 2000, 100):
            myBins.append(k)
        for l in range(2000, 4000+2000, 200):
            myBins.append(l)

        ROOT.gStyle.SetNdivisions(5, "X")
        startBlind       = 150  # 135 v. sensitive to bin-width!
        endBlind         = 2500 #v. sensitive to bin-width!
        kwargs["rebinX"] = myBins
        units            = "GeV/c^{2}"
        binWmin, binWmax = GetBinWidthMinMax(myBins)
        kwargs["ylabel"] = "Events / %.0f-%.0f %s" % (binWmin, binWmax, units)
        kwargs["log"]    = True
        kwargs["xlabel"] = "m_{jjbb} (%s)"  % units
        kwargs["cutBox"] = {"cutValue": 500.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": endBlind, "ymin": ymin, "ymaxfactor": ymaxF}
        if "AllSelections" in hName:
            kwargs["blindingRangeString"] = "%s-%s" % (startBlind, endBlind) #ale
            kwargs["moveBlindedText"]     = {"dx": -0.22, "dy": +0.08, "dh": -0.12}
    if "TetrajetBjetPt" in hName:
        units            = "GeV/c"
        kwargs["ylabel"] = "Events / %.0f " + units
        kwargs["xlabel"] = "p_{T} (%s)"  % units

    if "TetrajetBjetEta" in hName:
        kwargs["ylabel"] = "Events / %.2f"
        kwargs["xlabel"] = "#eta"
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -2.5, "xmax": +2.5, "ymin": ymin, "ymaxfactor": ymaxF}

    if opts.unblind:
        key = "blindingRangeString"
        if key in kwargs.keys():
            del kwargs[key]

    if kwargs["log"]==True or opts.logY == True:
        kwargs["opts"]["ymaxfactor"] = 8.0

    return kwargs
    
def ApplyBlinding(myObject, blindedRange = []):
    '''
    myObject must be an instance of:
    h=histograms.Histo(rootHisto, "Label")
    and the rooHistos is an instance of:
    rootHisto = p.histoMgr.getHisto("HistoName").getRootHisto()
    '''

    if len(blindedRange) != 2:
        msg = "Blinded range list requires exactly 2 values (got %s)" % len(blindedRange) 
        raise Exception(ShellStyles.ErrorStyle() + msg + ShellStyles.NormalStyle())

    # Definitions
    myMin = None
    myMax = None
    myHisto = myObject.getRootHisto()

    # For-loop: All histogram bins
    for i in range (1, myHisto.GetNbinsX()+1):
        myUpEdge  = myHisto.GetXaxis().GetBinUpEdge(i)
        myLowEdge = myHisto.GetXaxis().GetBinLowEdge(i)
        
        # Define conditions
        c1 = (myLowEdge >= blindedRange[0] and myLowEdge <= blindedRange[1]) 
        c2 = (myUpEdge  >= blindedRange[0] and myUpEdge  <= blindedRange[1])
        c3 = (myLowEdge <= blindedRange[0] and myUpEdge  >= blindedRange[1])

        # Blind if any edge of the current bin is inside the blinded range or if bin spans over the blinded range
        if ( c1 or c2 or c3):
            if myMin == None or myLowEdge < myMin:
                myMin = myLowEdge
            if myMax == None or myUpEdge > myMax:
                myMax = myUpEdge
            #  Blind data by setting bin content to -1.0
            myHisto.SetBinContent(i, -1.0)
            myHisto.SetBinError(i, 0.0)

    if myMin == None:
        return None
    
    # Prepare blinding string for printing on canvas
    myMinFormat = "%" + "d"
    myMaxFormat = "%" + "d"
    if abs(myMin) < 1.0 and abs(myMin) > 0.00000001:
        myMinFormat = "%%.%df" % (abs(int(log10(myMin)))+1)
    if abs(myMax) < 1.0  and abs(myMax) > 0.00000001:
        myMaxFormat = "%%.%df" % (abs(int(log10(myMax)))+1)
    bString = myMinFormat%myMin+"-"+myMaxFormat%myMax
    return bString


def PlotHistogram(dsetMgr, histoName, opts):

    # Get kistogram argumetns
    kwargs = GetHistoKwargs(histoName, opts)
    saveName = histoName.replace(opts.folder + "/", "")

    # Create the plotting object (Data, "FakeB")
    p1 = plots.DataMCPlot(dsetMgr, histoName, saveFormats=[])

    # Copy dataset manager before changing datasets. Keep only EWK (GenuineB) datasets
    datasetMgr = dsetMgr.deepCopy()
    datasetMgr.selectAndReorder(GetListOfEwkDatasets())
        
    # Create the MCPlot for the EWKGenuineB histograms
    histoNameGenuineB = histoName.replace(opts.folder, opts.folder + "EWKGenuineB")
    p2 = plots.MCPlot(datasetMgr, histoNameGenuineB, normalizeToLumi=opts.intLumi, saveFormats=[])
    # plots.drawPlot(p2, saveName + "_tmp", **{"rebiX": kwargs["rebinX"]})
    # p2.histoMgr.forEachHisto(lambda h: h.getRootHisto().RebinX(kwargs["rebinX"]))

    # Add the datasets to be included in the plot
    myStackList = []

    # Data-driven FakeB background
    if not opts.useMC:
        hFakeB  = p1.histoMgr.getHisto("FakeB").getRootHisto()
        hhFakeB = histograms.Histo(hFakeB, "FakeB", legendLabel="Fake-b")
        hhFakeB.setIsDataMC(isData=False, isMC=True)
        myStackList.append(hhFakeB)
    else:
        hQCD  = p1.histoMgr.getHisto("QCD").getRootHisto()
        hhQCD = histograms.Histo(hQCD, "QCD", legendLabel="QCD")
        hhQCD.setIsDataMC(isData=False, isMC=True)
        myStackList.append(hhQCD)

    # EWK GenuineB background (Replace all EWK histos with GenuineB histos)
    ewkNameList  = GetListOfEwkDatasets()
    ewkHistoList = []
    # For-loop: All EWK datasets 
    for dataset in ewkNameList:
        h = p2.histoMgr.getHisto(dataset).getRootHisto()
        hh = histograms.Histo(h, dataset,  plots._legendLabels[dataset])
        hh.setIsDataMC(isData=False, isMC=True)
        myStackList.append(hh)

    # Collision data
    hData  = p1.histoMgr.getHisto("Data").getRootHisto()
    hhData = histograms.Histo(hData, "Data")
    hhData.setIsDataMC(isData=True, isMC=False)
    myStackList.insert(0, hhData)

    # Signal
    hSignal  = p1.histoMgr.getHisto(opts.signal).getRootHisto()
    hhSignal = histograms.Histo(hSignal, opts.signal, plots._legendLabels[opts.signal])
    hhSignal.setIsDataMC(isData=False, isMC=True)
    myStackList.insert(1, hhSignal)

    # Create the final plot by passing the histogram list
    p3 = plots.DataMCPlot2(myStackList, saveFormats=[])
    p3.setLuminosity(opts.intLumi)
    p3.setDefaultStyles()

    # Apply blinding of data in Signal Region (After creating the plot)
    if "blindingRangeString" in kwargs:
        startBlind = float(kwargs["blindingRangeString"].split("-")[1])
        endBlind   = float(kwargs["blindingRangeString"].split("-")[0])
        plots.partiallyBlind(p3, maxShownValue=startBlind, minShownValue=endBlind, invert=True, moveBlindedText=kwargs["moveBlindedText"])

    # Draw and save the plot
    plots.drawPlot(p3, saveName, **kwargs)
    SavePlot(p3, saveName, os.path.join(opts.saveDir, opts.optMode), saveFormats = [".png"])
    return

def PrintPSet(selection, dsetMgr):
    selection = "\"%s\":"  % (selection)
    thePSets = dsetMgr.getAllDatasets()[0].getParameterSet()

    # First drop everything before the selection
    thePSet_1 = thePSets.split(selection)[-1]

    # Then drop everything after the selection
    thePSet_2 = thePSet_1.split("},")[0]

    # Final touch
    thePSet = selection + thePSet_2

    Print(thePSet, True)
    return

def getHisto(dsetMgr, datasetName, histoName):
    Verbose("getHisto()", True)

    h1 = dsetMgr.getDataset(datasetName).getDatasetRootHisto(histoName)
    h1.setName(datasetName)
    return h1

def replaceQCD(dMgr1, dMgr2, newName, newLabel="FakeB"):
    '''
    Replaces the QCD dataset with 
    a data-driven pseudo-dataset
    '''
    
    # Define variables
    oldName    = "QCD"
    newDataset = dMgr2.getDataset(newName)
    dMgr1.append(newDataset)
    newDataset.setName(newLabel)

    # Get list of dataset names
    names = dMgr1.getAllDatasetNames()

    # Return the index in the list of the first dataset whose name is oldName.
    index = names.index(oldName)

    # Remove the dataset at the given position in the list, and return it. 
    names.pop(index)
    
    # Insert the dataset to the given position  (index) of the list
    names.insert(index, newDataset.getName())

    # Remove from the dataset manager the oldName
    dMgr1.remove(oldName)
    names.pop(names.index(newDataset.getName()))
    
    # Select and reorder Datasets. 
    # This method can be used to either select a set of dataset.Dataset objects. reorder them, or both.
    dMgr1.selectAndReorder(names)
    return

def SavePlot(plot, plotName, saveDir, saveFormats = [".png", ".pdf"]):
    # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats, 0):
        saveNameURL  = saveName + ext
        saveNameURL  = saveNameURL.replace("/publicweb/a/aattikis/", "http://home.fnal.gov/~aattikis/")
        if opts.url:
            Verbose(saveNameURL, False) #i==0)
            opts.saveDir = os.path.dirname(saveNameURL) + "/"
        else:
            Verbose(saveName + ext, False) #i==0)
            opts.saveDir = os.path.dirname(saveName) + "/"
        plot.saveAs(saveName, formats=saveFormats)
    return


#================================================================================================ 
# Main
#================================================================================================ 
if __name__ == "__main__":
    '''
    https://docs.python.org/3/library/argparse.html
 
    name or flags...: Either a name or a list of option strings, e.g. foo or -f, --foo.
    action..........: The basic type of action to be taken when this argument is encountered at the command line.
    nargs...........: The number of command-line arguments that should be consumed.
    const...........: A constant value required by some action and nargs selections.
    default.........: The value produced if the argument is absent from the command line.
    type............: The type to which the command-line argument should be converted.
    choices.........: A container of the allowable values for the argument.
    required........: Whether or not the command-line option may be omitted (optionals only).
    help............: A brief description of what the argument does.
    metavar.........: A name for the argument in usage messages.
    dest............: The name of the attribute to be added to the object returned by parse_args().
    '''
    
    # Default Settings
    ANALYSISNAME = "Hplus2tbAnalysis"
    USEMC        = False
    SEARCHMODE   = "80to1000"
    DATAERA      = "Run2016"
    OPTMODE      = None
    BATCHMODE    = True
    INTLUMI      = -1.0
    SIGNALMASS   = 800
    SIGNAL       = None
    URL          = False
    SAVEDIR      = "/publicweb/a/aattikis/"
    VERBOSE      = False
    GRIDX        = False
    UNBLIND      = False
    GRIDY        = False
    LOGX         = False
    LOGY         = False
    FOLDER       = "ForDataDrivenCtrlPlots"

    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab1", dest="mcrab1", action="store", 
                      help="Path to the multicrab directory with the Hplus2tbAnalysis")

    parser.add_option("-n", "--mcrab2", dest="mcrab2", action="store", 
                      help="Path to the multicrab directory with the data-driven pseudo-datasets (e.g. FakeB)")

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--analysisName", dest="analysisName", type="string", default=ANALYSISNAME,
                      help="Override default analysisName [default: %s]" % ANALYSISNAME)

    parser.add_option("--useMC", dest="useMC", action="store_true", default=USEMC,
                      help="Use all backgrounds from MC. Do not use data-driven background estimation [default: %s]" % USEMC)

    parser.add_option("--unblind", dest="unblind", action="store_true", default=UNBLIND,
                      help="Switch off blinging of data in the signal region [default: %s]" % UNBLIND)

    parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
                      help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--signalMass", dest="signalMass", type=int, default=SIGNALMASS, 
                     help="Mass value of signal to use [default: %s]" % SIGNALMASS)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX,
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY,
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)

    parser.add_option("--logX", dest="logX", action="store_true", default=LOGX,
                      help="Set x-axis to logarithm scale [default: %s]" % LOGX)

    parser.add_option("--logY", dest="logY", action="store_true", default=LOGY,
                      help="Set y-axis to logarithm scale [default: %s]" % LOGY)

    parser.add_option("--folder", dest="folder", action="store_true", default=FOLDER,
                      help="Folder inside the ROOT files where all histograms for plotting are located [default: %s]" % FOLDER)


    (opts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    if opts.mcrab1 == None:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        parser.print_help()
        #print __doc__
        sys.exit(1)
    elif opts.mcrab2 == None:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        parser.print_help()
        #print __doc__
        sys.exit(1)
    else:
        mcrabDir = rchop(opts.mcrab1, "/")
        if len(mcrabDir.split("/")) > 1:
            mcrabDir = mcrabDir.split("/")[-1]
        opts.saveDir += mcrabDir + "/DataDriven/"


    # Sanity check
    allowedMass = [180, 200, 220, 250, 300, 350, 400, 500, 800, 1000, 1500, 2000, 2500, 3000, 5000, 7000]
    if opts.signalMass!=0 and opts.signalMass not in allowedMass:
        Print("Invalid signal mass point (=%.0f) selected! Please select one of the following:" % (opts.signalMass), True)
        for m in allowedMass:
            Print(m, False)
        sys.exit()
    else:
        opts.signal = "ChargedHiggs_HplusTB_HplusToTB_M_%i" % opts.signalMass

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotDataDriven.py: Press any key to quit ROOT ...")
