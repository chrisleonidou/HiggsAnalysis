#!/usr/bin/env python
'''
DESCRIPTION:
Generate all TH1 generated from the Kinematics analyzer (GEN-level info).


USAGE:
./plot_1d.py -m <pseudo_mcrab> [opts]


EXAMPLES:
./plot_1d.py -m <peudo_mcrab> -o "" --url --normaliseToOne


LAST USED:
./plot_1d.py --url --normaliseToOne -m  KinematicsHToHW_17h50m_16Apr2019

'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
import re
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
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux
import HiggsAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HiggsAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles

#================================================================================================ 
# Shell Definitions
#================================================================================================ 
ss = ShellStyles.SuccessStyle()
ns = ShellStyles.NormalStyle()
ts = ShellStyles.NoteStyle()
ls = ShellStyles.HighlightStyle()
hs = ShellStyles.HighlightAltStyle()
es = ShellStyles.ErrorStyle()
styleList = [styles.Style(24, ROOT.kBlack)] + styles.getStyles()

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
    aux.Print(msg, printHeader)
    return

def natural_sort_key(s):
    _nsre = re.compile('([0-9]+)')
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]   

def GetLumi(datasetsMgr):
    Verbose("Determininig Integrated Luminosity")
    
    lumi = 0.0
    for d in datasetsMgr.getAllDatasets():
        if d.isMC():
            continue
        else:
            lumi += d.getLuminosity()
    Verbose("Luminosity = %s (pb)" % (lumi), True)
    return lumi

def GetDatasetsFromDir(opts):
    Verbose("Getting datasets")
    
    if (not opts.includeOnlyTasks and not opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode, 
                                                        analysisName=opts.analysisName,
                                                        optimizationMode=opts.optMode)
    elif (opts.includeOnlyTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        includeOnlyTasks=opts.includeOnlyTasks,
                                                        optimizationMode=opts.optMode)
    elif (opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        excludeTasks=opts.excludeTasks,
                                                        optimizationMode=opts.optMode)
    else:
        raise Exception("This should never be reached")
    return datasets
    

def GetAllDatasetsFromDir(opts):
    datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                    dataEra=opts.dataEra,
                                                    searchMode=opts.searchMode,
                                                    analysisName=opts.analysisName,
                                                    optimizationMode=opts.optMode)
    return datasets

def main(opts):

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    
    optModes = [""]

    if opts.optMode != None:
        optModes = [opts.optMode]

    # For-loop: All optimisation modes
    for opt in optModes:
        opts.optMode = opt

        # Quick and dirty way to get total int lumi
        allDatasetsMgr = GetAllDatasetsFromDir(opts)
        allDatasetsMgr.updateNAllEventsToPUWeighted()
        # allDatasetsMgr.loadLuminosities() # from lumi.json
        plots.mergeRenameReorderForDataMC(allDatasetsMgr) 
        opts.intLumi = GetLumi(allDatasetsMgr)

        # Setup & configure the dataset manager 
        datasetsMgr = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted()
        # datasetsMgr.loadLuminosities() # from lumi.json
        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()

        # Set/Overwrite cross-sections
        print "dsets = " , datasetsMgr.getAllDatasetNames()
        for d in datasetsMgr.getAllDatasets():
                datasetsMgr.getDataset(d.getName()).setCrossSection(1.0)

        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr) 
        
        # Print dataset information
        datasetsMgr.PrintInfo()

        # Re-order datasets in ascending mass 
        newOrder = []
        for i, d in enumerate(datasetsMgr.getAllDatasetNames()):
            newOrder.append(d)
        newOrder.sort(key=natural_sort_key)
        datasetsMgr.selectAndReorder(newOrder)

        # Do the topSelection histos
        folder     = opts.folder
        histoPaths = []
        histoList  = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(folder)
        hList      = [h for h in histoList if "_Vs_" not in h]
        histoPaths = [os.path.join(folder, h) for h in hList]
        for h in histoPaths:
            PlotMC(datasetsMgr, h, opts.intLumi)
    return


def GetHistoKwargs(histo, opts):
    '''
    Dictionary with 
    key   = histogramName
    value = kwargs
    '''
    
    if opts.normaliseToOne:
        yLabel = "Arbitrary Units"
    else:
        yLabel = "Events"
    logY       = False
    yMaxFactor = 3 #1.2 chris

    # Create with default values
    kwargs = {
        "xlabel"           : "x-label",
        "ylabel"           : yLabel,
        "rebinX"           : 1,
        "rebinY"           : 1,
        "ratioYlabel"      : "Data/MC",
        "ratio"            : True, 
        "stackMCHistograms": True,
        "ratioInvert"      : False, 
        "addMCUncertainty" : False, 
        "addLuminosityText": True,
        "addCmsText"       : True,
        "cmsExtraText"     : "Preliminary",
        "opts"             : {"ymin": 0.0, "ymaxfactor": yMaxFactor},
        "opts2"            : {"ymin": 0.0, "ymax": 3.0}, #2
        "log"              : logY,
        "moveLegend"       : {"dx": -0.35, "dy": 0.0, "dh": -0.1}, #chris dx{-0,1 -> -0.7} , dh {-0,1->0.5}
        "cutBox"           : {"cutValue": 0.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        }

    ROOT.gStyle.SetNdivisions(10, "X")
    if "genHT" in histo:
        units            = "GeV"
        format           = "%0.0f " + units
        kwargs["xlabel"] = "H_{T} (%s)" % units
        kwargs["ylabel"] = yLabel + "/ %s " % format
        kwargs["opts"]   = {"xmin": 500.0, "xmax": 2000.0}
        kwargs["log"]    = True
        ROOT.gStyle.SetNdivisions(10, "X")
    elif histo.lower().endswith("met_et"):
        units            = "GeV"
        format           = "%0.0f " + units
        kwargs["xlabel"] = "E_{T}^{miss} (%s)" % units
        kwargs["ylabel"] = yLabel + "/ %s " % format
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 300.0}
        kwargs["log"]    = True
    elif histo.lower().endswith("_n"):
        units            = ""
        format           = "%0.0f " + units
        kwargs["xlabel"] = "multiplicity"
        kwargs["ylabel"] = yLabel + " / %.1f " +  units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +16.0}
    elif histo.lower().endswith("_pt"):
        units            = "GeV/c"
        format           = " / %0.0f " + units
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = yLabel + format
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +400.0}
        if "1" in histo:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +800.0}
        elif "2" in histo:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +500.0}
        elif "3" in histo:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +400.0}
        elif "4" in histo:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +300.0}
        elif "5" in histo:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +250.0}
        elif "6" in histo:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +200.0}
        else:
            pass

        if "gbb_BQuark" in histo:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +300.0}
        if "gtt_tbW_Wqq" in histo:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +500.0}
        if "tbH_BQuark" in histo or "tbH_HPlus" in histo or "tbH_TQuark" in histo:
            kwargs["log"]    = True
            
    elif histo.lower().endswith("_eta"):
        units            = ""
        kwargs["xlabel"] = "#eta"
        kwargs["ylabel"] = yLabel + " / %.1f " + units
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": +5.0}
        kwargs["log"]    = True
    elif histo.lower().endswith("_rap"):
        units            = ""
        kwargs["xlabel"] = "#omega"
        kwargs["ylabel"] = yLabel + " / %.1f " + units
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -5.0, "xmax": +5.0}
        kwargs["log"]    = True
    elif histo.lower().endswith("_phi"):
        units            = "rads"
        kwargs["xlabel"] = "#phi (%s)" % units
        kwargs["ylabel"] = yLabel + "/ %.1f " + units
    elif histo.lower().endswith("_mass"):
        units            = "GeV/c^{2}"
        format           = "%0.0f " + units
        kwargs["xlabel"] = "mass (%s)" % units
        kwargs["ylabel"] = yLabel + " / %.1f " +  units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +500.0}
        kwargs["log"]    = True
        if "MaxDiJetMass" in histo:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +3000.0}
            ROOT.gStyle.SetNdivisions(8, "X")
        if "tt_mass" in histo.lower():
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +3000.0}
            ROOT.gStyle.SetNdivisions(8, "X")
    elif histo.lower().endswith("_deta"):        
        units            = ""
        kwargs["xlabel"] = "#Delta#eta"
        kwargs["ylabel"] = yLabel + " / %.1f " + units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +10.0}
        kwargs["log"]    = True
    elif histo.lower().endswith("_dphi"):
        units            = "rads"
        kwargs["xlabel"] = "#Delta#phi"
        kwargs["ylabel"] = yLabel + " / %.1f " + units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +3.2}
        kwargs["log"]    = False
    elif  histo.lower().endswith("_drap"):
        units            = ""
        kwargs["xlabel"] = "#Delta#omega"
        kwargs["ylabel"] = yLabel + " / %.1f " + units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +10.0}
        kwargs["log"]    = True
    elif histo.lower().endswith("_dr"):
        units            = ""
        kwargs["xlabel"] = "#DeltaR"
        kwargs["ylabel"] = yLabel + " / %.1f " + units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +10.0}
        kwargs["log"]    = True
    elif  histo.lower().endswith("_drrap"):
        units            = ""
        kwargs["xlabel"] = "#DeltaR_{#omega}"
        kwargs["ylabel"] = yLabel + " / %.1f " + units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +10.0}
        kwargs["log"]    = True
    else:
        units = ""
        kwargs["ylabel"] = yLabel + " / %.0f " + units

    # General Options
    if opts.normaliseToOne:
        yMin = 1e-4
    else:
        yMin = 1e0

    if kwargs["log"] == True:
        yMaxFactor = 2 #2chris for eta 16
    else:
        yMaxFactor = 1.2 #1.2 for eta 16 for dphi 1.6

    # Finalise and return
    kwargs["opts"]["ymaxfactor"] = yMaxFactor
    kwargs["opts"]["ymin"]       = yMin
    return kwargs
  
def PlotMC(datasetsMgr, histo, intLumi):

    kwargs = {}
    if opts.normaliseToOne:
        p = plots.MCPlot(datasetsMgr, histo, normalizeToOne=True, saveFormats=[], **kwargs)
    else:
        p = plots.MCPlot(datasetsMgr, histo, normalizeToLumi=intLumi, saveFormats=[], **kwargs)

    # Get histogram<->kwargs dictionary 
    kwargs = GetHistoKwargs(histo, opts)

    # Plot customised histogram
    plots.drawPlot(p, 
                   histo,  
                   xlabel       = kwargs.get("xlabel"),
                   ylabel       = kwargs.get("ylabel"),
                   log          = kwargs.get("log"),
                   rebinX       = kwargs.get("rebinX"), 
                   cmsExtraText = "Preliminary", 
                   #createLegend = {"x1": 0.62, "y1": 0.75, "x2": 0.92, "y2": 0.92},
                   moveLegend   = kwargs.get("moveLegend"),
                   opts         = kwargs.get("opts"),
                   opts2        = {"ymin": 0.6, "ymax": 3}, #chris 1.4
                   cutBox       = kwargs.get("cutBox"),
                   )

    # Customise styling
    if 0:
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kSolid))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerSize(0))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerStyle(6))
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kSolid))
    
    # Save plot in all formats    
    saveName = histo.split("/")[-1]
    if opts.folder == "":
        savePath = os.path.join(opts.saveDir, opts.optMode)
    else:
        savePath = os.path.join(opts.saveDir, histo.split("/")[0], opts.optMode)
    SavePlot(p, saveName, savePath, opts.saveFormats)
    return


def SavePlot(plot, saveName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    Verbose("Saving the plot in %s formats: %s" % (len(saveFormats), ", ".join(saveFormats) ) )
    
    # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
        
    savePath = os.path.join(saveDir, saveName)

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats, 1):
        saveNameURL = savePath + ext
        saveNameURL = saveNameURL.replace("/publicweb/c/cleonido/", "http://home.fnal.gov/~cleonido/")
        if opts.url:
            Print(saveNameURL, i==0)
        else:
            Print(savePath + ext, i==0)
        plot.saveAs(savePath, formats=saveFormats)
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
    ANALYSISNAME = "GenHToHW"
    SEARCHMODE   = "80to1000"
    DATAERA      = "Run2016"
    OPTMODE      = ""
    BATCHMODE    = True
    INTLUMI      = -1.0
    SUBCOUNTERS  = False
    LATEX        = False
    URL          = False
    NOERROR      = True
    SAVEDIR      = "/publicweb/c/cleonido/"  #"/uscms_data/d3/cleonido/workspace/CMSSW_8_0_30/src/HiggsAnalysis/NtupleAnalysis/src/GenHToHW/work" #None #"/publicweb/a/aattikis/" + ANALYSISNAME
    VERBOSE      = False
    HISTOLEVEL   = "Vital" # 'Vital' , 'Informative' , 'Debug'
    NORMALISE    = False
    FOLDER       = "TH1" #"TH2"
    GRIDX        = False
    GRIDY        = False
    SAVEFORMATS  = "png"

    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--analysisName", dest="analysisName", type="string", default=ANALYSISNAME,
                      help="Override default analysisName [default: %s]" % ANALYSISNAME)

    parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
                      help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)
    
    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("--histoLevel", dest="histoLevel", action="store", default = HISTOLEVEL,
                      help="Histogram ambient level (default: %s)" % (HISTOLEVEL))

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("-n", "--normaliseToOne", dest="normaliseToOne", action="store_true", 
                      help="Normalise the histograms to one? [default: %s]" % (NORMALISE) )

    parser.add_option("--folder", dest="folder", type="string", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )

    parser.add_option("-s", "--saveFormats", dest="saveFormats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)     

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX,
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)
    
    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY,
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)      

    (opts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    if opts.mcrab == None:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        parser.print_help()
        #print __doc__
        sys.exit(1)

    # Create save formats
    if "," in opts.saveFormats:
        opts.saveFormats = opts.saveFormats.split(",")
    else:
        opts.saveFormats = [opts.saveFormats]
    opts.saveFormats = ["." + s for s in opts.saveFormats]

    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="")
        Verbose("Save directory set to %s" % (ls + opts.saveDir + ns), True)

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plot_1d.py: Press any key to quit ROOT ...")
