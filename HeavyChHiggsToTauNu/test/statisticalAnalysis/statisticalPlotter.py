from HiggsAnalysis.HeavyChHiggsToTauNu.tools.statisticalFunctions import *

import ROOT
from ROOT import TTree, gROOT, TGraph, TCanvas, TMultiGraph, TLegend, TAxis, TLatex
from HiggsAnalysis.HeavyChHiggsToTauNu.tools.histograms import *
from HiggsAnalysis.HeavyChHiggsToTauNu.tools.tdrstyle import *
import HiggsAnalysis.HeavyChHiggsToTauNu.tools.styles as styles
from array import array

### Apply TDR style
style = TDRStyle()
gROOT.Reset()
ROOT.gROOT.SetBatch(True)

def getLegend(x):
#    lege = TLegend(x, 0.7, 0.95, 0.93)
    lege = TLegend(0.2, 0.48, 0.55, 0.68)
    lege.SetFillStyle(0)
    lege.SetBorderSize(0)
    lege.SetTextSize(0.03)
    return lege

def writeText( myText, y ):
    text = TLatex()
    text.SetTextAlign(12)
    text.SetTextSize(0.04)
    text.SetNDC()
    text.DrawLatex(0.2,y,myText)
    return 0

def fillDataTeva( ):
    ## Tevatron 1/fb results, D0
    dataTevaMhMax = {
        "mass":array( 'd' ),
        "tanb":array( 'd' )
        }
    dataTevaMhMax["mass"].append(90)
    dataTevaMhMax["tanb"].append(30.8)
    dataTevaMhMax["mass"].append(100)
    dataTevaMhMax["tanb"].append(33.5)
    dataTevaMhMax["mass"].append(120)
    dataTevaMhMax["tanb"].append(50)
    dataTevaMhMax["mass"].append(140)
    dataTevaMhMax["tanb"].append(103.5)
    dataTevaMhMax["mass"].append(145)
    dataTevaMhMax["tanb"].append(147.6)
    return dataTevaMhMax

def getGraphTevatron():
    dataTevaMhMax = fillDataTeva()
    myGraph = TGraph(len(dataTevaMhMax["mass"]),dataTevaMhMax["mass"],dataTevaMhMax["tanb"])
    myGraph.SetLineStyle(2)
    myGraph.SetLineWidth(2)
    return myGraph
    
def main():

    luminosity = 100

    massPoints = {
        "PFTauCutBased": {
	    90:  MassPoint(luminosity*0.1359,20,luminosity*(0.3601+0.2010)),
	    100: MassPoint(luminosity*0.1262,20,luminosity*(0.3601+0.2010)),
	    120: MassPoint(luminosity*0.0943,20,luminosity*(0.3601+0.2010)),
	    140: MassPoint(luminosity*0.0381,20,luminosity*(0.3601+0.2010)),
	    160: MassPoint(luminosity*0.00833,20,luminosity*(0.3601+0.2010))
	},
	"PFTauTaNCBased": {
            90:  MassPoint(luminosity*0.148,20,luminosity*(0.1297+0.2192)),
            100: MassPoint(luminosity*0.1349,20,luminosity*(0.1297+0.2192)),
            120: MassPoint(luminosity*0.1015,20,luminosity*(0.1297+0.2192)),
            140: MassPoint(luminosity*0.0520,20,luminosity*(0.1297+0.2192)),
            160: MassPoint(luminosity*0.0095,20,luminosity*(0.1297+0.2192))
	}
    }

    mus  = [-1000,-200,200,1000] # mu parameters
    mHps  = massPoints["PFTauCutBased"].keys() # H+ masses
    mHps.sort()
#    mus  = [200]
#    mHps = [120]

    data = {}
        
    nSigma = 5
    clSigma = 1.95996
    sysError = 0.1
    for selection in massPoints.keys() :
	print selection
        data[selection] = {}
	for mu in mus :
	    print "mu = ",mu
	    tanbExclNoErr   = []
	    tanbExclWErr    = []
	    tanbReachNoErr  = []
	    tanbReachWErr   = []
	    tanbReachTheory = []
            nPoints = 5
            data[selection][mu] = {}
            data[selection][mu]["ExclNoErr"] = {
                "mass":array( 'd' ),
                "tanb":array( 'd' )
                }
            data[selection][mu]["ExclWErr"] = {
                "mass":array( 'd' ),
                "tanb":array( 'd' )
                }
            data[selection][mu]["ReachNoErr"] = {
                "mass":array( 'd' ),
                "tanb":array( 'd' )
                }
            data[selection][mu]["ReachWErr"] = {
                "mass":array( 'd' ),
                "tanb":array( 'd' )
                }
            data[selection][mu]["ReachTheory"] = {
                "mass":array( 'd' ),
                "tanb":array( 'd' )
                }
            massArray, tanbExclNoErrArray, tanbExclWErrArray, tanbReachNoErrArray, tanbReachWErrArray, tanbReachTheoryArray = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ),array( 'd' )
	    for mass in mHps :
		tanbTheoryReach = tanbForTheoryLimit(mass,mu)
		tanbReachTheory.append(tanbTheoryReach)
#		print massPoints[selection][mass].nSignal,massPoints[selection][mass].tanbRef,massPoints[selection][mass].nBackgr,mass,mu
#		print signif(massPoints[selection][mass].nSignal,massPoints[selection][mass].nBackgr,0),signalAtNsigma(massPoints[selection][mass].nBackgr,0,5)
		tanbAt5sigmaNoErr = tanbAtNsigma(massPoints[selection][mass].nSignal,massPoints[selection][mass].tanbRef,massPoints[selection][mass].nBackgr,0,mass,mu,nSigma)
		tanbAt5sigmaWErr = tanbAtNsigma(massPoints[selection][mass].nSignal,massPoints[selection][mass].tanbRef,massPoints[selection][mass].nBackgr,sysError,mass,mu,nSigma)
		tanbAt95CLNoErr = tanbAtNsigma(massPoints[selection][mass].nSignal,massPoints[selection][mass].tanbRef,massPoints[selection][mass].nBackgr,0,mass,mu,clSigma)
		tanbAt95CLWErr = tanbAtNsigma(massPoints[selection][mass].nSignal,massPoints[selection][mass].tanbRef,massPoints[selection][mass].nBackgr,sysError,mass,mu,clSigma)
		print "mass,th-reach,reach,reach(sys),excl,excl(sys) = ",mass,tanbTheoryReach,tanbAt5sigmaNoErr,tanbAt5sigmaWErr,tanbAt95CLNoErr,tanbAt95CLWErr
		tanbReachNoErr.append(tanbAt5sigmaNoErr)
		tanbReachWErr.append(tanbAt5sigmaWErr)
		tanbExclNoErr.append(tanbAt95CLNoErr)
		tanbExclWErr.append(tanbAt95CLWErr)
                massArray.append(mass)
                tanbExclNoErrArray.append(tanbAt5sigmaNoErr)
                tanbExclWErrArray.append(tanbAt5sigmaWErr)
                tanbReachNoErrArray.append(tanbAt95CLNoErr)
                tanbReachWErrArray.append(tanbAt95CLWErr)
                tanbReachTheoryArray.append(tanbTheoryReach)

                if tanbAt5sigmaNoErr>0:
                    data[selection][mu]["ExclNoErr"]["mass"].append(mass)
                    data[selection][mu]["ExclNoErr"]["tanb"].append(tanbAt5sigmaNoErr)
                if tanbAt5sigmaWErr>0:
                    data[selection][mu]["ExclWErr"]["mass"].append(mass)
                    data[selection][mu]["ExclWErr"]["tanb"].append(tanbAt5sigmaWErr)
                if tanbAt95CLNoErr>0:    
                    data[selection][mu]["ReachNoErr"]["mass"].append(mass)
                    data[selection][mu]["ReachNoErr"]["tanb"].append(tanbAt95CLNoErr)
                if tanbAt95CLWErr>0:
                    data[selection][mu]["ReachWErr"]["mass"].append(mass)
                    data[selection][mu]["ReachWErr"]["tanb"].append(tanbAt95CLWErr)
                if tanbTheoryReach>0:
                    data[selection][mu]["ReachTheory"]["mass"].append(mass)
                    data[selection][mu]["ReachTheory"]["tanb"].append(tanbTheoryReach)
                    
# Plot everything in one canvas
#    c1 = TCanvas( 'c1', 'tanbReach', 200, 10, 700, 500 )
#    graphExclNoErr   = TGraph(nPoints,massArray,tanbExclNoErrArray)
#    graphExclWErr    = TGraph(nPoints,massArray,tanbExclWErrArray)
#    graphReachNoErr  = TGraph(nPoints,massArray,tanbReachNoErrArray)
#    graphReachWErr   = TGraph(nPoints,massArray,tanbReachWErrArray)
#    graphReachTheory = TGraph(nPoints,massArray,tanbReachTheoryArray)
#    graphExclNoErr.Draw('ALP')   
#    graphExclWErr.Draw('LP')    
#    graphReachNoErr.Draw('LP')      
#    graphReachWErr.Draw('LP')       
#    graphReachTheory.Draw('LP')
#    addCmsPreliminaryText()
#    c1.Update()
#    c1.SaveAs(".png")

    ## Graph: different mu values
    ## plot theory reach & 5sigmaNoErr & Tevatron exclusion
    for selection in massPoints.keys() :
        color = 1
        print 'muPlot'+selection
        c1 = TCanvas( 'muPlot'+selection, 'tanbReach'+selection, 200, 10, 700, 500 )
        multi = TMultiGraph()
        lege = getLegend(0.6)
        writeText("L = "+str(luminosity)+" pb^{-1}", 0.9)
        writeText("m_{H}^{max} scenario",0.84)
        writeText("t#rightarrowbH#pm#rightarrowb#tau#nu#rightarrowhadrons + #nu", 0.78)
        writeText("no syst. errors", 0.72)
        for mu in mus :
            yvalues, massArray = array( 'd' ), array( 'd' )
            print "mu = ",mu
            for mass in mHps :
                tanbTheoryReach = tanbForTheoryLimit(mass,mu)
                tanbAt5sigmaNoErr = tanbAtNsigma(massPoints[selection][mass].nSignal,massPoints[selection][mass].tanbRef,massPoints[selection][mass].nBackgr,0,mass,mu,nSigma)
                # Skip points with tanB = -1 
                if tanbAt5sigmaNoErr != -1 :
                    yvalues.append( tanbAt5sigmaNoErr )
                    massArray.append(mass)
                print "graph1, theory & 5sigma & mass: ",tanbTheoryReach,tanbAt5sigmaNoErr, mass
            if len(massArray)>0:
                graphMus = TGraph(len(massArray),massArray,yvalues)
                graphMus.SetLineColor(color)
                graphMus.SetMarkerColor(color)
                multi.Add(graphMus,"lp")
                lege.AddEntry(graphMus,"5 sigma, mu = "+str(mu),"l")
            color = color + 1
        # Tevatron result
        graphTeva = getGraphTevatron()
        multi.Add(graphTeva,"lp")
        lege.AddEntry(graphTeva,"Tevatron 1fb^{-1} exclusion","l")
        lege.Draw()
        multi.Draw("a")
        c1.Update()
        multi.GetYaxis().SetRangeUser(0,200)
        multi.GetYaxis().SetTitle("tan(#beta)")
        multi.GetXaxis().SetTitle("M_{H^{#pm}} [GeV/c^{2}]")
        addCmsPreliminaryText()
        c1.SaveAs(".png")

## Plot comparison of reach and exclusion for 2 selections
## fix mu=200       
## make one plot without errors, one with errors
    setOfSetOfGraphs = [ ["ExclNoErr","ReachNoErr"],
                         ["ExclWErr","ReachWErr"]]
    for index, setOfGraphs in enumerate(setOfSetOfGraphs):
        c1 = TCanvas("exclusions"+str(index),"exclusions"+str(index),200, 10, 700, 500 )
        multi = TMultiGraph()
        lege = getLegend(0.5)
        mymu = 200
        color = 1
        for selection in massPoints.keys() :
            for oneGraph in setOfGraphs:
                graph = TGraph(len(data[selection][mymu][oneGraph]["mass"]),
                               data[selection][mymu][oneGraph]["mass"],
                               data[selection][mymu][oneGraph]["tanb"])
                graph.SetLineColor(color)
                graph.SetMarkerColor(color)
                lege.AddEntry(graph,selection+", "+oneGraph,"l")
                multi.Add(graph,"lp")
                color = color + 1
        graphTeva = getGraphTevatron()
        multi.Add(graphTeva,"lp")
        lege.AddEntry(graphTeva,"Tevatron 1fb^{-1} exclusion","l")
        multi.Draw("a")
        lege.Draw()
        addCmsPreliminaryText()
        writeText("L = "+str(luminosity)+" pb^{-1}", 0.9)
        writeText("m_{H}^{max} scenario",0.84)
        writeText("t#rightarrowbH#pm#rightarrowb#tau#nu#rightarrowhadrons + #nu", 0.78)
        multi.GetYaxis().SetRangeUser(0,200)
        multi.GetYaxis().SetTitle("tan(#beta)")
        multi.GetXaxis().SetTitle("M_{H^{#pm}} [GeV/c^{2}]")
        c1.SaveAs(".png")

    ############################### EXECUTION ###############################
    ### Script execution can be paused like this, it will continue after
    ### user has given some input (which must include enter)
#    raw_input("Hit enter to continue") ### keep canvas open until you hit enter

#    print signif(50,50,0)
#    print signalAtNsigma(50,0,5)



main()
