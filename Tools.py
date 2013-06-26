import ROOT
import os
from ROOT import kRed, kBlue

def CanvasMaker(histogram1, histogram2, legend, norm=False, save=False, Canvasname='', histogram_Title='', Text=''):
    if histogram1.Integral() > 0 and norm == True:
        histogram1.Scale(1/histogram1.Integral())
    if histogram2.Integral() > 0 and norm == True:
        histogram2.Scale(1/histogram2.Integral())
    canvas = ROOT.TCanvas(Canvasname, Canvasname, 800, 600)
    canvas.cd()
    histogram1.SetTitle(histogram_Title)
    histogram2.SetTitle(histogram_Title)
    if histogram1.GetMaximum() > histogram2.GetMaximum():
        histogram1.Draw()
        if histogram2.Integral() > 0:
            histogram2.Draw("same")
    else:
        histogram2.Draw()
        if histogram1.Integral() > 0:
            histogram1.Draw("same")
    if Text != '':
        CanvasText=ROOT.TPaveText(0.5,0.6,0.9,0.7,"NDC")
        CanvasText.AddText(Text)
        CanvasText.Draw()
    legend.Draw()
    if save == True:
        canvas.SaveAs("{0}.png".format(histogram_Title))
