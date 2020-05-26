import numpy
import pandas
import root_numpy
import matplotlib.pyplot as plt
from data import Dataset, CP_divide, L0trigg_selection_TOS
import ROOT
from ROOT import TH1D, TCanvas, TFile, TTree, TPad, TChain, TDirectoryFile, TGraph, TF1, TGraphErrors, TLine, TGaxis, gStyle, TLegend, TLorentzVector
ROOT.gROOT.SetBatch(True)
from ang_functions import add_angvar_todata


saveLocation = '/vols/lhcb/amascell/asymmetry_plots/'
plotName = 'flatq2_2016_100bins_allWeights_newBDT'

# load tree from file
filePath = '/vols/lhcb/amascell/rootFiles/flatq2_2016_newBDT.root'
# filePath = '/vols/lhcb/amascell/tmp/flatq2_2016_0.root'
inFile = TFile(filePath)
inTree = inFile.Get('DecayTree')

observables = {'phi': [[-3., 3.], '#phi'],
               'costhetal': [[-1., 1.], 'cos(#theta_{l})'],
               'costhetak': [[-1., 1.], 'cos(#theta_{k})'],
               'B0_MM': [[5170, 5400], 'M_{B^{0}} [MeV]'],
               'Kstar_MM': [[700, 1300], 'M_{K^{*}} [MeV]'],
               # 'mu_P': [[5000, 150000], 'P_{#mu} [MeV/c]'],
               # 'mu_PT': [[0, 10000], 'P_{T#mu} [MeV/c]'],
               # 'mu_TRACK_CHI2NDOF': [[0.3, 2.9], '#mu #chi^{2}/DOF'],
               # 'mu_ETA': [[1.7, 5.1], 'eta_{#mu}']
               }

B0_hist = {}
B0bar_hist = {}
ratio_hist = {}

# Define B0 and B0bar histograms for each observable
for n, obs in enumerate(observables.keys()):
    B0_hist[obs] = TH1D('B0_hist_{}'.format(obs), '', 100, observables[obs][0][0], observables[obs][0][1])
    B0_hist[obs].Sumw2()
    B0bar_hist[obs] = TH1D('B0bar_hist_{}'.format(obs), '', 100, observables[obs][0][0], observables[obs][0][1])
    B0bar_hist[obs].Sumw2()
    #ratio hist to be done later after filling the histograms

print('Reading the tree...')
for ie in range(0, inTree.GetEntries()):
    inTree.GetEntry(ie)
    if (ie+1)%100000 == 0:
        print('{} events processed'.format(ie+1))
    if getattr(inTree, 'newBDT') > 0.82:
        weight = getattr(inTree, 'track_weight') * getattr(inTree, 'L0_eff_2016') * getattr(inTree, 'Reweights_JpsiK') * getattr(inTree, 'q2_weight')
        for obs in observables.keys():
            obs_B0 = obs
            obs_B0bar = obs
            # if obs=='mu_P':
            #     obs_B0 = 'mu_plus_P'
            #     obs_B0bar = 'mu_minus_P'
            # if obs=='mu_PT':
            #     obs_B0 = 'mu_plus_PT'
            #     obs_B0bar = 'mu_minus_PT'
            # if obs=='mu_TRACK_CHI2NDOF':
            #     obs_B0 = 'mu_plus_TRACK_CHI2NDOF'
            #     obs_B0bar = 'mu_minus_TRACK_CHI2NDOF'
            # if obs=='mu_ETA':
            #     obs_B0 = 'mu_plus_ETA'
            #     obs_B0bar = 'mu_minus_ETA'
            if getattr(inTree, 'B0_ID')>0:
                B0_hist[obs].Fill(getattr(inTree, obs_B0), weight)
            else:
                B0bar_hist[obs].Fill(getattr(inTree, obs_B0bar), weight)

for n, obs in enumerate(observables.keys()):

    B0_hist[obs].Scale(1. / B0_hist[obs].Integral(), 'width')
    B0bar_hist[obs].Scale(1. / B0bar_hist[obs].Integral(), 'width')
    B0_hist[obs].GetYaxis().SetTitle('frequency')
    B0_hist[obs].GetXaxis().SetTitle(obs)
    B0bar_hist[obs].GetYaxis().SetTitle('frequency')
    B0bar_hist[obs].GetXaxis().SetTitle(obs)
    # B0_hist[obs].SetStats(0)
    B0bar_hist[obs].SetStats(0)
    B0_hist[obs].SetLineColor(4)
    B0bar_hist[obs].SetLineColor(2)

    gStyle.SetOptStat(0)
    gStyle.SetOptFit(1111)
    gStyle.SetStatW(0.1)
    gStyle.SetStatH(0.1)

    ratio_hist[obs] = (B0_hist[obs]).Clone('ratio_hist')
    ratio_hist[obs].Divide(B0bar_hist[obs])
    ratio_hist[obs].SetTitle('')
    ratio_hist[obs].SetLineColor(1)
    y = ratio_hist[obs].GetYaxis()
    y.SetTitle('ratio B^{0}/#bar{B^{0}}')
    y.SetRangeUser(0.5, 1.5)
    x = ratio_hist[obs].GetXaxis()
    x.SetTitle(observables[obs][1])
    ratio_hist[obs].Fit('pol0')
    #ratio_hist.GetYaxis.SetTitle('ratio')

    # Draw histograms and legend
    drawCanv = TCanvas('obs_{}'.format(obs), '{}'.format(obs), 600, 600)
    drawCanv.Divide(1, 2)
    drawCanv.cd(1)
    B0_hist[obs].Draw('hist')
    B0bar_hist[obs].Draw('hist same')
    leg = TLegend(.75, .7, 0.95, .9, '')
    leg.SetFillColor(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(B0_hist[obs], 'B^{0} final states', 'L')
    leg.AddEntry(B0bar_hist[obs], '#bar{B^{0}} final states', 'L')
    leg.Draw('same')
    drawCanv.cd(2)
    ratio_hist[obs].Draw('E')
    drawCanv.GetPad(1).SetMargin(0.1, 0.05, 0.0, 0.1)
    drawCanv.GetPad(2).SetMargin(0.1, 0.05, 0.1, 0.01)

    # Save canvas to PDF file
    if len(observables)==1:
        drawCanv.Print(saveLocation + plotName + '.pdf')
    elif n==0:
        drawCanv.Print(saveLocation + plotName + '.pdf(')
    else:
        if n == len(observables)-1:
            drawCanv.Print(saveLocation + plotName + '.pdf)')
        else:
            drawCanv.Print(saveLocation + plotName + '.pdf')

