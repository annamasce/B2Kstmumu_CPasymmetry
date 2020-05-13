import numpy
import pandas
import root_numpy
import matplotlib.pyplot as plt
from data import Dataset
import ROOT
from ROOT import TH1D, TCanvas, TFile, TTree, TPad, TChain, TDirectoryFile, TGraph, TF1, TGraphErrors, TLine, TGaxis, gStyle, TLegend
ROOT.gROOT.SetBatch(True)
from ang_functions import add_angvar_todata


def make_single_plot(data, var):
    plt.figure()
    plt.hist(data[var], weights=data['sWeight'])
    plt.show()


if __name__ == '__main__':
    dataPath = '/home/anna/master_thesis/data/whole_run2/'
    saveLocation = '/home/anna/master_thesis/files_asymmetry/'
    plotName = '0508_wholeRun2'
    fileType = 'sideband'
    years = [2016, 2017, 2018]
    # fileType = 'sideband'

    data_object = Dataset(dataPath, years=years, fileType=fileType, L0trigger=True)
    # data = data_object.get_data()
    print('Length data after trigger selection:', len(data_object.data))
    data_B0, data_B0bar = data_object.CP_divide()
    data_B0['phi'] = 0.
    print('Adding angular variables to B0 data')
    add_angvar_todata(data_B0)
    print('Adding angular variables to B0bar data')
    add_angvar_todata(data_B0bar)

    if fileType=='Jpsi-sWeight':
        weights_B0 = data_B0['sWeight'].values
        weights_B0bar = data_B0bar['sWeight'].values

    observables = {'phi': [[-3., 3.], '#phi'],
                   'costhetal': [[-1., 1.], 'cos(#theta_{l})'],
                   'costhetak': [[-1., 1.], 'cos(#theta_{k})'],
                   'B0_MM': [[5300, 7000], 'M_{B^{0}} [MeV]'],
                   'Kstar_MM': [[800, 990], 'M_{K^{*}} [MeV]'],
                   'mu_P': [[5000, 140000], 'P_{#mu} [MeV/c]'],
                   'mu_PT': [[0, 10000], 'P_{T#mu} [MeV/c]']
                   }

    if fileType=='Jpsi-sWeight':
        observables['B0_MM'][0] = [5170, 5400]
        observables['Kstar_MM'][0] = [700, 1600]

    for n, obs in enumerate(observables.keys()):

        B0_hist = TH1D('B0_hist_{}'.format(obs), '', 50, observables[obs][0][0], observables[obs][0][1])
        B0_hist.Sumw2()
        B0bar_hist = TH1D('B0bar_hist_{}'.format(obs), '', 50, observables[obs][0][0], observables[obs][0][1])
        B0bar_hist.Sumw2()
        #ratio hist to be done later after filling the histograms

        obs_B0 = obs
        obs_B0bar = obs
        if obs=='mu_P':
            obs_B0 = 'mu_plus_P'
            obs_B0bar = 'mu_minus_P'
        if obs=='mu_PT':
            obs_B0 = 'mu_plus_PT'
            obs_B0bar = 'mu_minus_PT'

        for i, event in enumerate(data_B0[obs_B0]):
            if fileType == 'Jpsi-sWeight':
                B0_hist.Fill(event, weights_B0[i])
            else:
                B0_hist.Fill(event)

        for i, event in enumerate(data_B0bar[obs_B0bar]):
            if fileType == 'Jpsi-sWeight':
                B0bar_hist.Fill(event, weights_B0bar[i])
            else:
                B0bar_hist.Fill(event)

        B0_hist.GetYaxis().SetTitle('counts')
        B0_hist.GetXaxis().SetTitle(obs)
        # B0_hist.SetStats(0)
        B0bar_hist.SetStats(0)
        B0_hist.SetLineColor(4)
        B0bar_hist.SetLineColor(2)

        gStyle.SetOptStat(0)
        gStyle.SetOptFit(1111)
        gStyle.SetStatW(0.1)
        gStyle.SetStatH(0.1)

        ratio_hist = (B0_hist).Clone('ratio_hist')
        ratio_hist.Divide(B0bar_hist)
        ratio_hist.SetTitle('')
        ratio_hist.SetLineColor(1)
        y = ratio_hist.GetYaxis()
        y.SetTitle('ratio B^{0}/#bar{B^{0}}')
        y.SetRangeUser(0.5, 1.5)
        x = ratio_hist.GetXaxis()
        x.SetTitle(observables[obs][1])
        ratio_hist.Fit('pol0')
        #ratio_hist.GetYaxis.SetTitle('ratio')

        # Draw histograms and legend
        drawCanv = TCanvas('obs_{}'.format(obs), '{}'.format(obs), 600, 600)
        drawCanv.Divide(1, 2)
        drawCanv.cd(1)
        B0bar_hist.Draw('hist')
        B0_hist.Draw('hist same')
        leg = TLegend(.75, .7, 0.95, .9, '')
        leg.SetFillColor(0)
        leg.SetTextSize(0.04)
        leg.AddEntry(B0_hist, 'B^{0} final states', 'L')
        leg.AddEntry(B0bar_hist, '#bar{B^{0}} final states', 'L')
        leg.DrawClone('same')
        drawCanv.cd(2)
        ratio_hist.Draw('E')
        drawCanv.GetPad(1).SetMargin(0.1, 0.05, 0.0, 0.1)
        drawCanv.GetPad(2).SetMargin(0.1, 0.05, 0.1, 0.01)

        if n==0:
            drawCanv.Print(saveLocation + plotName + '_' + fileType + '.pdf(')
            # drawCanv.Print(saveLocation + plotName + '_' + fileType + '.svg')
        else:
            if n == len(observables)-1:
                drawCanv.Print(saveLocation + plotName + '_' + fileType + '.pdf)')
            else:
                drawCanv.Print(saveLocation + plotName + '_' + fileType + '.pdf')

