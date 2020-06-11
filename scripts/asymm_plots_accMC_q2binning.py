import ROOT
from ROOT import TCanvas, gStyle, TLegend, TH1
import matplotlib
ROOT.gROOT.SetBatch(True)

matplotlib.use('Agg')
import sys

dataset_files = {
    'flatq2_2016': '/vols/lhcb/amascell/rootFiles/acceptanceMC/flatq2_2016_newBDT.root',
    'try': '/vols/lhcb/amascell/tmp/flatq2_2016_0.root'
}
binning_scheme = {
    # 0: (0.1, 0.98),
    # 1: (1.1, 2.5),
    # 2: (2.5, 4),
    # 3: (4, 6),
    # 4: (6, 8),
    # 5: (15, 17),
    # 6: (17, 19),
    # 7: (11, 12.5),
    '_all': (0.1, 19.0)
}
observables = {'phi': [[-3., 3.], '#phi'],
               'costhetal': [[-1., 1.], 'cos(#theta_{l})'],
               'costhetak': [[-1., 1.], 'cos(#theta_{k})'],
               'B0_MM': [[5170, 5400], 'M_{B^{0}} [MeV]'],
               'Kstar_MM': [[700, 1300], 'M_{K^{*}} [MeV]'],
               }

set_to_run = sys.argv[1]
nbins = 100

print('INFO: Processing sample %s' % set_to_run)
fileName = dataset_files[sys.argv[1]]
treeName = 'DecayTree'
d = ROOT.RDataFrame(treeName, fileName)
d_weights = d.Define('totalWeight', 'track_weight * L0_eff_2016 * Reweights_JpsiK * q2_weight')
TH1.SetDefaultSumw2()
for _bin in binning_scheme:

    bin_range = binning_scheme[_bin]
    cutq2 = 'q2 > {} && q2 < {}'.format(bin_range[0], bin_range[1])
    cutB0 = "B0_ID > 0."
    cutB0bar = "B0_ID < 0."
    cutBDT = 'newBDT > 0.82'
    # observables['q2'] = [[0, 19.], 'q^{2}']
    B0_entries = d_weights.Filter(cutq2).Filter(cutB0).Filter(cutBDT)
    B0bar_entries = d_weights.Filter(cutq2).Filter(cutB0bar).Filter(cutBDT)

    for n, obs in enumerate(observables):
        B0_hist = B0_entries.Histo1D(
            ROOT.RDF.TH1DModel('B0_{}_bin{}'.format(obs, _bin), '%s < q^{2} < %s GeV^{2}' % (bin_range[0], bin_range[1]),
                               nbins, observables[obs][0][0], observables[obs][0][1]), obs, 'totalWeight')
        # B0_hist.Sumw2()
        B0bar_hist = B0bar_entries.Histo1D(
            ROOT.RDF.TH1DModel('B0bar_{}_bin{}'.format(obs, _bin), '%s < q^{2} < %s GeV^{2}' % (bin_range[0], bin_range[1]),
                               nbins, observables[obs][0][0], observables[obs][0][1]), obs, 'totalWeight')
        # B0bar_hist.Sumw2()

        B0_hist = B0_hist.Clone()
        B0bar_hist = B0bar_hist.Clone()

        B0_hist.Scale(1. / B0_hist.Integral())
        B0bar_hist.Scale(1. / B0bar_hist.Integral())
        B0_hist.GetYaxis().SetTitle('frequency')
        B0_hist.GetXaxis().SetTitle(obs)
        B0bar_hist.GetYaxis().SetTitle('frequency')
        B0bar_hist.GetXaxis().SetTitle(obs)
        if obs=='phi' or obs=='costhetal' or obs=='costhetak':
            y_B0 = B0_hist.GetYaxis()
            low_lim = 0.0001
            up_lim = (1. / nbins) * 2.5
            y_B0.SetRangeUser(low_lim, up_lim)
        # B0_hist.SetStats(0)
        B0bar_hist.SetStats(0)
        B0_hist.SetLineColor(4)
        B0bar_hist.SetLineColor(2)

        gStyle.SetOptStat(0)
        gStyle.SetOptFit(1111)
        gStyle.SetStatW(0.1)
        gStyle.SetStatH(0.1)

        ratio_hist = B0_hist.Clone('ratio_{}_bin'.format(obs, _bin))
        ratio_hist.Divide(B0bar_hist)
        ratio_hist.SetTitle('')
        ratio_hist.SetLineColor(1)
        y = ratio_hist.GetYaxis()
        y.SetTitle('ratio B^{0}/#bar{B^{0}}')
        y.SetRangeUser(0., 2.)
        x = ratio_hist.GetXaxis()
        x.SetTitle(observables[obs][1])
        ratio_hist.Fit('pol0')

        # Draw histograms and legend
        drawCanv = TCanvas('obs{}_bin{}'.format(obs, _bin), '{}_bin{}'.format(obs, _bin), 600, 600)
        drawCanv.Divide(1, 2)
        drawCanv.cd(1)
        B0_hist.Draw('hist')
        B0bar_hist.Draw('hist same')
        leg = TLegend(.75, .7, 0.95, .9, '')
        leg.SetFillColor(0)
        leg.SetTextSize(0.04)
        leg.AddEntry(B0_hist, 'B^{0} final states', 'L')
        leg.AddEntry(B0bar_hist, '#bar{B^{0}} final states', 'L')
        leg.Draw('same')
        drawCanv.cd(2)
        ratio_hist.Draw('E')
        drawCanv.GetPad(1).SetMargin(0.1, 0.05, 0.0, 0.1)
        drawCanv.GetPad(2).SetMargin(0.1, 0.05, 0.1, 0.01)

        if n == 0:
            drawCanv.Print('/vols/lhcb/amascell/asymmetry_plots/q2bins_allweights/bin{}.pdf('.format(_bin))
        else:
            if n == len(observables) - 1:
                drawCanv.Print('/vols/lhcb/amascell/asymmetry_plots/q2bins_allweights/bin{}.pdf)'.format(_bin))
            else:
                drawCanv.Print('/vols/lhcb/amascell/asymmetry_plots/q2bins_allweights/bin{}.pdf'.format(_bin))
