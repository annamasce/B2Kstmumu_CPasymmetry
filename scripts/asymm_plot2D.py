import ROOT
from ROOT import TCanvas, gStyle, TLegend, TH1, TFile
import matplotlib
import array

ROOT.gROOT.SetBatch(True)

matplotlib.use('Agg')
import sys

dataset_files = {
    'flatq2_2016': '/vols/lhcb/amascell/rootFiles/acceptanceMC/flatq2_2016_newBDT.root',
    'try': '/vols/lhcb/amascell/tmp/flatq2_2016_0.root',
    'Jpsi_run1': '/vols/lhcb/amascell/rootFiles/data/B2KstJpsi_run1_sWeight_newBDT.root',
    'Jpsi_run2': '/vols/lhcb/amascell/rootFiles/data/B2KstJpsi_run2_sWeight_newBDT.root',
    'sideband_run1': '/vols/lhcb/amascell/rootFiles/data/B2Kstmumu_sideband_run1.root',
    'sideband_run2': '/vols/lhcb/amascell/rootFiles/data/B2Kstmumu_sideband_run2.root'
}
# observables = {'K': {'PX': ['K_PX', [-8000., 8000.], 'P_{X}(K) [MeV/c]'],
#                      'PZ': ['K_PZ', [0., 100000.], 'P_{Z}(K) [Mev/c]']},
#
#                'Pi': {'PX': ['Pi_PX', [-8000., 8000.], 'P_{X}(#pi) [MeV/c]'],
#                       'PZ': ['Pi_PZ', [0., 100000.], 'P_{Z}(#pi) [MeV/c]']},
#                }

observables = {'K': {'P': ['K_P', [100, 200000.], 'P(K) [MeV/c]'],
                     'ETA': ['K_ETA', [1.7, 5.0], '#eta(K)']},
               }

polarity = {1: ['Positive', 'Polarity>0.'],
            -1: ['Negative', 'Polarity<0.'],
            0: ['Average', 'Polarity>0 || Polarity<0.']}

set_to_run = sys.argv[1]
nbins = 20

print('INFO: Processing sample %s' % set_to_run)
fileName = dataset_files[sys.argv[1]]
treeName = 'DecayTree'
d = ROOT.RDataFrame(treeName, fileName)
if set_to_run in ['flatq2_2016', 'try']:
    print('Applying MC weights')
    d_weights = d.Define('weight', 'track_weight * L0_eff_2016 * Reweights_JpsiK * q2_weight')
elif 'Jpsi' in set_to_run:
    print('Applying Jpsi weights')
    d_weights = d.Define('weight', 'sWeight')
else:
    d_weights = d

dict_K = d_weights.AsNumpy(['K_P'])
print(dict_K, type(dict_K))

TH1.SetDefaultSumw2()

cutB0 = "B0_ID > 0."
cutB0bar = "B0_ID < 0."
B0_entries = d_weights.Filter(cutB0)
B0bar_entries = d_weights.Filter(cutB0bar)

for obs in observables:
    par_dict = observables[obs]

    print('Generating {} histograms'.format(obs))
    obs1 = 'P'
    obs2 = 'ETA'
    # Split by magnet polarity
    for n, pol in enumerate(polarity):
        print('Polarity {}'.format(polarity[pol][0]))
        B0_entries_pol = B0_entries.Filter(polarity[pol][1])
        B0bar_entries_pol = B0bar_entries.Filter(polarity[pol][1])
        if 'sideband' not in set_to_run:
            B0_hist = B0_entries_pol.Histo2D(
                ROOT.RDF.TH2DModel('B0_{}_pol{}'.format(obs, pol), '',
                                   nbins, par_dict[obs1][1][0], par_dict[obs1][1][1], nbins, par_dict[obs2][1][0],
                                   par_dict[obs2][1][1]), par_dict[obs1][0], par_dict[obs2][0], 'weight')
            B0bar_hist = B0bar_entries_pol.Histo2D(
                ROOT.RDF.TH2DModel('B0bar_{}_pol{}'.format(obs, pol), '',
                                   nbins, par_dict[obs1][1][0], par_dict[obs1][1][1], nbins, par_dict[obs2][1][0],
                                   par_dict[obs2][1][1]), par_dict[obs1][0], par_dict[obs2][0], 'weight')
        else:
            B0_hist = B0_entries_pol.Histo2D(
                ROOT.RDF.TH2DModel('B0_{}_pol{}'.format(obs, pol), '',
                                   nbins, par_dict[obs1][1][0], par_dict[obs1][1][1], nbins, par_dict[obs2][1][0],
                                   par_dict[obs2][1][1]), par_dict[obs1][0], par_dict[obs2][0])
            B0bar_hist = B0bar_entries_pol.Histo2D(
                ROOT.RDF.TH2DModel('B0bar_{}_pol{}'.format(obs, pol), '',
                                   nbins, par_dict[obs1][1][0], par_dict[obs1][1][1], nbins, par_dict[obs2][1][0],
                                   par_dict[obs2][1][1]), par_dict[obs1][0], par_dict[obs2][0])
        B0_hist = B0_hist.Clone()
        B0bar_hist = B0bar_hist.Clone()

        B0_hist.Scale(1. / B0_hist.Integral())
        B0bar_hist.Scale(1. / B0bar_hist.Integral())

        B0_hist.GetYaxis().SetTitle(par_dict[obs2][2])
        B0_hist.GetXaxis().SetTitle(par_dict[obs1][2])
        B0_hist.SetStats(0)

        diff_hist = B0_hist.Clone('diff_{}_pol{}'.format(obs, pol))
        diff_hist.Add(B0bar_hist, -1.)
        sum_hist = B0_hist.Clone('sum_{}_pol{}'.format(obs, pol))
        sum_hist.Add(B0bar_hist, 1.)

        ratio_hist = diff_hist.Clone('ratio_{}_pol{}'.format(obs, pol))
        ratio_hist.Divide(sum_hist)
        ratio_hist.SetTitle(polarity[pol][0] + ' Polarity')

        # Write raw asymmetry histogram to file (only for average polarity)
        if pol==0:
            f = TFile('/vols/lhcb/amascell/asymmetry_plots/2Dplots/track_asymm/{}_{}_plots.root'.format(set_to_run, obs), 'RECREATE')
            B0_hist.Write()
            B0bar_hist.Write()
            ratio_hist.Write()
            diff_hist.Write()
            sum_hist.Write()

        ratio_hist.SetMinimum(-1.)
        ratio_hist.SetMaximum(+1.)

        # Draw histograms
        drawCanv = TCanvas('canvas_{}'.format(obs), '', 600, 600)
        ratio_hist.Draw('COLZ')
        drawCanv.GetPad(0).SetMargin(0.15, 0.15, 0.1, 0.1)

        if n==0:
            drawCanv.Print('/vols/lhcb/amascell/asymmetry_plots/2Dplots/track_asymm/allpol_{}_{}.pdf('.format(set_to_run, obs))
        elif n== len(polarity)-1:
            drawCanv.Print('/vols/lhcb/amascell/asymmetry_plots/2Dplots/track_asymm/allpol_{}_{}.pdf)'.format(set_to_run, obs))
        else:
            drawCanv.Print('/vols/lhcb/amascell/asymmetry_plots/2Dplots/track_asymm/allpol_{}_{}.pdf'.format(set_to_run, obs))