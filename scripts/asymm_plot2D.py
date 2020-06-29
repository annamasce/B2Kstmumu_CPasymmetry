import ROOT
from ROOT import TCanvas, gStyle, TLegend, TH1
import matplotlib

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
binning_scheme = {
    0: (0.1, 0.98),
    1: (1.1, 2.5),
    2: (2.5, 4),
    3: (4, 6),
    4: (6, 8),
    5: (15, 17),
    6: (17, 19),
    7: (11, 12.5),
    '_all': (0.1, 19.0)
}
observables = {'K': {'PT': ['K_PT', [0., 12000.], 'K P_{T}'],
               'ETA': ['K_ETA', [1.7, 5.1], 'K #eta']},

               'Pi': {'PT': ['Pi_PT', [0., 10000.], '#pi P_{T}'],
               'ETA': ['Pi_ETA', [1.7, 5.1], '#pi #eta']},
               }

set_to_run = sys.argv[1]
nbins = 50

print('INFO: Processing sample %s' % set_to_run)
fileName = dataset_files[sys.argv[1]]
treeName = 'DecayTree'
d = ROOT.RDataFrame(treeName, fileName)
if set_to_run in ['flatq2_2016', 'try']:
    print('Applying MC weights')
    d_weights = d.Define('totalWeight', 'track_weight * L0_eff_2016 * Reweights_JpsiK * q2_weight')
if 'Jpsi' in set_to_run:
    print('Applying Jpsi weights')
    d_weights = d.Define('totalWeight', 'sWeight')
else:
    d_weights = d
TH1.SetDefaultSumw2()

cutB0 = "B0_ID > 0."
cutB0bar = "B0_ID < 0."
# cutBDT = 'newBDT > 0.82'
B0_entries = d_weights.Filter(cutB0)
B0bar_entries = d_weights.Filter(cutB0bar)

for obs in observables:
    par_dict = observables[obs]

    print('Generating {} histograms'.format(obs))
    obs1 = 'PT'
    obs2 = 'ETA'
    if 'sideband' not in set_to_run:
        B0_hist = B0_entries.Histo2D(
            ROOT.RDF.TH2DModel('B0_{obs}'.format(obs), '',
                               nbins, par_dict[obs1][1][0], par_dict[obs1][1][1], nbins, par_dict[obs2][1][0],
                               par_dict[obs2][1][1]), par_dict[obs1][0], par_dict[obs2][0], 'totalWeight')
        B0bar_hist = B0bar_entries.Histo2D(
            ROOT.RDF.TH2DModel('B0bar_{}'.format(obs), '',
                               nbins, par_dict[obs1][1][0], par_dict[obs1][1][1], nbins, par_dict[obs2][1][0],
                               par_dict[obs2][1][1]), par_dict[obs1][0], par_dict[obs2][0], 'totalWeight')
    else:
        print('Weight = 1 for the sideband')
        B0_hist = B0_entries.Histo2D(
            ROOT.RDF.TH2DModel('B0_{}'.format(obs), '',
                               nbins, par_dict[obs1][1][0], par_dict[obs1][1][1], nbins, par_dict[obs2][1][0],
                               par_dict[obs2][1][1]), par_dict[obs1][0], par_dict[obs2][0])
        B0bar_hist = B0bar_entries.Histo2D(
            ROOT.RDF.TH2DModel('B0bar_{}'.format(obs), '',
                               nbins, par_dict[obs1][1][0], par_dict[obs1][1][1], nbins, par_dict[obs2][1][0],
                               par_dict[obs2][1][1]), par_dict[obs1][0], par_dict[obs2][0])
    B0_hist = B0_hist.Clone()
    B0bar_hist = B0bar_hist.Clone()

    if 'flatq2' in set_to_run:
        B0_hist.Scale(1. / B0_hist.Integral())
        B0bar_hist.Scale(1. / B0bar_hist.Integral())

    B0_hist.GetYaxis().SetTitle(par_dict[obs2][2])
    B0_hist.GetXaxis().SetTitle(par_dict[obs1][2])
    B0_hist.SetStats(0)

    print('Computing difference')
    diff_hist = B0_hist.Clone('diff_{}'.format(obs))
    diff_hist.Add(B0bar_hist, -1)

    # Draw histograms
    drawCanv = TCanvas('canvas_{}'.format(obs), '', 600, 600)
    diff_hist.Draw('COLZ')
    drawCanv.GetPad(0).SetMargin(0.1, 0.15, 0.1, 0.1)

    drawCanv.Print('/vols/lhcb/amascell/asymmetry_plots/2Dplots/{}_{}.pdf'.format(set_to_run, obs))