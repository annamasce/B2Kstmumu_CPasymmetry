import pandas
from scripts.data_processing.data import CP_divide
import ROOT
from ROOT import TH1D, TCanvas, gStyle, TLegend

ROOT.gROOT.SetBatch(True)

binning_scheme = {
        0: (0.1, 0.98),
        1: (1.1, 2.5),
        2: (2.5, 4),
        3: (4, 6),
        4: (6, 8),
        5: (15, 17),
        6: (17, 19),
        7: (11, 12.5)
    }

if __name__ == '__main__':
    dataPath = '/home/anna/master_thesis/data/whole_run2/'
    saveLocation = '/home/anna/master_thesis/files_asymmetry/0526/'
    plotName = '0526_q2binning_angles_pol1'
    fileType = 'sideband'
    years = [2016, 2017, 2018]
    # fileType = 'sideband'

    # data_object = Dataset(dataPath, years=years, fileType=fileType, L0trigger=True)
    # data = data_object.get_data()
    #
    # # plt.figure()
    # # plt.hist((numpy.logical_or(data['B0_L0MuonDecision_TIS']==True, data['B0_L0DiMuonDecision_TIS']==True)).astype('float64'))
    # # plt.show()
    #
    # # data = data[data['Polarity']<0].reset_index(drop=True)
    # # print(data.head(10))
    # # plt.figure()
    # # plt.hist(data['Polarity'])
    # # plt.show()
    #
    # print('Length data after trigger selection:', len(data))
    # data_B0, data_B0bar = CP_divide(data)
    # data_B0['phi'] = 0.
    # print('Adding angular variables to B0 data')
    # add_angvar_todata(data_B0)
    # print('Adding angular variables to B0bar data')
    # add_angvar_todata(data_B0bar)

    # File containing new BDT results and right angular variables already
    bkg = pandas.read_pickle('/home/anna/master_thesis/data/withBDT/B2Kstmumu_sideband_wholeRun2.pkl')
    # bkg = pandas.read_pickle('/home/anna/master_thesis/data/withBDT/B2KstJpsi_wholeRun2_sWeight_wL0.pkl')
    print(bkg.head(5))
    data_B0, data_B0bar = CP_divide(bkg)

    observables = {'phi': [[-3., 3.], '#phi'],
                   'costhetal': [[-1., 1.], 'cos(#theta_{l})'],
                   'costhetak': [[-1., 1.], 'cos(#theta_{k})'],
                   # 'B0_MM': [[5300, 7000], 'M_{B^{0}} [MeV]'],
                   # 'Kstar_MM': [[800, 990], 'M_{K^{*}} [MeV]'],
                   # 'mu_P': [[5000, 150000], 'P_{#mu} [MeV/c]'],
                   'mu_PT': [[0, 10000], 'P_{T#mu} [MeV/c]'],
                   # 'mu_TRACK_CHI2NDOF': [[0.3, 2.9], '#mu #chi^{2}/DOF'],
                   # 'mu_ETA': [[1.7, 5.1], 'eta_{#mu}']
                   }

    if fileType=='Jpsi-sWeight':
        observables['B0_MM'][0] = [5170, 5400]
        observables['Kstar_MM'][0] = [700, 1600]

    for _bin in binning_scheme:

        bin_range = binning_scheme[_bin]
        data_B0_bin = data_B0.loc[(data_B0['q2']>bin_range[0]) & (data_B0['q2']<bin_range[1])].reset_index(drop=True)
        # plt.figure()
        # plt.hist(data_B0_bin['q2'].values, bins=10)
        # plt.show()
        data_B0bar_bin = data_B0bar.loc[(data_B0bar['q2'] > bin_range[0]) & (data_B0bar['q2'] < bin_range[1])].reset_index(drop=True)

        if fileType == 'Jpsi-sWeight':
            weights_B0 = data_B0_bin['sWeight'].values
            weights_B0bar = data_B0bar_bin['sWeight'].values

        for n, obs in enumerate(observables.keys()):

            B0_hist = TH1D('B0_{}_bin{}'.format(obs, _bin), '%s < q^{2} < %s GeV^{2}' %(bin_range[0], bin_range[1]), 30, observables[obs][0][0], observables[obs][0][1])
            B0_hist.Sumw2()
            B0bar_hist = TH1D('B0bar_{}_bin{}'.format(obs, _bin), '%s < q^{2} < %s GeV^{2}' %(bin_range[0], bin_range[1]), 30, observables[obs][0][0], observables[obs][0][1])
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
            if obs=='mu_TRACK_CHI2NDOF':
                obs_B0 = 'mu_plus_TRACK_CHI2NDOF'
                obs_B0bar = 'mu_minus_TRACK_CHI2NDOF'
            if obs=='mu_ETA':
                obs_B0 = 'mu_plus_ETA'
                obs_B0bar = 'mu_minus_ETA'


            for i, event in enumerate(data_B0_bin[obs_B0]):
                if fileType == 'Jpsi-sWeight':
                    B0_hist.Fill(event, weights_B0[i])
                else:
                    B0_hist.Fill(event)

            for i, event in enumerate(data_B0bar_bin[obs_B0bar]):
                if fileType == 'Jpsi-sWeight':
                    B0bar_hist.Fill(event, weights_B0bar[i])
                else:
                    B0bar_hist.Fill(event)

            B0_hist.Scale(1. / B0_hist.Integral(), 'width')
            B0bar_hist.Scale(1. / B0bar_hist.Integral(), 'width')
            B0_hist.GetYaxis().SetTitle('frequency')
            B0_hist.GetXaxis().SetTitle(obs)
            B0bar_hist.GetYaxis().SetTitle('frequency')
            B0bar_hist.GetXaxis().SetTitle(obs)
            # B0_hist.SetStats(0)
            B0bar_hist.SetStats(0)
            B0_hist.SetLineColor(4)
            B0bar_hist.SetLineColor(2)

            gStyle.SetOptStat(0)
            gStyle.SetOptFit(1111)
            gStyle.SetStatW(0.1)
            gStyle.SetStatH(0.1)

            ratio_hist = (B0_hist).Clone('ratio_{}_bin'.format(obs, _bin))
            ratio_hist.Divide(B0bar_hist)
            ratio_hist.SetTitle('')
            ratio_hist.SetLineColor(1)
            y = ratio_hist.GetYaxis()
            y.SetTitle('ratio B^{0}/#bar{B^{0}}')
            y.SetRangeUser(0., 3.)
            x = ratio_hist.GetXaxis()
            x.SetTitle(observables[obs][1])
            ratio_hist.Fit('pol0')
            #ratio_hist.GetYaxis.SetTitle('ratio')

            # Draw histograms and legend
            drawCanv = TCanvas('obs{}_bin{}'.format(obs, _bin), '{}_bin{}'.format(obs, _bin), 600, 600)
            drawCanv.Divide(1, 2)
            drawCanv.cd(1)
            B0bar_hist.Draw('hist')
            B0_hist.Draw('hist same')
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

            if _bin==0 and n==0:
                drawCanv.Print(saveLocation + plotName + '_' + fileType + '.pdf(')
                # drawCanv.Print(saveLocation + plotName + '_' + fileType + '.svg')
            else:
                if _bin==len(binning_scheme)-1 and n == len(observables)-1:
                    drawCanv.Print(saveLocation + plotName + '_' + fileType + '.pdf)')
                else:
                    drawCanv.Print(saveLocation + plotName + '_' + fileType + '.pdf')

