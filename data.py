import numpy
import pandas
import root_numpy
import matplotlib.pyplot as plt

def L0trigg_selection(data):
    selection = numpy.logical_or(data['B0_L0MuonDecision_TOS'] == 1, data['B0_L0DiMuonDecision_TOS'] == 1)
    selected_data = data.loc[selection].reset_index()
    # print(len(selected_data))
    return selected_data

def MergeDataframes(data):
    '''

    :param data: tuple of dataframes
    :return: merged dataframe
    '''

    data_merged = data[0]
    if len(data)>0:
        for i in range(1, len(data)):
            data_merged = pandas.concat([data_merged, data[i]], ignore_index=True)

    return data_merged

def CP_divide(data):
    data_B0 = data[data['B0_ID']>0]
    data_B0bar = data[data['B0_ID']<0]
    return data_B0, data_B0bar

def make_asymm_plots(data_B0, data_B0bar, var_name):
    fig, ax = plt.subplots(2,1)
    ax[0].plot(data_B0[var_name])
    ax[1].plot(data_B0bar[var_name])

class Dataset(object):

    def __init__(self, dataPath, years=(2016, 2017, 2018), fileType='Jpsi-sWeight', L0trigger=True):
        # Initialization
        self.path = dataPath
        self.years = years
        self.type = fileType
        self.L0trigger = L0trigger
        self.data = self.get_data()
        # self.length = length
    def get_data(self):
        data = []
        if self.type == 'Jpsi-sWeight':
            for year in self.years:
                filePath = self.path + 'B2KstJpsi_' + str(year) + '_sWeight_wL0.root'
                data.append(pandas.DataFrame(root_numpy.root2array(filePath, treename='DecayTree')))
        if self.type == 'sideband':
            for year in self.years:
                filePath = self.path + 'B2Kstmumu_sideband_' + str(year) + '.root'
                data.append(pandas.DataFrame(root_numpy.root2array(filePath, treename='DecayTree')))
        if self.L0trigger:
            return L0trigg_selection(MergeDataframes(data))
        else:
            return MergeDataframes(data)

    def CP_divide(self):
        data_B0 = self.data[self.data['B0_ID'] > 0].reset_index()
        data_B0bar = self.data[self.data['B0_ID'] < 0].reset_index()
        return data_B0, data_B0bar

if __name__ == '__main__':

    dataPath = '/home/anna/master_thesis/data/whole_run2/'
    # Jpsi_file1 = '/home/anna/master_thesis/data/whole_run2/B2KstJpsi_2016_sWeight_wL0.root'
    # Jpsi_data1 = pandas.DataFrame(root_numpy.root2array(Jpsi_file1, treename='DecayTree'))
    # print(len(Jpsi_data1))
    # Jpsi_file2 = '/home/anna/master_thesis/data/whole_run2/B2KstJpsi_2017_sWeight_wL0.root'
    # Jpsi_data2 = pandas.DataFrame(root_numpy.root2array(Jpsi_file2, treename='DecayTree'))
    # print(len(Jpsi_data2))
    #
    # data = MergeDataframes([Jpsi_data1, Jpsi_data2])
    # print('Total length data:', len(data))
    # print(data.head(10))
    # L0trigg_selection(data)

    data_object = Dataset(dataPath, years=(2016, 2017, 2018))
    #data = data_object.get_data()
    print('Length data after trigger selection:', len(data_object.data))
    datasets = data_object.CP_divide()
    print(len(datasets[0]))
    plt.hist(datasets[0]['B0_ID'])
    plt.show()

    print(len(datasets[1]))
    plt.hist(datasets[1]['B0_ID'])
    plt.show()


