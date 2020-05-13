import ROOT
import math
from ROOT import TLorentzVector
import pandas
import root_numpy

def calc_costhetal(vecA, vecB, mother, ID):

    vecA_copy = TLorentzVector(vecA.X(), vecA.Y(), vecA.Z(), vecA.T())
    vecB_copy = TLorentzVector(vecB.X(), vecB.Y(), vecB.Z(), vecB.T())
    mother_copy = TLorentzVector(mother.X(), mother.Y(), mother.Z(), mother.T())
    vecAB = vecA+vecB
    # Speed of vA+Vb ref. frame in the lab frame
    v_boost = vecAB.BoostVector() # v_boost = p_x/E, p_y/E, p_z/E, with p momentum of vA+vB
    vecA_copy.Boost(-v_boost) # Boost from lab to vA+vB reference frame
    vecB_copy.Boost(-v_boost)
    mother_copy.Boost(-v_boost)
    vec3_A = vecA_copy.Vect()
    vec3_B = vecB_copy.Vect()
    vec3_mother = mother_copy.Vect()

    if ID > 0:
        ct = math.cos(vec3_A.Angle(-vec3_mother))
    else:
        ct = math.cos(vec3_B.Angle(-vec3_mother))

    return ct

def calc_costhetak(vecA, vecB, mother):

    vecA_copy = TLorentzVector(vecA.X(), vecA.Y(), vecA.Z(), vecA.T())
    vecB_copy = TLorentzVector(vecB.X(), vecB.Y(), vecB.Z(), vecB.T())
    mother_copy = TLorentzVector(mother.X(), mother.Y(), mother.Z(), mother.T())
    vecAB = vecA + vecB
    # Speed of vA+Vb ref. frame in the lab frame
    v_boost = vecAB.BoostVector() # v_boost = p_x/E, p_y/E, p_z/E, with p momentum of vA+vB
    vecA_copy.Boost(-v_boost) # Boost from lab to vA+vB reference frame
    vecB_copy.Boost(-v_boost)
    mother_copy.Boost(-v_boost)
    vec3_A = vecA_copy.Vect()
    vec3_B = vecB_copy.Vect()
    vec3_mother = mother_copy.Vect()
    ct = math.cos(vec3_A.Angle(-vec3_mother))
    return ct

def calc_phi(vecA, vecB, vecC, vecD, mother, ID):
    vecA_copy = TLorentzVector(vecA.X(), vecA.Y(), vecA.Z(), vecA.T())
    vecB_copy = TLorentzVector(vecB.X(), vecB.Y(), vecB.Z(), vecB.T())
    vecC_copy = TLorentzVector(vecC.X(), vecC.Y(), vecC.Z(), vecC.T())
    vecD_copy = TLorentzVector(vecD.X(), vecD.Y(), vecD.Z(), vecD.T())
    mother_copy = TLorentzVector(mother.X(), mother.Y(), mother.Z(), mother.T())
    vecAB = vecA + vecB
    vecCD = vecC + vecD
    v_boost = mother_copy.BoostVector() # Speed of mother rest frame in the lab frame
    vecA_copy.Boost(-v_boost)
    vecB_copy.Boost(-v_boost)
    vecC_copy.Boost(-v_boost)
    vecD_copy.Boost(-v_boost)
    vecCD.Boost(-v_boost)
    vec3CD = vecCD.Vect()
    normalAB = vecA_copy.Vect().Cross(vecB_copy.Vect()) # Cross gives the vector product
    normalCD = vecC_copy.Vect().Cross(vecD_copy.Vect())

    if ID>0:
        phi = normalCD.Angle(normalAB)
        if((normalAB.Cross(normalCD)).Dot(vec3CD) < 0.0):
            phi = -phi
    else:
        phi = normalCD.Angle(-normalAB)
        if ((normalAB.Cross(normalCD)).Dot(vec3CD) < 0.0):
            phi = -phi
    return phi

def calc_angvar(mplus, mminus, kaon, pion, B, ID):
    costhetal = calc_costhetal(mplus, mminus, B, ID)
    costhetak = calc_costhetak(kaon, pion, B)
    phi = calc_phi(mplus, mminus, kaon, pion, B, ID)
    return costhetal, costhetak, phi

def add_angvar_todata(data):
    for index, row in data.iterrows():
        # if index%1000==0:
        #     print('Processing event', index)
        mu_plus = TLorentzVector(row['mu_plus_PX'], row['mu_plus_PY'], row['mu_plus_PZ'], row['mu_plus_PE'])
        mu_minus = TLorentzVector(row['mu_minus_PX'], row['mu_minus_PY'], row['mu_minus_PZ'], row['mu_minus_PE'])
        K = TLorentzVector(row['K_PX'], row['K_PY'], row['K_PZ'], row['K_PE'])
        Pi = TLorentzVector(row['Pi_PX'], row['Pi_PY'], row['Pi_PZ'], row['Pi_PE'])
        B = mu_plus + mu_minus + K + Pi
        ID = row['B0_ID']
        data['costhetal'].values[index], data['costhetak'].values[index], data['phi'].values[index] = calc_angvar(mu_plus, mu_minus, K, Pi, B, ID)
        #print(row['phi'])
    # print(data['phi'])

if __name__ == '__main__':

    # sigFile = '/home/anna/master_thesis/forBDT/B2KstJpsi_2016_sWeight_wL0.root'
    bkgFile = '/home/anna/master_thesis/forBDT/B2Kstmumu_sideband_2016.root'
    #
    # sig = pandas.DataFrame(root_numpy.root2array(sigFile, treename='DecayTree'))
    bkg = pandas.DataFrame(root_numpy.root2array(bkgFile, treename='DecayTree'))
    #
    bkg['phi'] = 0.
    print(bkg['phi'])
    add_angvar_todata(bkg)
    print(bkg['phi'])