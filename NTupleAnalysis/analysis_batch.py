#%%
import time
import numpy as np
from numba import njit
import boost_histogram as bh
import awkward as ak 
import uproot 
import matplotlib.pyplot as plt
import mplhep
import hist
mplhep.style.use("CMS")
plt.rcParams["figure.figsize"] = (12.5,10)

import vector
vector.register_awkward()

# NanoEvents
from coffea.nanoevents import NanoEventsFactory, BaseSchema 
from coffea.nanoevents import NanoAODSchema
# Processors
import coffea.processor as processor

num_files_bparking18 = 100
num_files_bparking23 = 10
num_files_scouting = 100

n_cores = 8

class MyProcessorBParking18(processor.ProcessorABC):
    def __init__(self):
        
        jpsi_axis = hist.axis.Regular(name="jpsi", label="SV mass(GeV)", bins=50, start=2.5, stop=3.5)
        lowmass_axis = hist.axis.Regular(name="lowmass", label="SV mass(GeV)", bins=101, start=0, stop=22)
        mass_axis = hist.axis.Regular(name="mass", label="Dimuon mass(GeV)", bins=101, start=0, stop=100)
        mudxysig_axis = hist.axis.Regular(name="mudxysig", label="Muon dxysig", bins=101, start=0, stop=100)
        svdxysig_axis = hist.axis.Regular(name="svdxysig", label="SV dxysig", bins=101, start=0, stop=100)
        mudxy_axis = hist.axis.Regular(name="mudxy", label="Muon dxy (cm)", bins=101, start=0, stop=100)
        svdxy_axis = hist.axis.Regular(name="svdxy", label="SV dxy (cm)", bins=101, start=0, stop=100)

        self.output = processor.dict_accumulator({
            'num_events': 0,
            'jpsi_hist': hist.Hist(jpsi_axis),
            'lowmass_hist': hist.Hist(lowmass_axis),
            'mass_hist': hist.Hist(mass_axis),
            'mudxysig_hist': hist.Hist(mudxysig_axis),
            'svdxysig_hist': hist.Hist(svdxysig_axis),
            'mudxy_hist': hist.Hist(mudxy_axis),
            'svdxy_hist': hist.Hist(svdxy_axis)
        })
        
    def process(self, events):
        
        MuonSV = events.muonSV
        Muons = events.Muon

        
        
        #Cut on Muons
        mu_cut = (np.abs(Muons.dxy/Muons.dxyErr) > 5)
        Muons = Muons[mu_cut]

        #Sort SVs to get only min chi2/ndof one
        min_args = ak.argsort(np.abs(MuonSV['chi2']/MuonSV['ndof']), ascending=True)
        MuonSV = MuonSV[min_args]

        mu_pairs = ak.combinations(Muons, 2)
        mu1, mu2 = ak.unzip(mu_pairs)
        mu_mass = np.sqrt(2*mu1.pt*mu2.pt*(np.cosh(mu1.eta - mu2.eta) - np.cos(mu1.phi - mu2.phi)))
        #dxy_cut = (np.abs(mu1.dxy) > 5) & (np.abs(mu2.dxy) > 5)

        #Cut on Muon SV's
        sv_cut = MuonSV["dxySig"] > 5
        MuonSV = MuonSV[sv_cut]
        
        
        #Get the first sorted MuonSV (remove nones) => MuonSV is now flattened
        MuonSV = ak.firsts(MuonSV)
        MuonSV = MuonSV[~ak.is_none(MuonSV)]
        

        self.output['num_events'] += len(events)
        self.output['jpsi_hist'].fill(MuonSV.mass)
        self.output['lowmass_hist'].fill(MuonSV.mass)
        self.output['mass_hist'].fill(ak.flatten(mu_mass))
        self.output['mudxysig_hist'].fill(ak.flatten(np.abs(Muons.dxy/Muons.dxyErr)))
        self.output['svdxysig_hist'].fill(MuonSV['dxySig'])
        self.output['mudxy_hist'].fill(ak.flatten(np.abs(Muons.dxy)))
        self.output['svdxy_hist'].fill(MuonSV['dxy'])

        return self.output
    
    def postprocess(self, accumulator):
        pass

class MyProcessorBParking23(processor.ProcessorABC):
    def __init__(self):

        jpsi_axis = hist.axis.Regular(name="jpsi", label="SV mass(GeV)", bins=50, start=2.5, stop=3.5)
        #jpsi_disp_axis = hist.axis.Regular(name="jpsi_disp", label="SV mass(GeV)", bins=50, start=2.5, stop=3.5)
        lowmass_axis = hist.axis.Regular(name="lowmass", label="SV mass(GeV)", bins=101, start=0, stop=22)
        mass_axis = hist.axis.Regular(name="mass", label="Dimuon mass(GeV)", bins=101, start=0, stop=100)
        mudxysig_axis = hist.axis.Regular(name="mudxysig", label="Muon dxysig", bins=101, start=0, stop=100)
        svdxysig_axis = hist.axis.Regular(name="svdxysig", label="SV dxysig", bins=101, start=0, stop=100)
        mudxy_axis = hist.axis.Regular(name="mudxy", label="Muon dxy (cm)", bins=101, start=0, stop=100)
        svdxy_axis = hist.axis.Regular(name="svdxy", label="SV dxy (cm)", bins=101, start=0, stop=100)

        self.output = processor.dict_accumulator({
            'num_events': 0,
            'jpsi_hist': hist.Hist(jpsi_axis),
            'jpsi_disp_hist': hist.Hist(jpsi_axis),
            'lowmass_hist': hist.Hist(lowmass_axis),
            'mass_hist': hist.Hist(mass_axis),
            'mudxysig_hist': hist.Hist(mudxysig_axis),
            'svdxysig_hist': hist.Hist(svdxysig_axis),
            'mudxy_hist': hist.Hist(mudxy_axis),
            'svdxy_hist': hist.Hist(svdxy_axis)
        })

    def process(self, events):
        HLTDecision = events.HLT
        MuonSV = events.SV
        Muons = events.Muon

        #Cut on displaced trigger
        DispCut = HLTDecision['DoubleMu4_LowMass_Displaced'] == True

        #Cut on Muons
        mu_cut = (np.abs(Muons.dxy/Muons.dxyErr) > 5)
        Muons = Muons[mu_cut]

        #Sort SVs to get only min chi2/ndof one
        min_args = ak.argsort(np.abs(MuonSV['chi2']/MuonSV['ndof']), ascending=True)
        MuonSV = MuonSV[min_args]

        mu_pairs = ak.combinations(Muons, 2)
        mu1, mu2 = ak.unzip(mu_pairs)
        mu_mass = np.sqrt(2*mu1.pt*mu2.pt*(np.cosh(mu1.eta - mu2.eta) - np.cos(mu1.phi - mu2.phi)))
        #dxy_cut = (np.abs(mu1.dxy) > 5) & (np.abs(mu2.dxy) > 5)

        #Cut on Muon SV's
        sv_cut = MuonSV["dxySig"] > 5
        MuonSV = MuonSV[sv_cut]

        #Get trigger cut SVs
        MuonSVDisp = MuonSV[DispCut]

        #Get the first sorted MuonSV (remove nones)
        MuonSV = ak.firsts(MuonSV)
        MuonSV = MuonSV[~ak.is_none(MuonSV)]
        MuonSVDisp = ak.firsts(MuonSVDisp)
        MuonSVDisp = MuonSVDisp[~ak.is_none(MuonSVDisp)]

        self.output['num_events'] += len(events)
        self.output['jpsi_hist'].fill(MuonSV.mass)
        self.output['jpsi_disp_hist'].fill(MuonSVDisp.mass)
        self.output['lowmass_hist'].fill(MuonSV.mass)
        self.output['mass_hist'].fill(ak.flatten(mu_mass))
        self.output['mudxysig_hist'].fill(ak.flatten(np.abs(Muons.dxy/Muons.dxyErr)))
        self.output['svdxysig_hist'].fill(MuonSV['dxySig'])
        self.output['mudxy_hist'].fill(ak.flatten(np.abs(Muons.dxy)))
        self.output['svdxy_hist'].fill(MuonSV['dxy'])

        return self.output
    
    def postprocess(self, accumulator):
        pass

class MyProcessorScouting(processor.ProcessorABC):
    def __init__(self):
        jpsi_axis = hist.axis.Regular(name="jpsi", label="SV mass(GeV)", bins=50, start=2.5, stop=3.5)
        lowmass_axis = hist.axis.Regular(name="lowmass", label="SV mass(GeV)", bins=101, start=0, stop=22)
        mass_axis = hist.axis.Regular(name="mass", label="Dimuon mass(GeV)", bins=101, start=0, stop=100)
        mudxysig_axis = hist.axis.Regular(name="mudxysig", label="Muon dxysig", bins=101, start=0, stop=100)
        svdxysig_axis = hist.axis.Regular(name="svdxysig", label="SV dxysig", bins=101, start=0, stop=100)
        mudxy_axis = hist.axis.Regular(name="mudxy", label="Muon dxy (cm)", bins=101, start=0, stop=100)
        svdxy_axis = hist.axis.Regular(name="svdxy", label="SV dxy (cm)", bins=101, start=0, stop=100)

        self.output = processor.dict_accumulator({
            'num_events': 0,
            'jpsi_hist': hist.Hist(jpsi_axis),
            'lowmass_hist': hist.Hist(lowmass_axis),
            'mass_hist': hist.Hist(mass_axis),
            'mudxysig_hist': hist.Hist(mudxysig_axis),
            'svdxysig_hist': hist.Hist(svdxysig_axis),
            'mudxy_hist': hist.Hist(mudxy_axis),
            'svdxy_hist': hist.Hist(svdxy_axis)
        })

    def process(self, events):

        EventPV = events.pVtx
        EventPV = ak.firsts(EventPV)
        MuonSV = events.sVtx
        Muons = events.Muon

        #Cut on Muons
        mu_cut = (np.abs(Muons.dxy/Muons.dxyerror) > 5)
        Muons = Muons[mu_cut]

        #Sort SVs to get only min chi2/ndof one
        min_args = ak.argsort(np.abs(MuonSV['chi2']/MuonSV['ndof']), ascending=True)
        MuonSV = MuonSV[min_args]

        mu_pairs = ak.combinations(Muons, 2)
        mu1, mu2 = ak.unzip(mu_pairs)
        mu_mass = np.sqrt(2*mu1.pt*mu2.pt*(np.cosh(mu1.eta - mu2.eta) - np.cos(mu1.phi - mu2.phi)))
        
        #Modification where no extra calcs needed
        #Cut on Muon SV's
        sv_cut = MuonSV["dxySig"] > 5
        MuonSV = MuonSV[sv_cut]
        
        #Get the first sorted MuonSV (remove nones)
        MuonSV = ak.firsts(MuonSV)
        MuonSV = MuonSV[~ak.is_none(MuonSV)]
        """
        #Calculation to get muon sv dxy
        sv_mult_cut = ak.num(MuonSV) > 0
        EventPV_filtered = EventPV[sv_mult_cut]

        #Returns the muon sv's we care about
        MuonSV_filtered = MuonSV[sv_mult_cut]

        #Calculate dxy
        dx = MuonSV_filtered['x'] - (ak.ones_like(MuonSV_filtered['x'])*EventPV_filtered['x']) 
        dy = MuonSV_filtered['y'] - (ak.ones_like(MuonSV_filtered['y'])*EventPV_filtered['y']) 

        dxerr = np.sqrt((MuonSV_filtered['xError'])**2 + (EventPV_filtered['xError'])**2)
        dyerr = np.sqrt((MuonSV_filtered['yError'])**2 + (EventPV_filtered['yError'])**2)

        dxy = np.sqrt(dx**2 + dy**2)
        dxyerr = np.sqrt((dx*dxerr)**2 + (dy*dyerr)**2)/dxy
        dxysig = dxy/dxyerr
        

        #Apply secondary vertex dxysig cut
        sv_cut = dxysig > 5
        dxy = dxy[sv_cut]
        dxysig = dxysig[sv_cut]
        MuonSV_filtered = MuonSV_filtered[sv_cut]

        #Pick out the first sorted SV from already cut values (remove nones)
        dxy = ak.firsts(dxy)
        dxysig = ak.firsts(dxysig)
        dxy = dxy[~ak.is_none(dxy)]
        dxysig = dxysig[~ak.is_none(dxysig)]

        #Pick out the first sorted SV from already cut values (remove nones)
        MuonSV = ak.firsts(MuonSV_filtered)
        MuonSV = MuonSV[~ak.is_none(MuonSV)]
        """
        """
        Quantities used in plotting
        dxy, dxysig, MuonSV_filtered and Muons

        Cuts applied:
            Muons: dxysig > 5
            MuonSV_filtered: dxysig > 5, first min chi2/ndof
        
        MuonSV_filtered is the unflattened array of muon sv's where zero sv
        events are removed
        """        

        self.output['num_events'] += len(events)
        self.output['jpsi_hist'].fill(MuonSV.mass)
        self.output['lowmass_hist'].fill(MuonSV.mass)
        self.output['mass_hist'].fill(ak.flatten(mu_mass))
        self.output['mudxysig_hist'].fill(ak.flatten(np.abs(Muons.dxy/Muons.dxyerror)))
        #self.output['svdxysig_hist'].fill(dxysig)
        self.output['svdxysig_hist'].fill(MuonSV['dxySig'])
        self.output['mudxy_hist'].fill(ak.flatten(np.abs(Muons.dxy)))
        #self.output['svdxy_hist'].fill(dxy)
        self.output['svdxy_hist'].fill(MuonSV['dxy'])


        return self.output
    
    def postprocess(self, accumulator):
        pass



#Process b parking 2018
start_bparking18 = time.time()

#Make the list
#bparking_list_file = open("fileset_bparking.txt", "r")
#bparking_list = bparking_list_file.read()
#bparking_list = bparking_list.split("\n")

#Manual method
bparking18_list = []

with open('fileset_bparking18_proc.txt', 'w') as file:
    for i in range(num_files_bparking18):
        fname = "/vols/cms/mc3909/bparkProductionAll_V1p0/ParkingBPH1_Run2018B-05May2019-v2_MINIAOD_v1p0_generationSync/output_{:d}.root".format((i+1))
        bparking18_list.append(fname)
        file.write(fname+"\n")


fileset_bparking18 = {"b-parking" : bparking18_list}


futures_run_bparking18 = processor.Runner(
    executor = processor.FuturesExecutor(compression=None, workers=n_cores),
    schema=NanoAODSchema
)

out_bparking18 = futures_run_bparking18(
    fileset_bparking18,
    treename='Events',
    processor_instance= MyProcessorBParking18()
)

end_bparking18 = time.time()
time_bparking18 = (end_bparking18 - start_bparking18)

print("\nB parking 2018 processed")

#Process b parking 2023
start_bparking23 = time.time()

bparking23_list = []

with open('fileset_bparking23_proc.txt', 'w') as file:
    for i in range(num_files_bparking23):
        fname = "/vols/cms/pb4918/StoreNTuple/BParking23/output_{:d}.root".format((i+1))
        bparking23_list.append(fname)
        file.write(fname+"\n")

fileset_bparking23 = {"b-parking" : bparking23_list}


futures_run_bparking23 = processor.Runner(
    executor = processor.FuturesExecutor(compression=None, workers=n_cores),
    schema=NanoAODSchema
)

out_bparking23 = futures_run_bparking23(
    fileset_bparking23,
    treename='Events',
    processor_instance= MyProcessorBParking23()
)

end_bparking23 = time.time()
time_bparking23 = (end_bparking23 - start_bparking23)

print("\nB parking 2023 processed")


#Process scouting
start_scouting = time.time()

#scouting_list_file = open("fileset_scouting.txt", "r")
#scouting_list = scouting_list_file.read()
#scouting_list = scouting_list.split("\n")

scouting_list = []

#List of failed CRAB jobs
skip_files = [17, 24, 40, 43]

with open('fileset_scouting_proc.txt', 'w') as file:
    for i in range(num_files_scouting):
        if (i+1) in skip_files: continue

        else:
            fname = "/vols/cms/pb4918/StoreNTuple/Scouting/2022FStoreExpanded/output_{:d}.root".format((i+1))
            scouting_list.append(fname)
            file.write(fname+"\n")

fileset_scouting = {"scouting" : scouting_list}



futures_run_scouting = processor.Runner(
    executor = processor.FuturesExecutor(compression=None, workers=n_cores),
    schema=NanoAODSchema
)

out_scouting = futures_run_scouting(
    fileset_scouting,
    treename='mmtree/tree',
    processor_instance= MyProcessorScouting()
)

end_scouting = time.time()
time_scouting = (end_scouting - start_scouting)

print("\nScouting processed")

event_num_bparking18 = out_bparking18['num_events']
event_num_bparking23 = out_bparking23['num_events']
event_num_scouting = out_scouting['num_events']


print("\nSummary:")
print("Number of b-parking events 2018 processed:", event_num_bparking18)
print("Processed in %.3f seconds"%(time_bparking18))
print("Number of b-parking events 2023 processed:", event_num_bparking23)
print("Processed in %.3f seconds"%(time_bparking23))
print("\nNumber of scouting events processed:", event_num_scouting)
print("Processed in %.3f seconds"%(time_scouting))
#%%
#Fitting
dxy_bins = np.linspace(0,50,101)
dxysig_bins = np.linspace(0,100,101)
jpsi_bins = np.linspace(2.5,3.5,50)
low_bins = np.linspace(0,22,100)
mass_bins = np.linspace(0,100,102)
#%%
#Plotting stuff 
fig, ax = plt.subplots()

jpsi_hist_bparking18 = out_bparking18['jpsi_hist']
jpsi_hist_bparking23 = out_bparking23['jpsi_hist']
jpsi_disp_hist_bparking23 = out_bparking23['jpsi_disp_hist']
jpsi_hist_scouting = out_scouting['jpsi_hist']


#Scale up to 1/fb
jpsi_hist_bparking18 = jpsi_hist_bparking18/(0.070226436/6)
jpsi_hist_bparking23 = jpsi_hist_bparking23/(0.496789061/8)
jpsi_disp_hist_bparking23 = jpsi_disp_hist_bparking23/(0.496789061/8)
jpsi_hist_scouting = jpsi_hist_scouting/0.034031063

mplhep.histplot(jpsi_hist_bparking18, color='blue', label='B-parking 2018', density=False)
mplhep.histplot(jpsi_hist_bparking23, color='green', label='B-parking 2023', density=False)
mplhep.histplot(jpsi_disp_hist_bparking23, color='purple', label='B-parking 2023 (Disp trigger)', density=False)
mplhep.histplot(jpsi_hist_scouting, color='red', label='Scouting', density=False)
ax.set_xlabel('Dimuon invariant mass for SV (GeV)')
ax.set_ylabel('Secondary vertices')
ax.text(0.70,0.70, "SV dxysig > 5", color='black', fontsize=18, ha='left', transform=ax.transAxes)
ax.text(0.70,0.65, "min $\\chi^{2}$/ndof", color='black', fontsize=18, ha='left', transform=ax.transAxes)
ax.legend(loc='upper left')
ax.set_yscale('log')
mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Scaled to 1 $\mathrm{fb}^{-1}$')
#mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Normalised')
fig.savefig("plots/sv_jpsi_comparison.png")


fig, ax = plt.subplots()

lowmass_hist_bparking18 = out_bparking18['lowmass_hist']
lowmass_hist_bparking23 = out_bparking23['lowmass_hist']
lowmass_hist_scouting = out_scouting['lowmass_hist']

#Scale up to 1/fb
lowmass_hist_bparking18 = lowmass_hist_bparking18/(0.070226436/6)
lowmass_hist_bparking23 = lowmass_hist_bparking23/(0.496789061/8)
lowmass_hist_scouting = lowmass_hist_scouting/0.034031063

mplhep.histplot(lowmass_hist_bparking18, color='blue', label='B-parking 2018', density=False)
mplhep.histplot(lowmass_hist_bparking23, color='green', label='B-parking 2023', density=False)
mplhep.histplot(lowmass_hist_scouting, color='red', label='Scouting', density=False)
ax.set_xlabel('Dimuon invariant mass for SV (GeV)')
ax.set_ylabel('Secondary vertices')
ax.text(0.70,0.75, "SV dxysig > 5", color='black', fontsize=18, ha='left', transform=ax.transAxes)
ax.text(0.70,0.70, "min $\\chi^{2}$/ndof", color='black', fontsize=18, ha='left', transform=ax.transAxes)
ax.legend(loc='upper right')
ax.set_yscale('log')
mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Scaled to 1 $\mathrm{fb}^{-1}$')
#mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Normalised')
fig.savefig("plots/sv_low_comparison.png")


fig, ax = plt.subplots()

mass_hist_bparking18 = out_bparking18['mass_hist']
mass_hist_bparking23 = out_bparking23['mass_hist']
mass_hist_scouting = out_scouting['mass_hist']

#Scale up to 1/fb
mass_hist_bparking18 = mass_hist_bparking18/(0.070226436/6)
mass_hist_bparking23 = mass_hist_bparking23/(0.496789061/8)
mass_hist_scouting = mass_hist_scouting/0.034031063

mplhep.histplot(mass_hist_bparking18, color='blue', label='B-parking 2018', density=False)
mplhep.histplot(mass_hist_bparking23, color='green', label='B-parking 2023', density=False)
mplhep.histplot(mass_hist_scouting, color='red', label='Scouting', density=False)
ax.set_xlabel('Dimuon invariant mass (GeV)')
ax.set_ylabel('Muon pairs')
ax.text(0.70,0.75, "Mu dxysig > 5", color='black', fontsize=18, ha='left', transform=ax.transAxes)
ax.legend(loc='upper right')
ax.set_yscale('log')
mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Scaled to 1 $\mathrm{fb}^{-1}$')
#mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Normalised')
fig.savefig("plots/dimu_mass_comparison.png")


fig, ax = plt.subplots()

svdxy_hist_bparking18 = out_bparking18['svdxy_hist']
svdxy_hist_bparking23 = out_bparking23['svdxy_hist']
svdxy_hist_scouting = out_scouting['svdxy_hist']

svdxy_hist_bparking18 = svdxy_hist_bparking18/(0.070226436/6)
svdxy_hist_bparking23 = svdxy_hist_bparking23/(0.496789061/8)
svdxy_hist_scouting = svdxy_hist_scouting/0.034031063

mplhep.histplot(svdxy_hist_bparking18, color='blue', label='B-parking 2018', density=False)
mplhep.histplot(svdxy_hist_bparking23, color='green', label='B-parking 2023', density=False)
mplhep.histplot(svdxy_hist_scouting, color='red', label='Scouting', density=False)
ax.set_xlabel('SV dxy (cm)')
ax.set_ylabel('Secondary vertices')
ax.text(0.70,0.70, "SV dxysig > 5", color='black', fontsize=18, ha='left', transform=ax.transAxes)
ax.text(0.70,0.65, "min $\\chi^{2}$/ndof", color='black', fontsize=18, ha='left', transform=ax.transAxes)
ax.legend(loc='upper right')
ax.set_yscale('log')
mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Scaled to 1 $\mathrm{fb}^{-1}$')
#mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Normalised')
fig.savefig("plots/sv_dxy_comparison.png")


fig, ax = plt.subplots()

svdxysig_hist_bparking18 = out_bparking18['svdxysig_hist']
svdxysig_hist_bparking23 = out_bparking23['svdxysig_hist']
svdxysig_hist_scouting = out_scouting['svdxysig_hist']

svdxysig_hist_bparking18 = svdxysig_hist_bparking18/(0.070226436/6)
svdxysig_hist_bparking23 = svdxysig_hist_bparking23/(0.496789061/8)
svdxysig_hist_scouting = svdxysig_hist_scouting/0.034031063

mplhep.histplot(svdxysig_hist_bparking18, color='blue', label='B-parking 2018', density=False)
mplhep.histplot(svdxysig_hist_bparking23, color='green', label='B-parking 2023', density=False)
mplhep.histplot(svdxysig_hist_scouting, color='red', label='Scouting', density=False)
ax.set_xlabel('SV dxysig')
ax.set_ylabel('Secondary vertices')
ax.text(0.70,0.70, "SV dxysig > 5", color='black', fontsize=18, ha='left', transform=ax.transAxes)
ax.text(0.70,0.65, "min $\\chi^{2}$/ndof", color='black', fontsize=18, ha='left', transform=ax.transAxes)
ax.legend(loc='upper right')
ax.set_yscale('log')
mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Scaled to 1 $\mathrm{fb}^{-1}$')
#mplhep.cms.label(data=True, label='Work In Progress', rlabel=r'Normalised')
fig.savefig("plots/sv_dxysig_comparison.png")

#%%
#Adding b-parking for just displaced trigger


#%%
#Define fitting functions for j-psi peak (taken from Mikael's icefit)
import iminuit
from iminuit import minimize
from iminuit.cost import LeastSquares
from scipy import interpolate
from scipy.optimize import curve_fit


def peak(x, mu, sigma, n, a, b):
    
    """
    Of form:
    y = Gauss + ax + b
    """
    
    gauss = (n / (np.abs(sigma) * np.sqrt(2*np.pi))) * np.exp(- 0.5 * ((x - mu)/sigma)**2)
    y = gauss + a*x + b

    return y

jpsi_hist_bparking = out_bparking18['jpsi_hist']
jpsi_hist_scouting = out_scouting['jpsi_hist']

#Scale up to 1/fb
mass_points = 0.5*(jpsi_bins[:-1] + jpsi_bins[1:])

jpsi_hist_bparking = jpsi_hist_bparking/((0.070226436/6))
jpsi_hist_scouting = jpsi_hist_scouting/0.001496566

initial_guess_bparking = [3.1, 0.2, 1e4, 0, 1e3]
initial_guess_scouting = [3.1, 0.2, 1e5, 0, 1e4]
popt_scouting, pcov_scouting = curve_fit(peak, mass_points, jpsi_hist_scouting, p0=initial_guess_scouting)
popt_bparking, pcov_bparking = curve_fit(peak, mass_points, jpsi_hist_bparking, p0=initial_guess_bparking)

#%%
#Summary info
print("\nSummary info: B-Parking")
print("J/Psi mass: {:.3f} +/- {:.3f} GeV".format(popt_bparking[0], popt_bparking[1]))
print("Yield: {:.2f} events".format(popt_bparking[2]))

print("\nSummary info: Scouting")
print("J/Psi mass: {:.3f} +/- {:.3f} GeV".format(popt_scouting[0], popt_scouting[1]))
print("Yield: {:.2f} events".format(popt_scouting[2]))

# %%
