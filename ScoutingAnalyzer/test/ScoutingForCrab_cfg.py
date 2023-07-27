# Note to self: I smell a memory leak
import FWCore.ParameterSet.Config as cms


# Define the process
process = cms.Process("LL")


# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )



process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(),
                            #fileNames = cms.untracked.vstring(fname),
                        )

# Load the standard set of configuration modules
#process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')

##--- l1 stage2 digis ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')

# Load the global tag for the data in the command line arguments!
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "123X_dataRun3_HLT_v14"

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("output.root")
                               )


#The triggers used for the SampleMuTrigNew file
L1Info = ['L1_DoubleMu_12_5', 'L1_DoubleMu_15_7', 'L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7', 'L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18', 'L1_DoubleMu4_SQ_OS_dR_Max1p2', 'L1_DoubleMu4p5_SQ_OS_dR_Max1p2']

# Make tree
process.mmtree = cms.EDAnalyzer('ScoutingAnalyzer',
                                triggerresults = cms.InputTag("TriggerResults::HLT"),
                                primaryVtx = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx", "HLT"),
                                secondaryVtx = cms.InputTag("hltScoutingMuonPacker", "displacedVtx", "HLT"),
                                pfparticles = cms.InputTag("hltScoutingPFPacker", "" ,"HLT"),
                                electrons = cms.InputTag("hltScoutingEgammaPacker", "" ,"HLT"),
                                muons = cms.InputTag("hltScoutingMuonPacker", "" ,"HLT"),
                                rho = cms.InputTag("hltScoutingPFPacker", "rho", "HLT"),
                                AlgInputTag = cms.InputTag("gtStage2Digis"),
                                l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1Seeds = cms.vstring(L1Info),
                                doL1 = cms.bool(True)
                            )

#Following process is to get Stage 2 Digis for L1 info
process.p = cms.Path( process.gtStage2Digis*process.mmtree )
#process.p = cms.Path( process.mmtree )