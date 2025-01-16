import FWCore.ParameterSet.Config as cms

# Define the process
process = cms.Process("MVA")

# Load standard sequences and RecoMTD modules
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# Define the input source (ROOT file with track data)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_14_2_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v1_STD_2026D110_PU-v1/2580000/09052bdb-baff-419a-92d6-08f7f4009336.root')  # Replace with your input file
)

#Set the number of events to process (optional, adjust as needed)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(3))

from RecoMTD.TimingIDTools.mtdTrackQualityMVAProducer_cfi import *

process.mtdTrackQualityMVA = mtdTrackQualityMVAProducer.clone()


# Path to run the module
process.p = cms.Path(process.mtdTrackQualityMVA)

# Optional: Set up output (you can adjust as needed)
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output_test.root')
    #fileName = cms.untracked.string('/eos/user/n/nstrautn/MTD_BDT/output_test.root')
)
process.outpath = cms.EndPath(process.output)


def customise(process):
    # Adding SimpleMemoryCheck service
    process.SimpleMemoryCheck = cms.Service(
        "SimpleMemoryCheck",
        ignoreTotal=cms.untracked.int32(1),
        oncePerEventMode=cms.untracked.bool(False)
    )
    
    # Adding Timing service
    process.Timing = cms.Service(
        "Timing",
        summaryOnly=cms.untracked.bool(True),
        excessiveTimeThreshold=cms.untracked.double(600)
    )
    
    # Adding timing summary to options
    if hasattr(process, "options"):
        process.options.wantSummary = cms.untracked.bool(True)
    else:
        process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True)
        )
    
    return process

def customiseWithTimeMemorySummary(process):
    return customise(process)

# Apply customisations
process = customiseWithTimeMemorySummary(process)