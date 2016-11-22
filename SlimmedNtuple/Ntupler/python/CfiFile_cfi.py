import FWCore.ParameterSet.Config as cms

TFileService = cms.Service("TFileService", fileName = cms.string("SlimmedNtuple.root") )

demo = cms.EDAnalyzer('Ntupler'
)
