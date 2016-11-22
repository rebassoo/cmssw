import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")



process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1001) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/tmp/rebassoo/miniAOD_PAT_1.root'
        #'file:/tmp/rebassoo/B2846F07-B046-E611-8CC0-0CC47A13CDA0.root'
        #'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90001/B2846F07-B046-E611-8CC0-0CC47A13CDA0.root'
        'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90001/E49F4520-B046-E611-B7DF-003048F5B2F0.root'
    ),
 skipEvents=cms.untracked.uint32(17000)
)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v4'  #

#process.demo = cms.EDAnalyzer('Ntupler'
#)

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.load("SlimmedNtuple.Ntupler.CfiFile_cfi") 

#process.p = cms.Path(process.demo*process.dump)
process.p = cms.Path(process.demo)
