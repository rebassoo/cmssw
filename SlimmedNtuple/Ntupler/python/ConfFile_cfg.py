import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1001) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/tmp/rebassoo/miniAOD_PAT_1.root'
        #'file:/tmp/rebassoo/B2846F07-B046-E611-8CC0-0CC47A13CDA0.root'
        #'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90001/B2846F07-B046-E611-8CC0-0CC47A13CDA0.root'
        #'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90001/E49F4520-B046-E611-B7DF-003048F5B2F0.root'
        'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90000/1C4F8054-D845-E611-8251-FA163E149DE8.root'
        #'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90002/8A922D2A-5C47-E611-B3E7-001A648F1BFA.root'
        #'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90000/A656B0C9-BF45-E611-9921-002590200A28.root'
        #'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/AOD/01Jul2016-v1/90000/CCB0A6F2-B145-E611-A8AB-B083FEC76567.root'
        #'file:/hadoop/cms/phedex/store/data/Run2016B/DoubleMuon/MINIAOD/01Jul2016-v1/90000/FEF32273-1647-E611-9AA3-848F69FD2D6F.root'
        #'root://cmsxrootd.fnal.gov//store/data/Run2016F/DoubleMuon/AOD/23Sep2016-v1/50000/9E0F8266-0590-E611-A14E-0242AC130002.root'
    )#,
   #skipEvents=cms.untracked.uint32(3420)
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


#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
#useAOD = True
#if useAOD == True :
#    dataFormat = DataFormat.AOD
#else :
#    dataFormat = DataFormat.MiniAOD

#switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
#my_id_modules = [
#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
#    ]

#add them to the VID producer
#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



#process.p = cms.Path(process.demo*process.dump)
#process.p = cms.Path(process.egmGsfElectronIDSequence * process.demo)
process.p = cms.Path(process.demo)
