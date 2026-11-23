
import CRABClient
from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'TTto4Q'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/data/dust/user/vashisht/CMSSW_17_0_0_pre3/src/Run3_LLStaus/Production/python/exonanoMC_bg_2024.py'
config.JobType.maxJobRuntimeMin = 10*60
config.JobType.maxMemoryMB = 3000
config.JobType.numCores = 4

config.section_("Data")
config.Data.inputDataset = '/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24DRPremix-140X_mcRun3_2024_realistic_v26-v2/AODSIM'
config.Data.inputDBS = 'global'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/user/svashish/LLSTauProduction_RUN3/TTto4Q'
config.Data.publication = True
config.Data.outputDatasetTag = 'NANOAOD'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
