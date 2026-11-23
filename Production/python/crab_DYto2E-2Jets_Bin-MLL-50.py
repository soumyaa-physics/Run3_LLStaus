
import CRABClient
from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'DYto2E-2Jets_Bin-MLL-50'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/data/dust/user/vashisht/CMSSW_17_0_0_pre3/src/Run3_LLStaus/Production/python/priv_exonanoMC_bg_2024.py'
config.JobType.maxJobRuntimeMin = 15*60
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 2

config.section_("Data")
config.Data.inputDataset = '/DYto2E-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24DRPremix-140X_mcRun3_2024_realistic_v26-v4/AODSIM'
config.Data.inputDBS = 'global'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
config.Data.totalUnits = 27212

config.Data.outLFNDirBase = '/store/user/svashish/LLSTauProduction_RUN3/DYto2E-2Jets_Bin-MLL-50'
config.Data.publication = True
config.Data.outputDatasetTag = 'NANOAOD'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
