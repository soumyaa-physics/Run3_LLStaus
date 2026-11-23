
import CRABClient
from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'Run3_2024C'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/data/dust/user/vashisht/CMSSW_17_0_0_pre3/src/Run3_LLStaus/Production/python/exoNanoData_PAT_NANO.py'
config.JobType.maxJobRuntimeMin = 15*60
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 4

config.section_("Data")
config.Data.inputDataset = '/Muon0/Run2024C-2024CDEReprocessing-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.lumiMask = '/data/dust/user/vashisht/CMSSW_17_0_0_pre3/src/Run3_LLStaus/Production/python/Cert_Collisions2024_378981_386951_Golden.json'

config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 11000    # ~480 min per job
config.Data.totalUnits = 1000000

config.Data.outLFNDirBase = '/store/user/svashish/LLSTauProduction_RUN3/Data_EXONanoAOD'
config.Data.publication = True
config.Data.outputDatasetTag = 'NANOAOD'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'

# to test: crab submit -c crab_sub_data.py --dryrun 