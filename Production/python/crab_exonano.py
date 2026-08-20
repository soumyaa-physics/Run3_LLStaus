# TEMPLATE used for automatic script submission of multiple datasets

# from CRABClient.UserUtilities import getUsernameFromCRIC

# username = svashish
import CRABClient
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'SUS-RunIIISummer2024_EXONANO-stau700_lsp1_ctau10m_v1'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/data/dust/user/vashisht/CMSSW_17_0_0_pre3/src/Run3_LLStaus/Production/python/exoNanoMC_2024_priv.py'

config.JobType.maxJobRuntimeMin = 10*60
config.JobType.maxMemoryMB = 3000
config.JobType.numCores = 4

config.section_("Data")

config.Data.inputDataset = '/SUS-Run3Summer2024-stau700_lsp1_ctau10m/svashish-RECO-77a50869adebf22fb98f8bb4cf16c171/USER'
config.Data.inputDBS = 'phys03'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/user/svashish/LLSTauProduction_RUN3/SUS-RunIIISummer2024-stau700_lsp1_ctau10m'
config.Data.publication = True
config.Data.outputDatasetTag = 'NANOAOD'

config.section_("Site")
config.Site.whitelist = ['T2_DE_DESY']
config.Site.storageSite = 'T2_DE_DESY'
