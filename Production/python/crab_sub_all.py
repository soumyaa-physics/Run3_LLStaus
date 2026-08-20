import subprocess

samples_file = "stau_exonano_samples.txt"

template = """
import CRABClient
from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = '{request_name}'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/data/dust/user/vashisht/CMSSW_17_0_0_pre3/src/Run3_LLStaus/Production/python/exoNanoMC_2024_priv.py'
config.JobType.maxJobRuntimeMin = 10*60
config.JobType.maxMemoryMB = 3000
config.JobType.numCores = 4

config.section_("Data")
config.Data.inputDataset = '{input_dataset}'
config.Data.inputDBS = 'phys03'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/user/svashish/LLSTauProduction_RUN3/{request_name}'
config.Data.publication = True
config.Data.outputDatasetTag = 'NANOAOD'

config.section_("Site")
config.Site.whitelist = ['T2_DE_DESY']
config.Site.storageSite = 'T2_DE_DESY'
"""

with open(samples_file) as f:
    for line in f:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        request_name, input_dataset = line.split("|", 1)

        cfg_name = f"crab_{request_name}.py"

        with open(cfg_name, "w") as cfg:
            cfg.write(
                template.format(
                    request_name=request_name,
                    input_dataset=input_dataset
                )
            )

        print(f"Submitting {request_name}")
        subprocess.run(["crab", "submit", "-c", cfg_name], check=True)