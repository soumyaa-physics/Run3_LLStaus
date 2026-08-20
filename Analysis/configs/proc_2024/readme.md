## adding background MC using DAS

i used the following steps for adding background MC from DAS:

made a list of all files in a dataset in lxplus:

`voms-proxy-init --voms cms --rfc`

`dasgoclient -query="file dataset=/WW_TuneCP5_13p6TeV_pythia8/Run3Winter24NanoAOD-JMENanoV14_133X_mcRun3_2024_realistic_v10-v2/NANOAODSIM" > ww_files.txt`

and added a pre-fix:

`sed 's|^|root://cms-xrd-global.cern.ch///|' ww_files.txt > ww_xrd.txt`

then pasted the files in naf:

`scp Run3_Winter2024_WW_TuneCP5_xrd.txt vashisht@naf-cms.desy.de:/Run3_LLStaus/configs/datasets/2024`



how i searched on DAS:
`dataset=/*WW_TuneCP5_13p6TeV*/Run3*2024*/*`


`dataset=/WW_TuneCP5_13p6TeV_pythia8/Run3Winter24NanoAOD-JMENanoV14_133X_mcRun3_2024_realistic_v10-v2/NANOAODSIM`


Structure of NANOAOD DATASETS in DAS:

`/PrimaryDataset/ProcessingString/DataTier`

Chose this for now:

`dataset=/Muon1/Run2024C-2024CDEReprocessing-v1/NANOAOD`-> i think this is the final one after all reprocessing?
`dataset=/Muon1/Run2024D-2024CDEReprocessing-v1/NANOAOD`
`dataset=/Muon1/Run2024E-2024CDEReprocessing-v1/NANOAOD`

i could not find the other files: F,G,H,I-v1 and I-v2


root://cms-xrd-global.cern.ch////store/mc/Run3Winter24NanoAOD/WW_TuneCP5_13p6TeV_pythia8/NANOAODSIM/JMENanoV14_133X_mcRun3_2024_realistic_v10-v2/120000/82bbdfae-9fa4-429d-8cff-eff66b5ffc7b.root
