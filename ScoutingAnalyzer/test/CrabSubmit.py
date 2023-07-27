from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Run3Scouting2022F'
config.General.workArea = 'Run3Scouting2022F'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ScoutingForCrab_cfg.py'
config.JobType.maxMemoryMB = 2500

config.Data.inputDataset = '/ScoutingPFRun3/Run2022F-v1/RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 100
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ppradeep/Run3Scouting2022F'
config.Data.publication = True
config.Data.outputDatasetTag = '2022F'

config.Site.storageSite = 'T2_UK_London_IC'