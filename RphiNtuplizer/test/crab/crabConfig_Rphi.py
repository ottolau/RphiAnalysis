if __name__ == '__main__':

# Usage : python crabConfig.py (to create jobs)
#         ./multicrab -c status -d <work area> (to check job status)

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    from CRABClient.UserUtilities import config
    config = config()
    
    from multiprocessing import Process

    # Common configuration

    config.General.workArea     = 'crab_projects_ntuples'
    config.General.transferLogs = False
    config.JobType.pluginName   = 'Analysis' # PrivateMC
    config.JobType.psetName     = 'run_data_102X_aod.py'
    #config.JobType.inputFiles   = []
    config.JobType.sendExternalFolder = True
    config.Data.inputDBS        = 'global'    
    config.Data.splitting       = 'LumiBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
    #config.Data.splitting       = 'Automatic' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)

    config.Data.totalUnits      = -1
    config.Data.publication     = False
    config.Data.lumiMask        = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

    config.Site.storageSite     = 'T3_US_FNALLPC'
    #config.JobType.maxMemoryMB = 4000
    config.Data.allowNonValidInputDataset = True

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    # dataset dependent configuration

    config.General.requestName = 'ParkingBPH1_Run2018A-05May2019-v1_AOD_11Jul19_BsPhiLL_LowPtElectron_SculptingStudy'

    config.Data.unitsPerJob    = 3
    config.Data.inputDataset   = '/ParkingBPH1/Run2018A-05May2019-v1/AOD'

    config.Data.outLFNDirBase  = '/store/user/klau/LowPtElectron_SculptingStudy'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()




