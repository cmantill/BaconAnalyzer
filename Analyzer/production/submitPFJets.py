#!/usr/bin/env python
import sys, commands, os, fnmatch
from optparse import OptionParser

samplesDict = {}
samplesDict['JetHTrereco_8X_withPF'] = {
    'JetHTRun2016B_07Aug17_ver1_v1_8X_withPF': 'rereco16',
    # 'JetHTRun2016B_07Aug17_ver2_v1_8X_withPF': 'rereco16',
    # 'JetHTRun2016C_07Aug17_v1_8X_withPF': 'rereco16',
    # 'JetHTRun2016D_07Aug17_v1_8X_withPF': 'rereco16',
    # 'JetHTRun2016E_07Aug17_v1_8X_withPF': 'rereco16',
    # 'JetHTRun2016F_07Aug17_v1_8X_withPF': 'rereco16',
    # 'JetHTRun2016G_07Aug17_v1_8X_withPF': 'rereco16',
    # 'JetHTRun2016H_07Aug17_v1_8X_withPF': 'rereco16',
}

def exec_me(command, dryRun=False):
    print command
    if not dryRun:
        os.system(command)
        
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                      help="Just print out commands to run")    
    parser.add_option("--monitor",default='',
                      help="Monitor mode (sub/resub/check directory of jobs)")
    parser.add_option('-s','--sample',dest="sample", default="",
                      help="samples to produce")
    parser.add_option('-e','--executable',dest="executable", default="runPFJets", 
                      help = "executable name")
    parser.add_option('-t','--tag',dest="tag", default = "zprimebits-v15.01", 
                      help = "tag, which is the same as folder") 
    parser.add_option('--production',dest="production", default = "15",
                      help="bacon production") 
    parser.add_option('-b','--batch',dest="sub", default = False, 
                      help = "use condor or batch system")
    parser.add_option("--njobs-per-file",dest="njobs_per_file",type='int',default=1,
                      help="Split into n jobs per file, will automatically produce submission scripts")
    parser.add_option("--nfiles-per-job", dest="nfiles_per_job", type='int', default=1,
                      help="Split into n files per job, will automatically produce submission scripts")    
    (options,args) = parser.parse_args()

    monitorOption = ''
    if options.monitor is not '':
        monitorOption = '--monitor %s'%options.monitor
    
    jsonPrompt16 = "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
    jsonRereco16 = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
    jsonPrompt17 = "Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
    jsonRereco17 = "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"
    jsonPrompt18 = "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"

    analysisDir = options.tag
    executable = options.executable
    samples = samplesDict[options.sample]

    if options.sub:
        eosOutDir = '/eos/cms/store/group/phys_exotica/dijet/dazsle'
        execPython = 'baconBatch.py'
        EOS = ''
        optionsDataMc = {
            'mc': "Output.root -a 6:subjob_i -a 7:%i -a 2:mc -a 3:none  -n 8000 -q 2nw4cores --njobs-per-file %d"%(options.njobs_per_file,options.njobs_per_file),
            'data': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonPrompt16,options.njobs_per_file),        
            'rereco': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonRereco16,options.njobs_per_file),
            }
    else:
        eosOutDir = '/store/user/cmantill/'
        execPython = 'baconCondor.py'
        EOS = ''#eos root://cmseos.fnal.gov'
        optionsDataMc = {
            'rereco16': "-a 2:Output.root -a 3:0 -a 4:data -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,options.nfiles_per_job),
            }

    exec_me('%s mkdir -p /eos/uscms/%s/%s'%(EOS,eosOutDir,analysisDir))  
    for label, isMc in samples.iteritems():
        exec_me('%s mkdir -p /eos/uscms/%s/%s/%s'%(EOS,eosOutDir,analysisDir,label))
        listLabel = '../lists/production%s/%s.txt'%(options.production,label)
        exec_me("python %s %s %s --list 1:%s --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s %s"%(execPython,executable,optionsDataMc[isMc],listLabel,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
