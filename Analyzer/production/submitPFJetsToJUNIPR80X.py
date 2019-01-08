#!/usr/bin/env python

import sys, commands, os, fnmatch
from optparse import OptionParser

samplesDict = {}
samplesDict['QCD'] = {       
    'QCD_HT1000to1500_13TeV_80X':'mc',
    #'QCD_HT1000to1500_13TeV_ext_80X':'mc',
    #'QCD_HT1500to2000_13TeV_80X':'mc',
    #'QCD_HT1500to2000_13TeV_ext_80X':'mc',
    #'QCD_HT2000toInf_13TeV_80X':'mc',
    #'QCD_HT2000toInf_13TeV_ext_80X':'mc',
    #'QCD_HT200to300_13TeV_80X':'mc',
    #'QCD_HT200to300_13TeV_ext_80X':'mc',
    'QCD_HT300to500_13TeV_80X':'mc',
    #'QCD_HT300to500_13TeV_ext_80X':'mc',
    'QCD_HT500to700_13TeV_80X':'mc',
    #'QCD_HT500to700_13TeV_ext_80X':'mc',
    #'QCD_HT50to100_13TeV_80X':'mc',
    #'QCD_HT700to1000_13TeV_80X':'mc',
    #'QCD_HT700to1000_13TeV_ext_80X':'mc',
    }

        
def exec_me(command, dryRun=False):
    print command
    if not dryRun:
        os.system(command)
        
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")    
    parser.add_option("--monitor",default='',help="Monitor mode (sub/resub/check directory of jobs)")
    parser.add_option('-s','--sample',dest="sample", default="All",
                      #choices=['All','Hbb','QCD','JetHT','SingleMuon','DMSpin0','TT','DY','W','Diboson','Triboson','SingleTop','VectorDiJet1Jet','VectorDiJet1Gamma','MC','Data'],
                      help="samples to produces")
    parser.add_option('--production',dest="production", default = "13",
                      help="bacon production")
    parser.add_option('-t','--tag',dest="tag", default = "qbertbits-v13.2", help = "tag, which is the same as folder") 
    parser.add_option("--njobs-per-file",dest="njobs_per_file",type='int',default=1,help="Split into n jobs per file, will automatically produce submission scripts")
    parser.add_option("--nfiles-per-job", dest="nfiles_per_job", type='int', default=1,
                      help="Split into n files per job, will automatically produce submission scripts")    
    (options,args) = parser.parse_args()

    monitorOption = ''
    if options.monitor is not '':
        monitorOption = '--monitor %s'%options.monitor
    
    jsonPrompt = "$PWD/../data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
    jsonRereco = "$PWD/../data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
    

    eosOutDir = '/eos/cms/store/group/phys_exotica/dijet/dazsle/'

    optionsDataMc = {
        'mc': "-a 2:Output.root -a 3:0 -a 4:mc -n 8000 --njobs-per-file %d --nfiles-per-job %d"%(options.njobs_per_file,options.nfiles_per_job),
    }
        
    analysisDir = options.tag
    executable = "runPFJetsToJUNIPR"
    execPython = 'baconCondor.py'

    samples = samplesDict[options.sample]

    exec_me('mkdir -p %s/%s'%(eosOutDir,analysisDir))  
    for label, isMc in samples.iteritems():
        exec_me('mkdir -p %s/%s/%s'%(eosOutDir,analysisDir,label))
        listLabel = '../lists/production%s/%s.txt'%(options.production,label)
        exec_me("python %s %s %s --list 1:%s --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s %s"%(execPython,executable,optionsDataMc[isMc],listLabel,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)


