#!/usr/bin/env python

import sys, commands, os, fnmatch
from optparse import OptionParser

samplesDict = {}
samplesDict['HH'] = {
   #'BulkGravToWW_narrow_M_4000_13TeV_madgraph': 'mc',
   #'BulkGravToWW_narrow_M_4000_13TeV_madgraph_herwigpp': 'mc', 
   #'BulkGravToWW_narrow_M_600_13TeV_madgraph': 'mc',
   #'BulkGravToWW_narrow_M_600_13TeV_madgraph_herwigpp': 'mc',	
   #'BulkGravToWW_narrow_M_1000_13TeV_madgraph': 'mc', 
   #'BulkGravToWW_narrow_M_1000_13TeV_madgraph_herwigpp': 'mc', 
   #'BulkGravToWW_narrow_M_2000_13TeV_madgraph': 'mc', 
   #'BulkGravToWW_narrow_M_2000_13TeV_madgraph_herwigpp': 'mc', 
   #'BulkGravToWW_narrow_M_3000_13TeV_madgraph': 'mc', 
   #'BulkGravToWW_narrow_M_3000_13TeV_madgraph_herwigpp': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_1000_13TeV_madgraph':'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_1200_13TeV_madgraph':'mc',
   'BulkGravTohhTohbbhbb_narrow_M_1400_13TeV_madgraph':'mc',
   'BulkGravTohhTohbbhbb_narrow_M_1600_13TeV_madgraph':'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_1800_13TeV_madgraph':'mc',
   'BulkGravTohhTohbbhbb_narrow_M_2000_13TeV_madgraph': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_2500_13TeV_madgraph': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_3000_13TeV_madgraph': 'mc', 
   'BulkGravTohhTohbbhbb_narrow_M_3500_13TeV_madgraph': 'mc', 
}
samplesDict['VV'] = {
        'ZprimeToWW_narrow_M_1000_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_1200_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_1400_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_1600_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_1800_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_2500_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_3000_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_3500_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_4000_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_4500_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_600_13TeV_madgraph':'mc',
        'ZprimeToWW_narrow_M_800_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_1000_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_1200_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_1400_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_1600_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_1800_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_2500_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_3000_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_3500_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_4000_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_4500_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_600_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_narrow_M_800_13TeV_madgraph':'mc',
        'ZprimeToWWToWlepWhad_width0p1_M_1200_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p1_M_2000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p1_M_3000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p1_M_4000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p1_M_800_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p2_M_1200_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p2_M_2000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p2_M_3000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p2_M_4000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p2_M_800_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p3_M_1200_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p3_M_2000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p3_M_3000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p3_M_4000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWWToWlepWhad_width0p3_M_800_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p1_M_1200_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p1_M_2000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p1_M_3000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p1_M_4000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p1_M_800_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p2_M_1200_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p2_M_2000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p2_M_3000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p2_M_4000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p2_M_800_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p3_M_1200_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p3_M_2000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p3_M_3000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p3_M_4000_13TeV_madgraph_pythia8':'mc',
        'ZprimeToWW_width0p3_M_800_13TeV_madgraph_pythia8':'mc',
}
samplesDict['TT'] = {
        'ZprimeToTTJet_M_1000_13TeV_amcatnlo_pythia8_80X':'mc',
        'ZprimeToTTJet_M_1250_13TeV_amcatnlo_pythia8_80X':'mc',
        'ZprimeToTTJet_M_1500_13TeV_amcatnlo_pythia8_80X':'mc',
        'ZprimeToTTJet_M_2000_13TeV_amcatnlo_pythia8_80X':'mc',
        'ZprimeToTTJet_M_2500_13TeV_amcatnlo_pythia8_80X':'mc',
        'ZprimeToTTJet_M_3000_13TeV_amcatnlo_pythia8_80X':'mc',
        'ZprimeToTTJet_M_3500_13TeV_amcatnlo_pythia8_80X':'mc',
        'ZprimeToTTJet_M_4000_13TeV_amcatnlo_pythia8_80X':'mc',
        'ZprimeToTTJet_M_500_13TeV_amcatnlo_pythia8_80X':'mc',
        'ZprimeToTTJet_M_750_13TeV_amcatnlo_pythia8_80X':'mc',
}
samplesDict['QCD'] = {       
    'QCD_HT1000to1500_13TeV_80X':'mc',
    'QCD_HT1000to1500_13TeV_ext_80X':'mc',
    'QCD_HT1500to2000_13TeV_80X':'mc',
    'QCD_HT1500to2000_13TeV_ext_80X':'mc',
    'QCD_HT2000toInf_13TeV_80X':'mc',
    'QCD_HT2000toInf_13TeV_ext_80X':'mc',
    'QCD_HT200to300_13TeV_80X':'mc',
    'QCD_HT200to300_13TeV_ext_80X':'mc',
    'QCD_HT300to500_13TeV_80X':'mc',
    'QCD_HT300to500_13TeV_ext_80X':'mc',
    'QCD_HT500to700_13TeV_80X':'mc',
    'QCD_HT500to700_13TeV_ext_80X':'mc',
    'QCD_HT50to100_13TeV_80X':'mc',
    'QCD_HT700to1000_13TeV_80X':'mc',
    'QCD_HT700to1000_13TeV_ext_80X':'mc',
    }

samplesDict['MC'] = dict(samplesDict['HH'].items() +
		         samplesDict['TT'].items() +
			 samplesDict['VV'].items() +
                         samplesDict['QCD'].items())

#samplesDict['Data'] = dict(samplesDict['JetHT'].items() +
#                           samplesDict['SingleMuon'].items())

samplesDict['All'] = dict(samplesDict['MC'].items())

for label, isMc in samplesDict['All'].iteritems():
    samplesDict[label] = {label: isMc}
        
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
    parser.add_option('-t','--tag',dest="tag", default = "qbertbits-v13.2", help = "tag, which is the same as folder") 
    parser.add_option("--njobs-per-file",dest="njobs_per_file",type='int',default=1,help="Split into n jobs per file, will automatically produce submission scripts")
    
    (options,args) = parser.parse_args()

    monitorOption = ''
    if options.monitor is not '':
        monitorOption = '--monitor %s'%options.monitor
    
    jsonPrompt = "$PWD/../data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
    jsonRereco = "$PWD/../data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
    
    xsec = 1

    eosOutDir = '/eos/cms/store/group/phys_exotica/monotop/qbert/v3/'

    
    optionsDataMc = {
        'mc': "Output.root --passSumEntries 5:Events -a 6:subjob_i -a 7:%i -a 2:mc -a 3:none  -n 8000 -q 2nw --njobs-per-file %d"%(options.njobs_per_file,options.njobs_per_file),
        'data': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonPrompt,options.njobs_per_file),        
        'rereco': "Output.root -a 5:1 -a 6:subjob_i -a 7:%i -a 2:data -a 3:%s -n 8000 -q 1nd --njobs-per-file %d"%(options.njobs_per_file,jsonRereco,options.njobs_per_file)
    }
        
    #analysisDir = "zprimebits-v11.051"
    analysisDir = options.tag
    executable = "runQbert"
    
    samples = samplesDict[options.sample]

    exec_me('mkdir -p %s/%s'%(eosOutDir,analysisDir))  
    for label, isMc in samples.iteritems():
        exec_me('mkdir -p %s/%s/%s'%(eosOutDir,analysisDir,label))
        if isMc in ['data','rereco']:
            exec_me("python baconBatch.py %s %s -a 4:%f --list 1:../lists/production13/%s.txt --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s  %s"%(executable,optionsDataMc[isMc],xsec,label,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
        else:            
            exec_me("python baconBatch.py %s %s -a 4:%f --list 1:../lists/production13/%s.txt --outdir $PWD/../%s/%s_%s --eosoutdir %s/%s/%s  %s"%(executable,optionsDataMc[isMc],xsec,label,analysisDir,label,isMc,eosOutDir,analysisDir,label,monitorOption),options.dryRun)
