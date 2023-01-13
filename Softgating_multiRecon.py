import os
import glob
import shutil

#Number of Softgating Frames
noFrames = 8

#Get the P-file name to run the recon command
targetPattern = "*.7"
for name in glob.glob(targetPattern):
    Pfile = name

#Get current working directory
parentDir = os.getcwd()
fullpath = os.path.join

#Run the recon command in loop for all 8 frames 
for i in range(8):
    directory = 'Frame_0{}'.format(i)
    path = os.path.join(parentDir, directory)
    os.mkdir(path)
    #os.chdir(path)
    command = '~/local/bin/pcvipr_recon_binary -f {} -pils -dat_plus_dicom -external_gating_weights_name RespWeight_soft_0{}.dat -external_gating_weights -resp_gate_signal bellows -resp_gate thresh -resp_gate_efficiency .99999 -external_timestamps -external_timestamps_name external_gating_file'.format(Pfile,i)
    os.system(command)
    #os.chdir(parentDir)
    for dirname, dirnames, filenames in os.walk(parentDir):
        for filename in filenames:
            source = fullpath(dirname, filename)
            if filename.endswith("dcm"):
               shutil.move(source, path)
            elif filename.endswith("txt"):
                shutil.move(source, path)    
    printCom = 'Done reconstructing Frame {}'.format(i)
    print(printCom)