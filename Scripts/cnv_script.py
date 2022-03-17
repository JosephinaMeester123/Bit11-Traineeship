#!/usr/bin/python3
################################################################################
# Load dependencies
import os, time
from subprocess import Popen, PIPE, STDOUT

################################################################################
# start time
t0 = time.time()

################################################################################
# Create a new folder
newdir = "/home/jmeester/WGS_Jeannette/Testing/"
if not os.path.exists(newdir):
        os.mkdir(newdir)
print("Changing directory to WGS_Jeannette/Testing folder...")
os.chdir(newdir)
print("\nCurrent working directory: ", os.getcwd())

################################################################################
# Iterate over the list of filenames and directories
cram_dir = "/home/jmeester/WGS_Jeannette/Cram/"
cram_list = os.listdir(cram_dir)
print("Files and folders in {} :".format(cram_dir))
if os.path.exists("/home/jmeester/WGS_Jeannette/Testing/sample_list.txt"):
        os.remove("/home/jmeester/WGS_Jeannette/Testing/sample_list.txt")
file_object=open("sample_list.txt", "a")
for line in cram_list:
        if line.endswith(".cram"):
                print(line)
                file_object.write(cram_dir + line + "\n")
file_object.close()

################################################################################
# Loop over files in sample_list.txt to create PBS script per sample and perform counting
print("Looping over lines in sample_list")
with open ("sample_list.txt", "r") as sample_list:
        for filename in sample_list:
                filename=filename.strip()
                sample=filename.split("/")[-1]
                sample=sample.split(".")[0]
                print("processing the following sample: {}".format(sample))
                pbs_scriptname=sample + "_pbs.sh"
                with open(pbs_scriptname, "w") as pbs_file:
                        pbs_file.write("#!/usr/bin/env bash \n")
                        pbs_file.write("#PBS -d /home/jmeester/WGS_Jeannette \n")
                        pbs_file.write("#PBS -l nodes=1:ppn=4 \n")
                        pbs_file.write("#PBS -l mem=12gb \n")
                        pbs_file.write("#PBS -l walltime=48:00:00 \n")
                        pbs_file.write("#PBS -N WGS_count \n")
                        pbs_file.write("#PBS -o /home/jmeester/WGS_Jeannette/Testing/" + sample + ".count.stdout.$PBS_JOBID \n")
                        pbs_file.write("#PBS -e /home/jmeester/WGS_Jeannette/Testing/" + sample + ".count.stderror.$PBS_JOBID \n")
                        pbs_file.write("#PBS -m abe \n")
                        pbs_file.write("#PBS -q batch \n")
                        pbs_file.write("#PBS -M josephina.meester@uantwerpen.be \n")
                        pbs_file.write("#PBS -V \n")
                        pbs_file.write("#PBS -A default \n")
                        #pbs_file.write("echo 'Start Time : ' `date` \n")
                        pbs_file.write("/opt/NGS/binaries/samtools/1.9/samtools view -b --input-fmt-option required_fields=0x100 -F 0x100 -T /opt/NGS/References/hs38DH/genome/hs38DH.fa " + filename +$
                Popen(["qsub", pbs_scriptname])
                #copy output to results folder
################################################################################
# end time
t1 = time.time()
# print running time
print("Total time running: {} seconds".format(str(t1-t0)))

