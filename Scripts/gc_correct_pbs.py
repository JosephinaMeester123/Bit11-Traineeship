#!/usr/bin/python3

################################################################################
# command + arguments:
# python /home/jmeester/Internship/Scripts/cnv_script.py input_dir output_dir emailaddress

# example:
# python /home/jmeester/Internship/Scripts/gc_correct_pbs.py /home/jmeester/Internship/countfiles/ /home/jmeester/Internship/GC/ josephina.meester@uantwerpen.be

################################################################################
# Load dependencies
import os, time, sys
from subprocess import Popen, PIPE, STDOUT

################################################################################
# start time
t0 = time.time()

################################################################################
# Get input arguments
input_dir=sys.argv[1]
output_dir=sys.argv[2]
email=sys.argv[3]

print(input_dir)
print(output_dir)
print(email)

################################################################################
# Create a new folder
if not os.path.exists(output_dir):
        os.mkdir(output_dir)
#print("Changing directory to output_dir")
os.chdir(output_dir)
print("Current working directory: " + os.getcwd())

################################################################################
# Iterate over the list of filenames and directories
countfiles = os.listdir(input_dir)
if os.path.exists("/home/jmeester/Internship/GC/sample_list.txt"):
        os.remove("/home/jmeester/Internship/GC/sample_list.txt")
file_object=open("sample_list.txt", "a")
for line in countfiles:
        if line.endswith(".count"):
                print(line)
                file_object.write(input_dir + line + "\n")
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
                pbs_scriptname=sample + "_pbs_gc.sh"
                with open(pbs_scriptname, "w") as pbs_file:
                        pbs_file.write("#!/usr/bin/env bash \n")
                        pbs_file.write("#PBS -d /home/jmeester \n")
                        pbs_file.write("#PBS -l nodes=1:ppn=4 \n")
                        pbs_file.write("#PBS -l mem=12gb \n")
                        pbs_file.write("#PBS -l walltime=48:00:00 \n")
                        pbs_file.write("#PBS -N gc_correct \n")
                        pbs_file.write("#PBS -o " + output_dir + sample + ".gc_correct.stdout.$PBS_JOBID \n")
                        pbs_file.write("#PBS -e " + output_dir + sample + ".gc_correct.stderror.$PBS_JOBID \n")
                        pbs_file.write("#PBS -m abe \n")
                        pbs_file.write("#PBS -q batch \n")
                        pbs_file.write("#PBS -M " + email + "\n")
                        pbs_file.write("#PBS -V \n")
                        pbs_file.write("#PBS -A default \n")
                        pbs_file.write("/opt/software/miniconda_rnaseq/bin/Rscript /home/jmeester/Internship/Scripts/gcCorrectWGS.R " + filename + " " + sample + " " + output_dir)
			pbs_file.write("\n")
		Popen(["qsub", pbs_scriptname])

################################################################################
#sleep 3 sec to have enough time to submit all jobs
time.sleep(3)

# end time
t1 = time.time()
# print running time
print("Total time running: {} seconds".format(str(t1-t0)))

