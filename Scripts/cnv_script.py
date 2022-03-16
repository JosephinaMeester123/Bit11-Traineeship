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
# Loop over files in sample_list.txt
print("Looping over lines in sample_list")
with open ("sample_list.txt", "r") as sample_list:
	for filename in sample_list:
		filename=filename.strip()
		sample=filename.split("/")[-1]
		sample=sample.split(".")[0]
		print("processing the following sample: {}".format(sample))
		with open(sample + "_count.txt", "w") as output_file:
        		cmd1=Popen(["/opt/NGS/binaries/samtools/1.9/samtools", "view", "-@", "4", "-b", "--input-fmt-option", "required_fields=0x100", "-F", "0x100", "-T",
                		"/opt/NGS/References/hs38DH/genome/hs38DH.fa", filename, "chr1:1-500000"], stdout=PIPE)
        		cmd2=Popen(["/opt/NGS/binaries/BedTools/2.28.0/bin/coverageBed", "-F", "0.5001", "-sorted", "-counts", "-b", "stdin", "-a", "/home/jmeester/Internship/windows1kb.bed", "-g",
                		"/home/jmeester/Internship/genome.txt"], stdin=cmd1.stdout, stdout=(output_file))

################################################################################
# end time
t1 = time.time()
# print running time
print("Total time running: {} seconds".format(str(t1-t0)))
