#!/usr/bin/perl


#############################
## CMG NIPT Pipeline V 2.0 ##
#############################


# - Created by Geert Vandeweyer & Philip Holmgren
# - Date : 2014-10-01
# - Goal : Reproduce validated && accredited workflow of Genomics Core Leuven
# - Used (as) identical (as possible) versions of tools and reference data.
# - Switched to indivdual jobs instead of job-arrays on 2014-04-30 : maui does not handle them very well.
# - 2014-09-11 : Rewrite job submission to decrease the total number of jobs on the queue : intermediate jobs to submit next batch/step.
# - 2014-10-01 : Version 2.0 : Upgrade algorithm version, implement new binning methods (fixed size), fetal fraction
# - 2015-05-12 : Removed GATK steps.
# - 2015-06-08 : Added SeqFF method to calculate ff from autosomal regions
# - 2015-06-16 : Adapt to work with NextSeq data
# - 2015-08-03 : Replaced samtools view | filter | wc -l  principle to bedtools to improve performance. 
# - 2015-12-25 : 	1) One pipeline for HiSeq/NextSeq
#				 	2) Inegrated GetStats.pl code as subroutine in main script (easier to maintain and manage failed samples)
#					3) Improved error handling (changed if with while loop)
#					4) Set Job limit so not to overload PBS queue (manager only looks at first 4096 jobs)
#					5) Various code cleanup

# Future changes logged in hg repo


## order :
# 1. map all fastq files 
#	- bwa 0.7.5a : aln -q 15 -t 8 <ref.fasta> -f <sai file> <fastq.gz>
#	- bwa 0.7.5a : samse -r '@RG\tID:MedGen-<BC>-<DATE>_NIPT_V_1\tSM:<SAMPLE\tPL:ILLUMINA\tCN:MedGen_UZA' <ref.fasta> <SAI> <fastq.gz> -f <SAM>
# 2. Merge Lanes & Sort (primary bam 
#	- picard 1.78 MergeSam I=<BAM> O=/tmp/<BAM>
# 3. Generate per Chromosome sorted BAM, dedup and filter for unwanted reads
#	- samtools 0.1.19 : samtools view -uS -r 1 <SAM> | <grep's> | samtools sort - <bam> 
#	- picard 1.78 markduplicates
#	- samtools 0.1.19 : samtools index <BAM>
# 4. Count Reads in Bins
# 	- count / chr. (bedtools | sort)
# 5. Merge results into raw count file (cat)
# 6. GC Correct
# 	- gcCorrectNIPT.Rscript (Rscript)
# 7. Sliding Window Counts (WinBinFixedSize)
#	- getWinBinFixedSize (Rscript)
# 8. Analyze data
#	- get fetal fraction
#	- Generate Plots => more specific sex-chr handling.
# 9. Generate report (not on HPC)
# 10. Wrap up (failed samples)
# 11. Backup data on LTS



$|++;

############
## MODULES  #
############
use Getopt::Std;
use Cwd 'abs_path';
use Time::Piece;
use Data::Dumper;
use XML::Simple;
use DBI;
use Digest::MD5 'md5_hex';
use File::Basename;
use File::Spec::Functions 'splitdir';
#no warnings 'experimental::smartmatch';


##########################
# COMMAND LINE ARGUMENTS #
##########################
# mandatory : datadir, configuration file
my %opts;
getopts('d:c:', \%opts);  # option are in %opts
# mandatory : datadir
if (!defined($opts{'d'}) || $opts{'d'} eq '' || !-d $opts{'d'}) {
	dieWithGrace('Datadir not specified or not accessible', 'Not Checked');
}
my $datadir = $opts{'d'};
if (!-d "$datadir") {
	dieWithGrace("DataDir $datadir does not exist", 'Not Checked');
}
$datadir = abs_path("$datadir");
my @dp = split(/\//,$datadir);
my $rundate = $dp[-1];
my $project = $dp[-1];
my $RunDate = 'n/a';
if ($dp[-1] =~ m/^(\d+)_.*/) {
	$rundate = $1;
	$RunDate = '20'.substr($rundate,0,2).'-'.substr($rundate,2,2).'-'.substr($rundate,4,2);
}

if (! defined($opts{'c'}) || $opts{'c'} eq '') {
	dieWithGrace("XML file is not given", $project);
}
my $configfile = $opts{'c'};
$configfile = abs_path("$configfile");
if (!-e $configfile) {
	dieWithGrace("XML config file '$configfile' does not exist", $project);
}
my $xml = new XML::Simple;
my $config = $xml->XMLin($configfile);


######################
# CONFIGURATION FILE #
######################
my $locations = $config->{locations};
my $basedir = $locations->{basedir};           # /home/nipt/NIPT_Validated/
my $refdatadir = $locations->{refdatadir};     # /opt/NGS/References/nipt/
my $basedatadir = $locations->{basedatadir};   # /home/nipt/Run_Data/
my $toolsdir = $locations->{ngstools};	       
my $softwaredir = $locations->{software};	
my $destinations = $config->{destinations};
my $sambadir = $destinations->{ict_uza};       # /home/ict_uza
my $users = $config->{users};	
my $email = $config->{email};
my $sender = $email->{sender};                 # NIPT Pipeline
my $addressed = $email->{addressed};           # katrien.janssens@uantwerpen.be
my $admin = $email->{admins};                  # geert.vandeweyer@uantwerpen.be,philip.holmgren@uantwerpen.be
my $debug = $config->{debug};                  # In case of debugging (1) sending data to transix server is off
my $rerun = $config->{rerun};                  # If 1: do perform reruns
my $mysql = $config->{mysql};		       # If 0: do not enter data in database
my $database = $config->{database};
my $userid = $database->{userid};
my $userpass = $database->{userpass};
my $db = $database->{dbname};
my $host = $database->{host};
my $mainversion = $config->{version};
my $niptqueue = $config->{queue};
my $niptaccount = $config->{hpcaccount};

#####################
# DETERMINE SEQTYPE #
#####################
my ($seqtype, $abbr);
if ($project =~ /^\d*_SNL/) {
	$seqtype = "HiSeq";
	$abbr = "HS";
}
elsif ($project =~ /^\d*_N/) {
	$seqtype = "NextSeq";
	$abbr = "NS";
}

###################
# HARDCODED FILES #
###################
my $ss_file = "$datadir/SampleSheet.csv";
my $redo_file = "$datadir/Redo.txt";
my $rerun_file = "$datadir/Rerun.txt";
my $status_file = "$datadir/Status.txt";
my $runinfo = "$datadir/RunInfo.xml";
my $runparam = "$datadir/RunParameters.xml";
my $runtimeoutput = "$datadir/RunTime.Output.txt";
my $samplesdone = "$basedatadir/Samples_Done.csv";
my $unusedcountsfile = "$datadir/Unused.counts.txt";
my $checksumfile = "$datadir/checksums.databuffer.md5";
my $checksumresult = "$datadir/checksums.result.stdout";
my $checksumerrorresult = "$datadir/checksums.result.stderr";
my $checksumLTSReport = "$datadir/checksums.LTS.log";
my $checksumLTSReportfile = "$datadir/checksums.LTS.Report.md5";
my $checksumLTSDatafile = "$datadir/checksums.LTS.Data.md5";
my $qcfolder = "$datadir/InterOp";
my $UBfolder = "$datadir/UnusedIndices";
my $txt = "$datadir/Run.QC.txt";
my $pdf = "$datadir/Run.QC.pdf";
my $cnvreport = "$datadir/CNVs.txt";
my $failedjobs = "$datadir/tmp_files/Failed_Jobs.txt";
my $fullproject = $abbr."_".$project;
my $ltsreportsdir = "$sambadir/CMG/NIPT2/Reports/$fullproject";
my $ltsdatadir = "$sambadir/CMG/NIPT2/Data/$fullproject/";
my $original_project;
my $original_datadir;

# get original project for ReDo runs
if (-e $redo_file) {
	$original_project = `cat $redo_file`;
	chomp($original_project);
	$original_datadir = "$basedatadir/$original_project";
	$mysql = 0;
}

#######################
# HARDCODED VARIABLES #
#######################
## Reference data (NextSeq refset are 75bp unless specified differently!)
my $ref;
my %referencedatadir = (
#	"ns_chip50"		=> "$basedir/Reference_Data/Reference_Data_NS_ChIP_50bp",
#	"ns_chip" 		=> "$basedir/Reference_Data/Reference_Data_NS_ChIP",
#	"ns_nano"   		=> "$basedir/Reference_Data/Reference_Data_NS_Nano",
#	"ns_nano_hamilton_old"	=> "$basedir/Reference_Data/Reference_Data_NS_Nano_Hamilton",
	"ns_nano_hamilton_xy"	=> "$basedir/Reference_Data/Reference_Data_NS_Nano_Hamilton_XYCorrected",
	"ns_nano_hamilton"	=> "$basedir/Reference_Data/Reference_Data_NS_Nano_Hamil_XY_tris2",
#	"ns_nano_hamilton"	=> "/home/niptdev/Run_Data/Ref13Map/full_map/ReferenceSet_TRISFINAL",
#	"hs_chip"   		=> "$basedir/Reference_Data/Reference_Data_HS_ChIP",
#	"hs_nano"   		=> "$basedir/Reference_Data/Reference_Data_HS_Nano"
);
my %referencekit = (
#		"ns_chip50"		=> "NextSeq (50bp reads)+ ChIP",
#		"ns_chip" 		=> "NextSeq  + ChIP",
#		"ns_nano"   		=> "NextSeq + Nano (Manuele prep)",
#		"ns_nano_hamilton_xy"	=> "NextSeq + Nano (Hamilton) (XY)",
		"ns_nano_hamilton"	=> "NS Nano (Hamil, XY, triscorr2)",
#		"ns_nano_hamilton_old"	=> "NextSeq + Nano (Hamilton prep)",
#		"hs_chip"   		=> "HiSeq (50bp reads) + ChIP",
#		"hs_nano"   		=> "HiSeq (50bp reads) + Nano"
);

## Samplesheet structure
my %ssind = (
	NextSeq => {
		sample => {
			8 => 0,
			10 => 0
		},
		bc1 => {
			8 => 5,
			10 => 5
		},
		bc2 => {
			8 => '',
			10 => 7
		},
		description => {
			8 => 7,
			10 => 9
		}
	},

	HiSeq => {
		sample => {
			9  => 1,
			11 => 1
		},
		bc1 => {
			9  => 6,
			11 => 6
		},
		bc2 => {
			9  => '',
			11 => 8
		},
		description => {
			9 => 8,
			11 => 10
		}
	}



);


## Calender
my %month = ();
$month{'Jan'} = '01';
$month{'Feb'} = '02';
$month{'Mar'} = '03';
$month{'Apr'} = '04';
$month{'May'} = '05';
$month{'Jun'} = '06';
$month{'Jul'} = '07';
$month{'Aug'} = '08';
$month{'Sep'} = '09';
$month{'Oct'} = '10';
$month{'Nov'} = '11';
$month{'Dec'} = '12';



####################
# PIPELINE VERSION #
####################
my $pversion = `cd $basedir && git log -1 --format=\%cd`;
my $pchangeset = `cd $basedir && git log -1 --format=\%h`;
chomp($pversion);
chomp($pchangeset);
$pversion =~ m/(\S{3})\s(\S{3})\s(\d{1,2})\s(\d{2}:\d{2}:\d{2})\s(\d{4})\s.*/;
$pversion = "$5-".$month{$2}."-$3";
$pversion .= "-$pchangeset";
my $version = "_NIPT_V".$mainversion."_".$pversion;



#########
# START #
#########
my $user = `id -un`;
chomp($user);
print "Running script as : '$user'\n";
my $starttimeObj = Time::Piece->new();
my $starttime = $starttimeObj->epoch();
my $date = $starttimeObj->strftime('%Y-%m-%d');


####################
# GATHER UNUSED BC #
####################
if(-e $UBfolder && -d $UBfolder){
	if(-e $unusedcountsfile){
		system("mv $unusedcountsfile $unusedcountsfile.bak");
	}
	
	# Make Unused.counts.txt
	print "Making Unused.counts.txt\n";
	my @unusedFastQFiles = `cd $UBfolder && ls *.fastq.gz`;
	chomp(@unusedFastQFiles);
	foreach $file (@unusedFastQFiles) {
		$file =~ /(.*)_S\d*_L\d*_R\d_/;
		my $pname = $1;
		my $nrl = `cd $UBfolder/ && zcat $file | wc -l `;
		chomp($nrl);
		my $nrr = $nrl / 4;
		system("echo '$pname,$file,$nrr' >> $unusedcountsfile");
	}
	
	if(	!-e "$datadir/unused_barcodes_S0_L001_R1_001.fastq.gz" ||
		!-e "$datadir/unused_barcodes_S0_L002_R1_001.fastq.gz" ||
		!-e "$datadir/unused_barcodes_S0_L003_R1_001.fastq.gz" ||
		!-e "$datadir/unused_barcodes_S0_L004_R1_001.fastq.gz"){

		# Make Unused fastq.gz
		print "Making unused barcode fastqs\n";
		
		system("zcat $UBfolder/unused_*_L001_R1_001.fastq.gz | gzip -c > $datadir/unused_barcodes_S0_L001_R1_001.fastq.gz");
		system("zcat $UBfolder/unused_*_L002_R1_001.fastq.gz | gzip -c > $datadir/unused_barcodes_S0_L002_R1_001.fastq.gz");
		system("zcat $UBfolder/unused_*_L003_R1_001.fastq.gz | gzip -c > $datadir/unused_barcodes_S0_L003_R1_001.fastq.gz");
		system("zcat $UBfolder/unused_*_L004_R1_001.fastq.gz | gzip -c > $datadir/unused_barcodes_S0_L004_R1_001.fastq.gz");
		
		if (!-e "$checksumfile") {
		        dieWithGrace("No checksumfile found, can't check FASTQ data integrity", "$project");
		} else {
			system("cd $datadir && md5sum unused_barcodes*.fastq.gz >> $checksumfile");
		}
	}
} 

########################
# CHECK DATA INTEGRITY #
########################
if (!-e "$checksumfile") {
	dieWithGrace("No checksumfile found, can't check FASTQ data integrity", "$project");
}

# Do checksum comparison
system("cd $datadir && md5sum -c $checksumfile 1> $checksumresult 2> $checksumerrorresult");

if (-s $checksumerrorresult) {
	my $output = `cat $checksumresult $checksumerrorresult`;
	dieWithGrace("POSSIBLE DATA CORRUPTION:\n$output\nCheck data manually!!", "$project");
}
else {
	print "\nCHECKSUM ANALYSIS: OK ($checksumresult)\n\n";
}


#####################
# PREPARE RUNFOLDER #
#####################
my $backupfolder;
# Delete and make all directories
my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
my $date_suffix = "$year-$mon-$day.$hour-$min-$sec";
## Make backup of old data
if (-e $status_file) {
	$backupfolder = "Old_analysis.$date_suffix";
	system("mkdir $datadir/$backupfolder");
	print "Backing up old data to $backupfolder\n";
	system("cd $datadir && mv CNVs.txt Run.QC.* Job_Output tmp_files tmp_binaries Results $backupfolder");
}
system("cd $datadir && mkdir Job_Output tmp_files tmp_binaries Results");

## set status.
system("echo 0 > $status_file");


###############
# GET SAMPLES #
###############
my %Samples = ();

if (!-e $ss_file) {
	dieWithGrace("No Samplesheet Found", "$project");
}
system("dos2unix '$ss_file'");

print "Loading samples:\n";

my $nrofastq; # for next step: PBS queue max calculation
my $start = 1; # switch
open IN, $ss_file;
while (<IN>) {
	my $line = $_;
	chomp($line);
	next if ($line eq "");
	next if (substr($line,0,1) eq "#");

	if ($line =~ /\[Data\]/) {
		$start = 0;
		# skip header
		<IN>;
		next;
	}
	# skip lines before data lines
	if ($start) {
		next;
	}

	my @sscolumns = split(',', $line, -1);
	my $nrcols = scalar(@sscolumns);
	my $sample = lc($sscolumns[$ssind{$seqtype}{sample}{$nrcols}]);
	my ($bc, $bc2);
	my $bc1 = lc($sscolumns[$ssind{$seqtype}{bc1}{$nrcols}]);
	if ($ssind{$seqtype}{bc2}{$nrcols}) {
		$bc2 = lc($sscolumns[$ssind{$seqtype}{bc2}{$nrcols}]);
		$bc = $bc1."-".$bc2;
	}
	else {
		$bc = $bc1;
	}

	$ref = lc($sscolumns[$ssind{$seqtype}{description}{$nrcols}]);

	# to grep fastq, bam files (we drop the 'R' in the next step)
	my $gsample = $sample;

	# rerun?
	if ($sample =~ /(\d+)r$/i) {
		my $sampler = $sample;
		$sample = $1;
		my $rerunsample = $sample;
		$Samples{$sample}{'rerun'} = 1;
		$Samples{$sample}{'rerunname'} = $sampler;
	}
	else {
		$Samples{$sample}{'rerun'} = 0;
	}

	# avoid processing samples twice with 2 lanes
	if ($Samples{$sample}{'seen'}) {
		next;
	}
	else {
		$Samples{$sample}{'seen'} = 1;
		# this replaces the failed sampleshash, turns 1 if sample fails
		$Samples{$sample}{'failed'} = 0;

	}

	print "- $sample\n";

	mkdir("$datadir/Job_Output/$sample");
	mkdir("$datadir/tmp_binaries/$sample");
	mkdir("$datadir/tmp_files/$sample");
	mkdir("$datadir/Results/$sample");

	# get input files
	my @fastqs = split(/\n/, `cd $datadir && ls | grep -i '^$gsample\_.*\.fastq\.gz'`);
	$nrofastq = scalar(@fastqs);
	my $fastq = join('|', @fastqs);
	$Samples{$sample}{'fastqs'} = $fastq;

	if (! $fastq) {
		dieWithGrace("No FASTQ files found for sample $sample\n", $project);
	}

	$Samples{$sample}{'bc'} = $bc;
	$Samples{$sample}{'bc1'} = $bc1;
	$Samples{$sample}{'bc2'} = $bc2;
	$Samples{$sample}{'fullid'} = "MedGen-$sample-$bc-$rundate".$version;
}

close IN;

# Stop project if reference set is not filled in correctly (choosing a default is probably not a good idea)
if (! ($ref ~~ [keys %referencedatadir])) {
		dieWithGrace("Reference set '$ref' is not valid", $project);
}

# Add unused barcode sample if fastqs are present
if(     -e "$datadir/unused_barcodes_S0_L001_R1_001.fastq.gz" &&
        -e "$datadir/unused_barcodes_S0_L002_R1_001.fastq.gz" &&
        -e "$datadir/unused_barcodes_S0_L003_R1_001.fastq.gz" &&
        -e "$datadir/unused_barcodes_S0_L004_R1_001.fastq.gz"){

	my $sample = "unused_barcodes";

	print "- $sample\n";

	$Samples{$sample}{'rerun'} = 0;
	$Samples{$sample}{'seen'} = 1;
	$Samples{$sample}{'failed'} = 0;

	my $bc1 = my $bc2 = "XXXXXXXX";
	my $bc = $bc1."-".$bc2;

        $Samples{$sample}{'bc'} = $bc;
        $Samples{$sample}{'bc1'} = $bc1;
        $Samples{$sample}{'bc2'} = $bc2;
        $Samples{$sample}{'fullid'} = "MedGen-$sample-$bc-$rundate".$version;

	my $fastq = join("|",(
				"unused_barcodes_S0_L001_R1_001.fastq.gz",
				"unused_barcodes_S0_L002_R1_001.fastq.gz",
				"unused_barcodes_S0_L003_R1_001.fastq.gz",
				"unused_barcodes_S0_L004_R1_001.fastq.gz"
			     )
			);

	$Samples{$sample}{'fastqs'} = $fastq;

	mkdir("$datadir/Job_Output/$sample");
	mkdir("$datadir/tmp_binaries/$sample");
	mkdir("$datadir/tmp_files/$sample");
	mkdir("$datadir/Results/$sample");

}




####################################
# Calculate max jobs for PBS queue #
####################################
## during first loop most jobs are submitted
##    - instant submission: bwa.aln, bwa.samse, mkpBAM; # jobs depend on # fastq files : (# of fastq * 2) + 1
##    - afterwards, in sequence: 24 Split -> 24 DeDup -> 24 Count
##  --> each sample will have a maximum of 24 non-completed jobs at the same time

## to not overload the PBS queue (hard max of 4096 jobs) the number of sample submissions must be limited, set to 50 samples simultaneous
my $initial_max_pbs = 1200;
my $per_sample = 24;
my $sample = $initial_max_pbs / $per_sample;
my $max_pbs = $sample * ((2 * $nrofastq) + 1);


############
# DATABASE #
############

my $connectionInfo="dbi:mysql:$db:$host";
my $dbh;
my $pid;
print "Connecting to database $db on $host\n";
$dbh = DBI->connect($connectionInfo,$userid,$userpass) ;
## retry on failed connection (server overload?)
$i = 0;
while ($i < 10 && ! defined($dbh)) {
	sleep 7;
	$i++;
	print "Connection to $host failed, retry nr $i/10\n";
	$dbh = DBI->connect($connectionInfo,$userid,$userpass) ;
}
$dbh->{mysql_auto_reconnect} = 1;
$dbh->{PrintError} = 1;
$dbh->{RaiseError} = 1;
$dbh->{HandleError} = sub{ dieWithGrace(shift(),$project); };

my $sampleCount = 0;

if($mysql){
	print "Entering Project $project\n";
	my $projquery = "INSERT INTO `Projects`";
	$projquery   .= "(name,date_run,date_analysis,
				date_analysis_start,rerun,path_archive_data,
				path_archive_reports,code_revision)";
	$projquery   .= "VALUES";
	$projquery   .= "(?,?,?,?,?,?,?,?)";


        (my $ltsstordatadir = $ltsdatadir) =~ s|seqpilot/{1,2}Backup_Data|storageshare|;
        (my $ltsstorreportsdir = $ltsreportsdir) =~ s|seqpilot/{1,2}Backup_Data|storageshare|;

	
	my $projsth = $dbh->prepare($projquery);


	my @projvalues = ($project,$RunDate,$date,
			 $starttimeObj->strftime("%y-%m-%d %H:%M:%S"),
			 -e $rerun_file ? 1 : 0,
			 $ltsstordatadir,
			 $ltsstorreportsdir,
			 $pversion);

	$projsth -> execute(@projvalues);
	$projsth -> finish;
	$pid = $dbh->last_insert_id( undef, undef, undef, undef );


	  #Samples
	my $samplequery = "INSERT INTO `Samples` ";
	$samplequery   .= "(dna_nr,bc1,bc2,pid,reference_set,blanco,unused) ";
	$samplequery   .= "VALUES ";
	$samplequery   .= "(?,?,?,$pid,?,?,?)";
	my $samplesth = $dbh->prepare($samplequery);

	  #Ref
	my $refquery = 'SELECT rid FROM `Referencesets` WHERE samplesheet_tag = ?';
	my $refsth = $dbh->prepare($refquery);
	$refsth->execute($ref);
	my ($refid) = $refsth->fetchrow_array();
	
	  #BC
	my $bcseqquery = 'SELECT bid FROM `Barcodes` WHERE sequence = ?';
	my $bcseqsth = $dbh->prepare($bcseqquery);
	print "Inserting samples:\n";
	foreach my $sample (keys %Samples){
	#insert samples that are not validation sample into DB

		next unless($sample =~ m{^[0-9]{5}$} 
				|| $sample =~ m{blanco|unused_barcodes}i);
		$sampleCount++;

		print "-$sample\n";
		#--------- Sample ---------#
	
		# BC 1
		$bcseqsth->execute(uc($Samples{$sample}{"bc1"}));
		my ($bc1id) = $bcseqsth->fetchrow_array();
			
		# BC 2
		$bcseqsth->execute(uc($Samples{$sample}{"bc2"}));
		my ($bc2id) = $bcseqsth->fetchrow_array();	

		my @samplevalues = ($sample,
				    $bc1id,
				    $bc2id,
				    $refid);
		if($sample =~ m{blanco}i){	
			$samplesth -> execute( (@samplevalues,1,0) );
		} elsif($sample =~ m{unused_barcodes}i){	
			$samplesth -> execute( (@samplevalues,0,1) );
		} else {
			$samplesth -> execute( (@samplevalues,0,0) );
		}
		$Samples{$sample}{'sid'} = $dbh->last_insert_id( undef, undef, undef, undef );
	}

	$bcseqsth	-> finish;
	$samplesth	-> finish;

}

#################
# RERUN SAMPLES #
#################
print "Checking for Rerun Samples\n";
my $newdatadir = $datadir."_ReRuns";
if ( -d $newdatadir) {
	system("mv $newdatadir $newdatadir.Old.$date_suffix");
}

# to rename L001|2 -> L005|6; not existing lane numbers to avoid file name collision
my %rename = ();
$rename{'1'} = '5';
$rename{'2'} = '6';
$rename{'3'} = '7';
$rename{'4'} = '8';


  #Rerunsample selection
my $rerunselectquery = "SELECT p.name,p.date_run,s.sid,p.path_archive_data ";
$rerunselectquery   .= "FROM Projects AS p ";
$rerunselectquery   .= "JOIN Samples AS s USING (pid) ";
$rerunselectquery   .= "WHERE p.rerun = 0 AND p.pid <> $pid AND s.dna_nr = ? ";
$rerunselectquery   .= "ORDER BY p.date_run DESC LIMIT 1;";
my $rerunselectsth = $dbh->prepare($rerunselectquery);

  #Rerun line insert
my $reruninsert1query = "INSERT INTO `Reruns` ";
$reruninsert1query   .= "(original_sample_1,original_sample_2) ";
$reruninsert1query   .= "VALUES ";
$reruninsert1query   .= "(?,?)";
my $reruninsert1sth = $dbh->prepare($reruninsert1query);


my $finalize_rerun = 0;
foreach $sample (keys %Samples) {
	if ($Samples{$sample}{'rerun'} && $rerun) {
		my $sampler = $Samples{$sample}{'rerunname'};
		
		$rerunselectsth -> execute ($sample);
		my ($rerunprojectname,$rerundaterun,$rerunsample1id,$rerunarchivepath) = $rerunselectsth->fetchrow_array();

		$project_original = "$basedatadir/$rerunprojectname";


		if (! $rerunsample1id) {
			print "\tCan't find $sample in a previous run; skipping sample $sample for rerun...\n";
		} else {

			if (! -d $newdatadir) {
				system("mkdir $newdatadir");
				print "\nProcess rerun samples:\n----------------------\n"; 

				print "- Create samplesheet header\n";
				my $headernumber = `grep -n 'Sample_ID' $datadir/SampleSheet.csv | cut -d: -f1`;
				chomp($headernumber);
				system("head -n $headernumber $datadir/SampleSheet.csv > $newdatadir/SampleSheet.csv");
			}
			print "- $sample:\n";

			$reruninsert1sth -> execute( $rerunsample1id, $Samples{$sample}{"sid"} );
			$rerunID = $dbh->last_insert_id(undef,undef,undef,undef);
			# insert into line with rerunID the composed sample ID when made
		
			open(my $rerunFH,">>","$newdatadir/Rerun.txt") 
				|| dieWithGrace("Couldn't open $newdatadir/Rerun.txt: $!",$project);
			print $rerunFH "$sample\t$rerunID\n";
			close($rerunFH);


			if (! -d $project_original) {
				print "\t$project_original doesn't exist anymore; taking from archive\n";
				$project_original = $rerunarchivepath;
				# symlink from archive
			}

			$finalize_rerun	= 1;
			## Make new samplesheet: header + exact lines original project + lines new project & delete the 'R' in sample name
			print "\t- Add to new samplesheet\n";
			system("grep -i $sampler $datadir/SampleSheet.csv | sed \"s/$sampler/$sample/gi\" >> $newdatadir/SampleSheet.csv");
			print "\t- Copy original data: $project_original --> $newdatadir\n";

			## find the original fastq files (awk for full paths) and copy
			#my @fqs_original = split(/\n/,`cd $project_original && ls | grep -i '^$sample\_.*\.fastq.gz' | awk -vpath=\$PWD/ '{print path\$1}'`);
			my @fqs_original = split(/\n/,`cd $project_original && ls | grep -i '^$sample\_.*\.fastq.gz'`);
			foreach my $fastq (@fqs_original) {
				print "\t  $fastq\n";
				system("cp $project_original/$fastq $newdatadir");
				system("grep -i $fastq $project_original/".basename($checksumfile)."  >> $newdatadir/".basename($checksumfile));

			}
			## copy and rename the fastqs from the current project
			my @fqs_new = split(/\|/, $Samples{$sample}{'fastqs'});
			print "\t- Change lane number for new FASTQ files and copy to $newdatadir\n";
			foreach my $fastq (@fqs_new) {
				# caps because original file
				$fastq =~ m/^(\S*)_(S\d*)_L00(\d{1})_(R.*)/;
				my $fastq_new = $sample."_".$2."_L00".$rename{$3}."_".$4;
				print "\t  $fastq --> $fastq_new \n";
				system("cp $datadir/$fastq $newdatadir/$fastq_new");
				system("grep -i $fastq $checksumfile | sed \"s/$fastq/$fastq_new/gi\" >> $newdatadir/".basename($checksumfile));
			}

		}
	}

}

$rerunselectsth -> finish;
$reruninsert1sth -> finish;

# is false if all the original project are already archived
if ( $finalize_rerun) {
	print "- Copy additional files (Unused counts, XMLs, SAV data) to $newdatadir\n";
	system("cd $datadir && rsync -a Unused.counts.txt InterOp Run*.xml $newdatadir");
	system("cd $datadir && grep -P 'InterOp|Run.*xml' $checksumfile >> $newdatadir/".basename($checksumfile));
	system("cd $newdatadir && md5sum ./SampleSheet.csv >> $newdatadir/".basename($checksumfile));

	system("cd $datadir && rsync -a transfer.databuffer.done $newdatadir");   #needed to start automatic run
}

#if this is rerun project
if(-e $rerun_file){
	open(my $rerunFH,"<","$rerun_file") 
		|| dieWithGrace("Couldn't read $rerun_file",$project);
	while(my $rerunline = <$rerunFH>){
		chomp($rerunline);
		my ($sample,$ID) = split("\t",$rerunline);
		$Samples{$sample}{'rerunID'} = $ID;
	}
	close($rerunFH);
}


###################
# PROCESS SAMPLES #
###################

print "\nStart Analysis:--------------\n";

foreach $sample (keys(%Samples)) {
	&CheckErrors(); # the error checking routine. Updates the failed samples hash
	#&CheckQueue();
	print " - Sample : $sample\n";

	print "Generating BAM File\n";
	my @fqs = split(/\|/,$Samples{$sample}{'fastqs'});
	my ($sai, $sam, $bam, $bai, $samseids, $mergein);
	my $chrfullbam = "$sample.\$CHR.full.bam";
	my $chrfullbai = "$sample.\$CHR.full.bai";
	my $chrbam = "$sample.\$CHR.bam";
	my $chrbai = "$sample.\$CHR.bai";
	my $chrdedupbam = "$sample.\$CHR.dedup.bam";
	my $chrdedupbai = "$sample.\$CHR.dedup.bai";
	my $primarybam = "$sample.bam";
	my $primarybai = "$sample.bai";
	my $finalbam = "$sample.final.bam";
	my $finalbai = "$sample.final.bai";
	foreach my $fastq (@fqs) {
		my $rgid = ''; # unique per lane.
		print "\t - $fastq\n";
		## extract lane and bc from fastq name: rest from %Samples
		$fastq =~ m/^(.*)_(s\d*)_l(\d{3})_r/i;
		my $lane = "L$3";
		my $bc = $Samples{$sample}{'bc'};
		# prepare filenames
		$rgid = "MedGen-$sample-$bc-$lane-$rundate".$version;
		$sai = $sam = $bam = $bai = $fastq;
		## the following are specific for the fastq file, as they are based on the name. no need for rgid-based naming.
		$sai =~ s/\.gz$/\.sai/;
		$sam =~ s/\.gz$/\.sam/;
		$bam =~ s/\.gz$/\.bam/;
		$bai =~ s/\.gz$/\.bai/;
		my $root = substr($bam,0,-4);

		## 1a. Map fastq : bwa aln
		my $scriptname = "$sample.bwa.aln.$lane";
		my $script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
		my $jobname = "NIPT.".$scriptname;
		open OUT, ">$script";
		#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
		&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 8, 12, $niptqueue, $niptaccount, undef);

		&PBSWriteCommand("OUT", "$toolsdir/bwa/0.7.5a/bin/bwa aln -q 15 -t 8 $refdatadir/human_g1k_v37.fasta -f /tmp/$project.$rgid.sai $datadir/$fastq", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (aln)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$rgid.sai $datadir/tmp_files/$sample/$sai", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy)", $failedjobs);

		&PBSWriteCommand("OUT", "rm -f /tmp/$project.$rgid.sai", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);
		&PBSWriteEnd("OUT");
		close OUT;

		my $alnid = `qsub $script`;
		chomp($alnid);
		while ($alnid !~ m/^\d+\..*/) {
			sleep 5;
			my $tmpfile = `mktemp`;
			chomp($tmpfile);
			my $return = `qstat -x 2>$tmpfile | grep $script`;
			chomp($return);
			my $errorreturn = `cat $tmpfile`;
			if (! ($return or $errorreturn)) {
				$alnid = `qsub $script`;
				chomp($alnid);
			}
			system("rm $tmpfile");
		}
		chomp($alnid);
		$alnid =~ s/(\d+)\..*/$1/;



		## 1b. Map fastq : bwa samse + convert to indexed bam
		$scriptname = "$sample.bwa.samse.$lane";
		$script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
		$jobname = "NIPT.".$scriptname;
		open OUT, ">$script";
		#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
		my $depend = ["#PBS -W depend=afterok:$alnid"];
		&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 1, 8, $niptqueue, $niptaccount, $depend);

		&PBSWriteCommand("OUT", "$toolsdir/bwa/0.7.5a/bin/bwa samse -r '\@RG\\tID:$rgid\\tSM:$sample\\tPL:ILLUMINA\\tCN:MedGen_UZA\\tDT:$date' $refdatadir/human_g1k_v37.fasta $datadir/tmp_files/$sample/$sai $datadir/$fastq -f /tmp/$project.$rgid.sam", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (samse)", $failedjobs);

		## 2014-10-06 : possibly another network writing issue with direct writing to $datadir upon sort. samtools reports truncated bam and spits out gigabytes of text on stderr.
		&PBSWriteCommand("OUT", "$toolsdir/samtools/0.1.19/samtools view -uS /tmp/$project.$rgid.sam | $toolsdir/samtools/0.1.19/samtools sort - /tmp/$project.$root", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (view/sort)", $failedjobs);

		&PBSWriteCommand("OUT", "$toolsdir/samtools/0.1.19/samtools index /tmp/$project.$bam /tmp/$project.$bai", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (index)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$bam $datadir/tmp_files/$sample/$bam", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy BAM)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$bai $datadir/tmp_files/$sample/$bai", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy BAI)", $failedjobs);

		&PBSWriteCommand("OUT", "rm -f /tmp/$project.$rgid.sa* /tmp/$project.$bam /tmp/$project.$bai", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);

		&PBSWriteEnd("OUT");
		close OUT;
		my $samseid = `qsub $script`;
		chomp($samseid);
		while ($samseid !~ m/^\d+\..*/) {
			sleep 5;
			my $tmpfile = `mktemp`;
			chomp($tmpfile);
			my $return = `qstat -x 2>$tmpfile | grep $script`;
			chomp($return);
			my $errorreturn = `cat $tmpfile`;
			if (! ($return or $errorreturn)) {
				$samseid = `qsub $script`;
				chomp($samseid);
			}
			system("rm $tmpfile");

		}

		$samseid =~ s/(\d+)\..*/$1/;
		$samseids .= "$samseid:";
		$mergein .= "INPUT=$datadir/tmp_files/$sample/$bam ";

	}



	## 2. Merge into 1 Primary BAM file
	my $scriptname = "$sample.Make.Primary.BAM";
	my $script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
	my $bamscript = $script;
	my $jobname = "NIPT.".$scriptname;
	open OUT, ">$script";
	#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
	my $depend = ["#PBS -W depend=afterok:".substr($samseids,0,-1)];
	&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 2, 6, $niptqueue, $niptaccount, $depend);

	# the merge commands
	if (scalar(@fqs) > 1) {
		print OUT "## merge using picard.\n";
		&PBSWriteCommand("OUT", "$softwaredir/java/sun-jre-7/bin/java -Xmx5g -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=/tmp/ -jar $toolsdir/picard/picard-tools-1.78/MergeSamFiles.jar VALIDATION_STRINGENCY=SILENT TMP_DIR=/tmp/ USE_THREADING=true SORT_ORDER=coordinate ASSUME_SORTED=false $mergein OUTPUT=/tmp/$project.$sample.$primarybam CREATE_INDEX=true", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (primary BAM merge)", $failedjobs);
	}

	else {
		## sort and index (and dedup) bam.
		&PBSWriteCommand("OUT", "$softwaredir/java/sun-jre-7/bin/java -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=/tmp/ -jar $toolsdir/picard/picard-tools-1.78/SortSam.jar $mergein O=/tmp/$project.$sample.$primarybam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=/tmp/", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (primary BAM sort)", $failedjobs);
		## index.
		&PBSWriteCommand("OUT", "$toolsdir/samtools/0.1.19/samtools index /tmp/$project.$sample.$primarybam /tmp/$project.$sample.$primarybai", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (primary BAM, index)", $failedjobs);
	}

	## Immediately make make final BAM file, just dedup
	# Do dedup
	&PBSWriteCommand("OUT", "$softwaredir/java/sun-jre-7/bin/java -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=/tmp/ -jar $toolsdir/picard/picard-tools-1.78/MarkDuplicates.jar QUIET=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true I=/tmp/$project.$sample.$primarybam O=/tmp/$project.$sample.final_step1.bam M=/tmp/$project.$sample.MarkDuplicates.metrics TMP_DIR=/tmp/", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (final BAM dedup)", $failedjobs);

	# Sort
	&PBSWriteCommand("OUT", "$toolsdir/samtools/0.1.19/samtools sort /tmp/$project.$sample.final_step1.bam /tmp/$project.$sample.final_step2", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (final BAM sort)", $failedjobs);

	# Rename
	&PBSWriteCommand("OUT", "mv /tmp/$project.$sample.final_step2.bam /tmp/$project.$sample.$finalbam", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (final BAM rename)", $failedjobs);

	# Index
	&PBSWriteCommand("OUT", "$toolsdir/samtools/0.1.19/samtools index /tmp/$project.$sample.$finalbam /tmp/$project.$sample.$finalbai", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (final BAM index)", $failedjobs);


	## Copy BAM and BAI files to datadir
	&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.$primarybam $datadir/tmp_files/$sample/$primarybam", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (copy primary BAM)", $failedjobs);
	&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.$primarybai $datadir/tmp_files/$sample/$primarybai", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (copy primary BAI)", $failedjobs);
	&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.$finalbam $datadir/tmp_files/$sample/$finalbam", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (copy final BAM)", $failedjobs);
	&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.$finalbai $datadir/tmp_files/$sample/$finalbai", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (copy final BAI)", $failedjobs);


	# Clean-up
	&PBSWriteCommand("OUT", "rm -f /tmp/$project.$sample.*.bam /tmp/$project.$sample.*.bai /tmp/$project.$sample.*.metrics", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);

	&PBSWriteCommand("OUT", "rm -f $datadir/tmp_files/$sample/$sample*.fastq.bam $datadir/tmp_files/$sample/$sample*.fastq.sai $datadir/tmp_files/$sample/$sample*.fastq.bai", 1);
	&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);


	## Next Step: split by chromosome (idea: launch pre-written PBS job at the end of the previous one)
	print OUT "\n## Submit subsequent jobs.\n";
	for (my $i = 1; $i<= 24; $i++) {
		my $qsubprint = PBSWriteSubmitNextStep($sample, "SplitByChr.$i");
		print OUT $qsubprint;
	}
	#my $qsubprint = PBSWriteSubmitNextStep($sample, "SplitByChr");
	#print OUT $qsubprint;
	&PBSWriteEnd("OUT");
	close OUT;


	for (my $i = 1; $i<= 24; $i++) {
		$scriptname = "$sample.SplitByChr.$i";
		$script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
		$jobname = "NIPT.".$scriptname;
		open OUT, ">$script";
		## 3. Generate per Chr sorted bam : split, dedup, filter & index
		#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
		&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 3, 4, $niptqueue, $niptaccount, undef);
		print OUT "## define chr array for array job assignment\n";
		print OUT "CHRS=(".join(" ",0 .. 22)." X Y)\n";
		print OUT "## split + sort\n";

		## Change chromosome
		print OUT 'CHR=${CHRS['.$i.']}'."\n";
		# split by chr + SortSam&index (to tmp. sortsam needed to keep results identical to old pipeline)
		&PBSWriteCommand("OUT", "$toolsdir/samtools/0.1.19/samtools view -bh $datadir/tmp_files/$sample/$primarybam \$CHR | $softwaredir/java/sun-jre-7/bin/java -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=/tmp/ -jar $toolsdir/picard/picard-tools-1.78/SortSam.jar O=/tmp/$project.$chrfullbam I=/dev/stdin SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=/tmp/ CREATE_INDEX=true", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (Split & Sort)", $failedjobs);

		&PBSWriteCommand("OUT", "$softwaredir/java/sun-jre-7/bin/java -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=/tmp/ -jar $toolsdir/picard/picard-tools-1.78/MarkDuplicates.jar QUIET=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true I=/tmp/$project.$chrfullbam O=/tmp/$project.$chrdedupbam M=/tmp/$project.$sample.\$CHR.MarkDuplicates.metrics TMP_DIR=/tmp/", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (Mark Duplicates)", $failedjobs);

		&PBSWriteCommand("OUT", "$toolsdir/samtools/0.1.19/samtools index /tmp/$project.$chrdedupbam /tmp/$project.$chrdedupbai", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (Index 1)", $failedjobs);

		&PBSWriteCommand("OUT", "($toolsdir/samtools/0.1.19/samtools view -H /tmp/$project.$chrdedupbam ; $toolsdir/samtools/0.1.19/samtools view -L $basedir/Reference_Data/genome.with.niptBlacklist.excluded.bed /tmp/$project.$chrdedupbam | egrep X0:i:1[^0-9] | fgrep X1:i:0 | fgrep XM:i:0 | fgrep XO:i:0 | fgrep XG:i:0 ) | $toolsdir/samtools/0.1.19/samtools view -uS - | $toolsdir/samtools/0.1.19/samtools sort - /tmp/$project.$sample.\$CHR", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (Filter)", $failedjobs);

		&PBSWriteCommand("OUT", "$toolsdir/samtools/0.1.19/samtools index /tmp/$project.$chrbam /tmp/$project.$chrbai", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (Index 2)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$chrbam $datadir/tmp_files/$sample/$chrbam", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy chromosome BAM)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$chrbai $datadir/tmp_files/$sample/$chrbai", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy chromosome BAI)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.\$CHR.MarkDuplicates.metrics $datadir/tmp_files/$sample/$sample.\$CHR.MarkDuplicates.metrics", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy MarkDuplicate metrics)", $failedjobs);

		&PBSWriteCommand("OUT", "rm /tmp/$project.$sample.\$CHR.*", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);

		# Submit next step: Chromosome counting
		print OUT "\n## submit next step\n";
		$qsubprint = PBSWriteSubmitNextStep($sample, "Count.chr$i");
		print OUT $qsubprint;
		&PBSWriteEnd("OUT");

		close OUT;
	}



	## 4. Chromosome counting, count uniquely mapped reads per chromosome
	my $rgid = $Samples{$sample}{'fullid'};
	for (my $i = 1; $i <= 24; $i++) {
		$scriptname = "$sample.Count.chr$i";
		$script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
		$jobname = "NIPT.".$scriptname;
		open OUT, ">$script";
		#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
		&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 1, 4, $niptqueue, $niptaccount, undef);
		print OUT "## define chr array for array job assignment\n";
		print OUT "CHRS=(0 ".join(" ",1 .. 22)." X Y)\n";
		print OUT 'CHR=${CHRS['.$i.']}'."\n";
		print OUT "## get correct split file\n";
		print OUT "SPLIT=`cd $basedir/Reference_Data/gc.split_chr && ls split.chr$i`\n";
		print OUT "## Run Counting \n";

		## Do preliminary clean-up
		&PBSWriteCommand("OUT", "rm -f /tmp/$project.$sample.\$SPLIT.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (prelim clean-up)", $failedjobs);

		## revised to use bedtools.
		&PBSWriteCommand("OUT", "$toolsdir/BedTools/2.17.0/bin/coverageBed -counts -abam $datadir/tmp_files/$sample/$chrbam -b $basedir/Reference_Data/gc.split_chr/\$SPLIT | awk -v rgid='$rgid' -v OFS=',' '{print rgid,\$4,\$5}' | sort -t',' -k2g,2 -k3g > /tmp/$project.$sample.\$SPLIT.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (count)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.\$SPLIT.csv $datadir/tmp_files/$sample/\$SPLIT.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy)", $failedjobs);

		&PBSWriteCommand("OUT", "rm /tmp/$project.$sample.\$SPLIT.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);

		## write out state file.
		print OUT "echo 1 > $datadir/tmp_files/$sample/counting.chr$i.done\n";
		&PBSWriteEnd("OUT");
		close OUT;
	}


	## all scripts are in place => submit the merging step (only at the end to prevent racing conditions)
	my $PrimaryMergeID = `qsub $bamscript`;
	chomp($PrimaryMergeID);
	while ($PrimaryMergeID !~ m/^\d+\..*/) {
		sleep 5;
		my $tmpfile = `mktemp`;
		chomp($tmpfile);
		my $return = `qstat -x 2>$tmpfile | grep $bamscript`;
		chomp($return);
		my $errorreturn = `cat $tmpfile`;
		if (! ($return or $errorreturn)) {
			$PrimaryMergeID = `qsub $bamscript`;
			chomp($PrimaryMergeID);
		}
		system("rm $tmpfile");

	}

	# give pbs-master time to process this round.
	sleep 3;
}



print "Submitted all Samples for mapping + chromosome counting...\n";

my %queue = %Samples;
while (keys(%queue) > 0) {
	# Check for errors and queue
	&CheckErrors();
	#&CheckQueue();

	foreach (keys(%queue)) {
		$sample = $_;
		my $rgid = $Samples{$sample}{'fullid'};

		## if sample is failed, skip processing and delete from queue.
		if ($Samples{$sample}{'failed'}) {
			delete($queue{$sample});
			next;
		}

		my $done = `ls $datadir/tmp_files/$sample/ | grep "counting\.chr.*\.done" | wc -l`;
		chomp($done);
		if ($done < 24) {
			# not ready, goto next sample.
			next;
		}

		# done, submit next steps, delete from queue.
		delete($queue{$sample});

		print " - Sample $sample submitted for merge counting\n";
		## 5. Count Bins: Merge Counting results
		my $scriptname = "$sample.MergeCount";
		my $script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
		my $mergescript = $script;
		my $jobname = "NIPT.".$scriptname;
		open OUT, ">$script";
		#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
		&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 1, 1, $niptqueue, $niptaccount, undef);

		## CHECK ADDED BY GEERT : Number of split files must be equal to 24 (1 per chromosome).
		print OUT "LC=`ls $datadir/tmp_files/$sample/split* | wc -l`\n";
		print OUT "if [[ x\$LC != x24 ]]; then \n";
		print OUT "  echo '$script (missing files)' >> $failedjobs\n";
		print OUT '  exit $?'."\n";
		print OUT "fi\n";
		print OUT "echo 'SAMPLE,CHR,BIN.START,BIN.END,BIN.GC.CONTENT,BIN.N.COUNT,$rgid.COUNT' > /tmp/$project.$sample.binned.50000.raw.counts.filtered.csv\n";

		&PBSWriteCommand("OUT", "cat $datadir/tmp_files/$sample/split* >> /tmp/$project.$sample.binned.50000.raw.counts.filtered.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (cat)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.binned.50000.raw.counts.filtered.csv $datadir/tmp_files/$sample/binned.50000.raw.counts.filtered.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy)", $failedjobs);

		&PBSWriteCommand("OUT", "rm /tmp/$project.$sample.binned.50000.raw.counts.filtered.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);

		print OUT "\n## submit next steps\n";
		$qsubprint = PBSWriteSubmitNextStep($sample, "seqFF");
		print OUT $qsubprint;
		&PBSWriteEnd("OUT");
		close OUT;



		## 8a. Analyze data : seqFF
		$scriptname = "$sample.seqFF";
		$script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
		$jobname = "NIPT.".$scriptname;
		open OUT, ">$script";
		#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
		## More memory needed with BAM file as input (4GB)
		&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 1, 4, $niptqueue, $niptaccount, undef);

		print OUT "# 7. Fetal Fraction SeqFF (new method Kim et al 2015)\n";
		&PBSWriteCommand("OUT", "$softwaredir/R/R-3.2.1/bin/Rscript $basedir/Binaries/SeqFF/SeqFF.R --i=$datadir/tmp_files/$sample/ --f=$sample.final.bam --d=/tmp/ --o=$project.$sample.fetal.fraction.SeqFF --t=bam --r=$basedir/Binaries/SeqFF/SeqFF.RData --b=$basedir/Binaries/SeqFF/SeqFF_bininfo.csv --c=$toolsdir/samtools/0.1.19/samtools", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (SEQFF)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.fetal.fraction.SeqFF $datadir/tmp_files/$sample/$sample.fetal.fraction.SeqFF", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy ff-seqff)", $failedjobs);

		&PBSWriteCommand("OUT", "rm -f /tmp/$project.$sample.fetal.fraction.SeqFF", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);

		## Submit next step: GCC & FF
		print OUT "\n## submit next step\n";
		my $qsubprint = PBSWriteSubmitNextStep($sample, "GCC_FF");
		print OUT $qsubprint;
		&PBSWriteEnd("OUT");
		close OUT;



		## 6. GC Correct Counts & Foetal Fraction
		$scriptname = "$sample.GCC_FF";
		$script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
		$jobname = "NIPT.".$scriptname;
		open OUT, ">$script";
		#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
		&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 1, 2, $niptqueue, $niptaccount, undef);

		&PBSWriteCommand("OUT", "seqFF=\$(grep 'SeqFF' $datadir/tmp_files/$sample/$sample.fetal.fraction.SeqFF | cut -d, -f2)",1);

		&PBSWriteCommand("OUT", "$softwaredir/R/R-3.2.1/bin/Rscript $basedir/Binaries/WorkFlow_Scripts/gcCorrectNIPT.Rscript $datadir/tmp_files/$sample/binned.50000.raw.counts.filtered.csv /tmp/$project.$sample.binned.50000.count.filtered.reads.csv /tmp/$project.$sample.chr.count.filtered.reads.csv /tmp/$project.$sample \$seqFF $referencedatadir{$ref}/", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (GC correct)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.binned.50000.count.filtered.reads.csv $datadir/Results/$sample/$sample.binned.50000.count.filtered.reads.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy binreads)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.chr.count.filtered.reads.csv $datadir/Results/$sample/$sample.chr.count.filtered.reads.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy chrreads)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.autoMedian $datadir/tmp_files/$sample/$sample.autoMedian", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy automedian)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.XMedian $datadir/tmp_files/$sample/$sample.XMedian", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy xmedian)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.YMedian $datadir/tmp_files/$sample/$sample.YMedian", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy ymedian)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.fetal.fraction.X $datadir/tmp_files/$sample/$sample.fetal.fraction.X", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy ffx)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.fetal.fraction.Y $datadir/tmp_files/$sample/$sample.fetal.fraction.Y", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy ffy)", $failedjobs);

		&PBSWriteCommand("OUT", "rm -rf /tmp/$project.$sample.binned.50000.count.filtered.reads.csv /tmp/$project.$sample.chr.count.filtered.reads.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);

		&PBSWriteCommand("OUT", "rm -rf /tmp/$project.$sample.fetal.fraction.? /tmp/$project.$sample.*Median", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up 2)", $failedjobs);

		## submit next steps : WinBin + Y-bin.
		print OUT "\n## submit next step\n";
		$qsubprint = PBSWriteSubmitNextStep($sample, "winBin_Ycounts");
		print OUT $qsubprint;
		&PBSWriteEnd("OUT");
		close OUT;



		## 7. Sliding Window Counts + Ycounts
		$scriptname = "$sample.winBin_Ycounts";
		$script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
		$jobname = "NIPT.".$scriptname;
		open OUT, ">$script";
		#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
		&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 1, 2, $niptqueue, $niptaccount, undef);

		print OUT "## CREATE SLIDING WINDOW COUNTS\n";
		&PBSWriteCommand("OUT", "$softwaredir/R/R-3.2.1/bin/Rscript $basedir/Binaries/WorkFlow_Scripts/getWinbinFixedSize.CMG.R $datadir/Results/$sample/$sample.binned.50000.count.filtered.reads.csv /tmp/$project.$sample.winbin.50000.100.csv 50000 100", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (winbin)", $failedjobs);

		print OUT "## COUNT Y-BINS\n";
		&PBSWriteCommand("OUT", "sed -e 's/\"//g' $datadir/Results/$sample/$sample.binned.50000.count.filtered.reads.csv | grep -f $basedir/Reference_Data/YSpecificBins.csv | cut -d',' -f 1-4,8 > /tmp/$project.$sample.ybins.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (Ycounts)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.winbin.50000.100.csv $datadir/Results/$sample/$sample.winbin.50000.100.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy winbins)", $failedjobs);

		&PBSWriteCommand("OUT", "cp /tmp/$project.$sample.ybins.csv $datadir/tmp_files/$sample/$sample.y_spec_bins.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy Ycounts)", $failedjobs);

		&PBSWriteCommand("OUT", "rm /tmp/$project.$sample.winbin.50000.100.csv /tmp/$project.$sample.ybins.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);

		print OUT "\n## submit next step\n";
		$qsubprint = PBSWriteSubmitNextStep($sample, "Plots");
		print OUT $qsubprint;
		&PBSWriteEnd("OUT");
		close OUT;



		## 8b. Analyze data : Generate Plots
		my $tmpdir = `mktemp -du`; #to results
		chomp($tmpdir);
		my $tmpdir2 = `mktemp -du`; #to tmp_files
		chomp($tmpdir2);
		$scriptname = "$sample.Plots";
		$script = "$datadir/tmp_binaries/$sample/".$scriptname.".sh";
		$jobname = "NIPT.".$scriptname;
		open OUT, ">$script";
		#my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
		&PBSWriteHeader("OUT", $sample, $datadir, $jobname, $admin, 1, 2, $niptqueue, $niptaccount, undef);

		&PBSWriteCommand("OUT", "mkdir $tmpdir $tmpdir2", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (mkdir tempdirs)", $failedjobs);

		&PBSWriteCommand("OUT", "$softwaredir/R/R-3.2.1/bin/Rscript $basedir/Binaries/WorkFlow_Scripts/getReportNIPT.R $tmpdir $referencedatadir{$ref}/referenceNormals.stats $datadir/Results/$sample/$sample.chr.count.filtered.reads.csv $datadir/Results/$sample/$sample.winbin.50000.100.csv $datadir/tmp_files/$sample/$sample.y_spec_bins.csv /$basedir/Reference_Data/diagnosticBand.b37.csv /$basedir/Reference_Data/cytoBand.b37.csv $datadir/tmp_files/$sample/$sample.fetal.fraction.SeqFF", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (Make plots)", $failedjobs);

		&PBSWriteCommand("OUT", 'GESLACHT=$(tail -n +2 '.$tmpdir.'/MedGen-$sample*.summary.csv | cut -d, -f12)', 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (Get gender)", $failedjobs);

		## find Maternal CNVs per sample
		&PBSWriteCommand("OUT", "$softwaredir/R/R-3.2.1/bin/Rscript $basedir/Binaries/WorkFlow_Scripts/Maternal.CNV.Calling.R $datadir/Results/$sample/$sample.binned.50000.count.filtered.reads.csv $sample \$GESLACHT $referencedatadir{$ref} /$tmpdir2/ -4 3 -3 3 8 8 8 8 /$basedir/Reference_Data/diagnosticBand.b37.csv /$basedir/Reference_Data/cytoBand.b37.csv", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (CNVs)", $failedjobs);

                ## find Fetal CNVs per sample
		&PBSWriteCommand("OUT", "$softwaredir/R/R-3.2.1/bin/Rscript $basedir/Binaries/WorkFlow_Scripts/Fetal.CNV.Calling.R $datadir/Results/$sample/$sample.binned.50000.count.filtered.reads.csv $sample \$GESLACHT $referencedatadir{$ref} /$tmpdir2/ -3 3 8 8 8 8", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (Fetal CNVs)", $failedjobs);

		my @fqs = split(/\|/,$Samples{$sample}{'fastqs'});
		&PBSWriteCommand("OUT", "cd $datadir && zcat @fqs | wc -l > /$tmpdir2/$sample.fastqlines.txt",1);
		&PBSWriteCheckFailedCommand("OUT", "$script (zcat lines fastq)", $failedjobs);

		&PBSWriteCommand("OUT", "cp $tmpdir/* $datadir/Results/$sample/", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy report results)", $failedjobs);

		&PBSWriteCommand("OUT", "cp -r $tmpdir2/* $datadir/tmp_files/$sample/", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (copy CNV and fastqlines)", $failedjobs);

		&PBSWriteCommand("OUT", "rm -rf $tmpdir $tmpdir2", 1);
		&PBSWriteCheckFailedCommand("OUT", "$script (clean-up)", $failedjobs);




		# write out $sample/ReadyForReports.txt file.
		print OUT "echo 1 > '$datadir/tmp_files/$sample/ReadyForReports.txt'\n";
		&PBSWriteEnd("OUT");
		close OUT;

		# all scripts in place : submit the Merging step.
		my $MCID = `qsub $mergescript`;
		chomp($MCID);
		while ($MCID !~ m/^\d+\..*/) {
			sleep 5;
			my $tmpfile = `mktemp`;
			chomp($tmpfile);
			my $return = `qstat -x 2>$tmpfile | grep $mergescript`;
			chomp($return);
			my $errorreturn = `cat $tmpfile`;
			if (! ($return or $errorreturn)) {
				$MCID = `qsub $mergescript`;
				chomp($MCID);
			}
			system("rm $tmpfile");

		}
	}
	# take a nap before the next iteration of the queue.
	sleep 30;

}

#############################
## WAIT FOR JOBS TO FINISH ##
#############################
## torque/pbs does not have the -sync options, hence this workaround:
my $finished = 0;
while ($finished == 0) {
	&CheckErrors();
	$finished = 1;
	%queue = %Samples;
	foreach $sample (keys(%queue)) {
		## if sample is failed, we don't have to wait for it.
		if ($Samples{$sample}{'failed'}) {
			delete($queue{$sample});
			next;
		}
		if (!-e "$datadir/tmp_files/$sample/ReadyForReports.txt") {
			$finished = 0;
			sleep 60;
			last;
		}
	}
}

## 9. Prepare Sample report TEX files
# get sample information
#  - sample name (from fastq)
#  - sample gender (estimated)
#  - sample lane (from fastq)
#  - sample barcode (from fastq)
#  - quality (estimated : sd )
# get run info : date, operator, ... (from RunInfo.xml)

&CheckErrors();
my @samplereports;
foreach $sample (keys(%Samples)) {
	## if sample is failed, we don't have to make a report.
	if ($Samples{$sample}{'failed'}) {
		next;
	}

	my $bc = $Samples{$sample}{'bc'};
	my $rgid = $Samples{$sample}{'fullid'};
	my $samplereport = "$datadir/Results/$sample/$rgid.pdf";

	## get qc values for Results/sample/rgID.summary.csv
	my $summary = `tail -n 1 $datadir/Results/$sample/$rgid.summary.csv`;
	chomp($summary);
	my ($sid,$in_run,$qcval,$testquality,$g21,$z21,$g18,$z18,$otherchr,$gother,$zother,$gender,$wbnumber,$gcnumber) = split(/,/,$summary);
	## get nr reads
	my $nrreads = `cut -f 3 -d ',' $datadir/Results/$sample/$sample.chr.count.filtered.reads.csv | tail -n +2 | paste -s -d+ | bc`;
	chomp($nrreads);	
	my $nrreadsY = `awk -F',' '{if(\$2~/Y/) print \$3}' $datadir/Results/$sample/$sample.chr.count.filtered.reads.csv`; ## these are raw reads.
	chomp($nrreadsY);
	my $nrfilteredY = `tail -n +2 $datadir/tmp_files/$sample/$sample.y_spec_bins.csv | cut -f5 -d, | paste -sd+ | bc`;
	chomp($nrfilteredY);
	my $nrfiltYbinsAtLeastOneRead = `tail -n +2 $datadir/tmp_files/$sample/$sample.y_spec_bins.csv | cut -f5 -d, | xargs -I {} echo '{} >= 1' | bc -l | paste -sd+ | bc`;
	chomp($nrfiltYbinsAtLeastOneRead);
	my $ffY = `cat $datadir/tmp_files/$sample/$sample.fetal.fraction.Y`;
	chomp($ffY);
	$ffY = sprintf("%.1f",$ffY*100);
	my $ffX = `cat $datadir/tmp_files/$sample/$sample.fetal.fraction.X`;
	chomp($ffX);
	$ffX = sprintf("%.1f",$ffX*100);
	my $ffSeqFF = `grep 'SeqFF' $datadir/tmp_files/$sample/$sample.fetal.fraction.SeqFF | cut -d, -f2`;
	chomp($ffSeqFF);
	$ffSeqFF = sprintf("%.1f",$ffSeqFF*100);
	# make input file with details
	open OUT, ">$datadir/tmp_files/$sample/$sample.input.tex";
	print OUT '\newcommand{\NettoReadCount}{$\numprint{'.$nrreads.'}$}'."\n";
	print OUT '\newcommand{\Logo}{'."$basedir/Report_Template/CMG.jpg".'}'."\n";
	if ($nrreads < 4000000) {
		print OUT '\newcommand{\nrreads}{\textcolor{red}{'.$nrreads.'}}'."\n";
	}
	else {
		print OUT '\newcommand{\nrreads}{'.$nrreads.'}'."\n";
	}	
	print OUT '\newcommand{\SampleNumber}{'.latexProofString($sample).'}'."\n";
	print OUT '\newcommand{\PipelineVersion}{'.$pversion.'}'."\n";
	if($gender eq "onbepaald"){
		print OUT '\newcommand{\SexFoetus}{\textcolor{red}{'.$gender.'}}'."\n";
	} else {	
		print OUT '\newcommand{\SexFoetus}{'.$gender.'}'."\n";
	}
	print OUT '\newcommand{\RefSet}{'.$referencekit{$ref}.'}'."\n";
	print OUT '\newcommand{\Adaptor}{'.uc($bc).'}'."\n";
	print OUT '\newcommand{\PngEighteen}{'."$datadir/Results/$sample/$rgid".'_chr18.png}'."\n";
	print OUT '\newcommand{\PngOthers}{'."$datadir/Results/$sample/$rgid".'_others.png}'."\n";
	if ($testquality eq 'geslaagd') {
		print OUT '\newcommand{\TestQuality}{'.$testquality.'}'."\n";
	}
	else {
		print OUT '\newcommand{\TestQuality}{\textcolor{red}{'.$testquality.'}}'."\n";
	}
	print OUT '\newcommand{\RunDate}{'.$RunDate.'}'."\n";
	print OUT '\newcommand{\AnalysisDate}{'.$date.'}'."\n";

	if ($qcval <= 1.5) {
		print OUT '\newcommand{\qcvalue}{'.$qcval.'}'."\n";
	}
	else {
		print OUT '\newcommand{\qcvalue}{\textcolor{red}{'.$qcval.'}}'."\n";
	}		
	print OUT '\newcommand{\PngTwentyOne}{'."$datadir/Results/$sample/$rgid".'_chr21.png}'."\n";
	print OUT '\newcommand{\PngThirteen}{'."$datadir/Results/$sample/$rgid".'_chr13.png}'."\n";
	print OUT '\newcommand{\PngAuto}{'."$datadir/Results/$sample/$rgid".'_auto.png}'."\n";
	if ($gender eq "meisje") {
		print OUT '\newcommand{\PngX}{'."$datadir/Results/$sample/$rgid".'_chrX_female.png}'."\n";
		print OUT '\newcommand{\PngY}{'."$datadir/Results/$sample/$rgid".'_chrY_female.png}'."\n";
	}
	elsif ($gender eq "jongen") {
		print OUT '\newcommand{\PngX}{'."$datadir/Results/$sample/$rgid".'_chrX_male.png}'."\n";
		print OUT '\newcommand{\PngY}{'."$datadir/Results/$sample/$rgid".'_chrY_male.png}'."\n";
	}
	## onbepaald geslacht => gebruik male (meest informatief voor XXY in twee cases).
	else {
		print OUT '\newcommand{\PngX}{'."$datadir/Results/$sample/$rgid".'_chrX_male.png}'."\n";
		print OUT '\newcommand{\PngY}{'."$datadir/Results/$sample/$rgid".'_chrY_male.png}'."\n";
	}

	print OUT '\newcommand{\nrreadsY}{'.$nrreadsY.'}'."\n";
	if(($nrfiltYbinsAtLeastOneRead <= 1 && $nrfilteredY < 2.1) || $nrfilteredY >= 10){
		print OUT '\newcommand{\nrfilteredY}{'.$nrfilteredY.'}'."\n";
	} else {
		print OUT '\newcommand{\nrfilteredY}{\textcolor{red}{'.$nrfilteredY.'}}'."\n";		
	}
	print OUT '\newcommand{\ffY}{'.$ffY.'\%}'."\n";
	print OUT '\newcommand{\ffX}{'.$ffX.'\%}'."\n";

	if ($ffSeqFF >= 4) {
	print OUT '\newcommand{\ffSeqFF}{'.$ffSeqFF.'\%}'."\n";
	}
	else {
	print OUT '\newcommand{\ffSeqFF}{\textcolor{red}{'.$ffSeqFF.'\%}}'."\n";
	}

	#######################
	#add CNV table command#
	#######################
	open CNV_IN, "<$datadir/tmp_files/$sample/$sample.CNVs_maternal.txt";
	<CNV_IN>; #header
	my @CNVdoc = <CNV_IN>;
	close CNV_IN;
		
	#remove whitespace lines:
	chomp(@CNVdoc);
	@CNVdoc = grep(/\S/,@CNVdoc);
	
	my $header = 'Staal & Positie & \#Bins & \#Informatieve Bins & Z-Mean & Read Ratio';
	
	#make command:
	print OUT '\\usepackage{longtable}'."\n";
	print OUT '\\definecolor{red}{RGB}{255,0,0}'."\n";
	print OUT '\\newcommand{\\cnvtable}{'."\n";
	print OUT '\\begin{longtable}{clcccc}'."\n";
	print OUT '\\hline'."\n";
	print OUT '\\multicolumn{6}{c}{Maternale CNV\'s}\\\\'."\n";
	print OUT "$header \\\\\n";
	print OUT "\\endhead\n";
	print OUT '\\hline'."\n";
	print OUT '\\multicolumn{6}{l}{\\footnotesize Beslissingsregels \\textbf{Maternale CNV\'s}:}\\\\'."\n";
	print OUT '\\multicolumn{6}{l}{\\footnotesize - Z-Mean $\\geq 3$ voor duplicaties}\\\\'."\n";
	print OUT '\\multicolumn{6}{l}{\\footnotesize - Z-Mean $\\leq -4$ voor deleties op autosomen, $\\leq -3$ op chromosoom X}\\\\'."\n";
	print OUT '\\multicolumn{6}{l}{\\footnotesize - Lengte $\\geq 8$ bins (rood indien aantal informatieve bins $<$ 8)}\\\\'."\n";
	print OUT '\\multicolumn{6}{l}{\\footnotesize   (Lengte $\\geq 6$ bins voor segmenten overlappend met het GJB6 gen)}\\\\'."\n";
	print OUT "\\endfoot\n";

	if(!@CNVdoc){
	        print OUT '\\multicolumn{6}{c}{Geen} \\\\'."\n";
	} else {
		foreach(@CNVdoc){
		        my @line = split(" ");
		        chomp(@line);
		        if($line[4]<8){
		        	print OUT "\\color{red}".latexProofString($line[0])." \& \\color{red}chr$line[1]:$line[2] \& \\color{red}$line[3] \& \\color{red}$line[4] \& \\color{red}$line[5] \& \\color{red}$line[6] \\\\\n";
		        } else {
		        	print OUT latexProofString($line[0])."\& chr$line[1]:$line[2] \& $line[3] \& $line[4] \& $line[5] \& $line[6] \\\\\n";
		        }
		}
	}

	print OUT '\\end{longtable}'."\n";
	print OUT '}'."\n";

	######################
	#add CNV plot command#
	######################
	print OUT '\\newcommand{\\cnvplots}{'."\n";
	if(@CNVdoc){
		print OUT '\\begin{large} \\begin{center} \\bf CNV plots \\end{center} \\end{large}'."\n";
		print OUT '\\noindent'."\n";
		foreach(@CNVdoc){
			my @line = split(" ");
			chomp(@line);
			print OUT '\\vspace{0.5cm}'."\n";
			print OUT '\\makebox[\textwidth][c]{'."\n";
			print OUT "\\includegraphics[width=\\textwidth, height=\\textheight, keepaspectratio]{{$datadir/tmp_files/$sample/$line[0].chr$line[1].$line[2].CNVPlot}.pdf}\n";
			print OUT "}"."\n";
		}
	}
	print OUT '}'."\n";

        ###############################
	# Add FETAL CNV table command #
	###############################
	open FET_CNV_IN, "<$datadir/tmp_files/$sample/fetal_cnvs_$sample/$sample.fetal_CNVs.txt";
	<FET_CNV_IN>; #header
	my @fetCNVdoc = <FET_CNV_IN>;
	close FET_CNV_IN;
		
	#remove whitespace lines:
	chomp(@fetCNVdoc);
	@fetCNVdoc = grep(/\S/,@fetCNVdoc);
	
	my $header = 'Staal & Positie & \#Bins & \#Inform Bins & Z-Mean & Z & Read Ratio';
	
	#make command:
	print OUT '\\newcommand{\\fetalcnvtable}{'."\n";
	print OUT '\\begin{longtable}{clccccc}'."\n";
	print OUT '\\hline'."\n";
	print OUT '\\multicolumn{7}{c}{Foetale CNV\'s}\\\\'."\n";
	print OUT "$header \\\\\n";
	print OUT "\\endhead\n";
	print OUT '\\hline'."\n";
	print OUT '\\multicolumn{7}{l}{\\footnotesize Beslissingsregels \\textbf{Foetale CNV\'s}:}\\\\'."\n";
	print OUT '\\multicolumn{7}{l}{\\footnotesize - Z $\\geq 3$ voor duplicaties}\\\\'."\n";
	print OUT '\\multicolumn{7}{l}{\\footnotesize - Z $\\leq -3$ voor deleties}\\\\'."\n";
	print OUT '\\multicolumn{7}{l}{\\footnotesize - Lengte $\\geq 8$ bins}\\\\'."\n";
	print OUT "\\endfoot\n";

	if(!@fetCNVdoc){
	        print OUT '\\multicolumn{7}{c}{Geen} \\\\'."\n";
	} else {
		foreach(@fetCNVdoc){
		        my @line = split(" ");
		        chomp(@line);
		        print OUT latexProofString($line[0])."\& chr$line[1]:$line[2] \& $line[3] \& $line[4] \& $line[5] \& $line[6] \& $line[7] \\\\\n";
		}
	}

	print OUT '\\end{longtable}'."\n";
	print OUT '}'."\n";

	close OUT;

	## Generate Report (gender specific)
	print("Making samplereport for $sample: ");
	system("cp $basedir/Report_Template/document.tex $datadir/tmp_files/$sample/$sample.tex");
	system('sed -i s$%INPUT_TEX%$'."$datadir/tmp_files/$sample/$sample.input.tex".'$ '."$datadir/tmp_files/$sample/$sample.tex");
	system("cd $datadir/tmp_files/$sample && pdflatex $sample.tex > $sample.latex.log1  && pdflatex $sample.tex > $sample.latex.log2");
	if (-e "$datadir/tmp_files/$sample/$sample.pdf") {
		system("cp $datadir/tmp_files/$sample/$sample.pdf $samplereport");
	}
	else {
		dieWithGrace("Failed to generate final Report", $project);
	}
	push(@samplereports, $samplereport);
	print("done\n");
}


## 10. Final wrap-up
($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$mon++;
$year += 1900;
my $timestamp = "$year.$mon.$day-$hour.$min.$sec";
my $still_failed = [];
my $finished_failed = [];


if (-e $redo_file) {
	## only one sample but this is easiest
	## copy data to original project
	&CheckErrors();
	print "Copying data back to original map\n";
	foreach $sample (keys %Samples) {
		foreach my $section ('Results', 'tmp_binaries', 'tmp_files', 'Job_Output') {
			my $to = "$original_datadir/$section/$sample";
			my $backup = "$original_datadir/$section/$sample.$timestamp";
			my $from = "$datadir/$section/$sample";
			# backup old data in original projects
			if (-d $to) {
				system("mv $to $backup");
			}
			system("mkdir -p $to && cd $from && cp * $to");
			system("cp $runtimeoutput $original_datadir/RunTime.Output_ReDo_$sample.txt");
		}
	}
	system("echo 1 > $status_file");

	my $endtime = Time::Piece->new()->epoch();
        my $now = $endtime - $starttime;
        printf("\n\nAnalysis Finished.\nRunning time:%02d:%02d:%02d\n",int($now/3600),int(($now % 3600)/60),int($now % 60));
	exit;
}

else {
	&CheckErrors();

	# Main run script will loop here until all redo runs are either finished or failed
	%queue = %Samples;
	while (keys(%queue) > 0) {
		foreach $sample (keys %queue) {
			my $redodir = $Samples{$sample}{'failed'}; # can also just be 0
			if ($redodir && -d $redodir) { # double check if dir exists
				if(!-e "$redodir/Status.txt"){
					next;
				}
				my $status = `cat $redodir/Status.txt`;
				chomp($status);
				if($status == 0){
					next;
				}
				if ($status == -1) {
					push(@$still_failed, $sample);
					print "REDO: $sample failed after re-analysis\n";
				}
				if ($status == 1) {
					push(@$finished_failed, $sample);
					print "REDO: $sample succeeded after re-analysis\n";
					$Samples{$sample}{'failed'} = 0;
					system("rm -rf $redodir");
				}
				delete($queue{$sample});
			}
			else {
				delete($queue{$sample});
			}
		}
		sleep(60);
	}

	# Generate Final GetStats report
	&GetStats($txt, $pdf);
	if (! ($txt or $pdf)) {
			dieWithGrace("Failed to generate Run-QC Report", $project);
	}

}

## 11. Backup data to LTS
if (! $debug) {
	print "Copy to $sambadir (SeqPilot server)\n";
	print "Copy reports:\n-------------\n";
	$command = "mkdir -p $ltsreportsdir && cd $datadir && rsync -av @samplereports $txt $pdf $cnvreport $ltsreportsdir";
	print "$command\n";
	system($command);

	# copy CBS plots to Reports
	$command = "mkdir $ltsreportsdir/CNVs && cd $datadir && rsync -av --include '*/' --include '*CBSplot.pdf' --include '*fetal_cnvs_*/*' --exclude '*' $datadir/tmp_files/ $ltsreportsdir/CNVs";
	print "$command\n";
	system($command);

	my @fastqs = split(/\n/, `ls $datadir/*fastq.gz`);
	$command = "mkdir -p $ltsdatadir && cd $datadir && rsync -av @fastqs $checksumfile $checksumresult $checksumerrorresult $runinfo $runparam $ss_file $unusedcountsfile $qcfolder $UBfolder $ltsdatadir";
	print "Backup data:\n------------\n";
	print "$command\n";
	system($command);

	# copy CNVplot data and CBS list to data
        $command = "mkdir $ltsdatadir/CNVs && cd $datadir && rsync -av --include '*/' --include '*CNVData.csv' --include '*CBS.SegmentList.csv' --include '*fetal_cnvs_*/*' --exclude '*' $datadir/tmp_files/ $ltsdatadir/CNVs";
	print "$command\n";
	system($command);
	
	# copy Bin data to data
	$command = "mkdir $ltsdatadir/BinData && cd $datadir && rsync -av $datadir/Results/*/*.binned.50000.count.filtered.reads.csv $ltsdatadir/BinData";
	print "$command\n";
	system($command);

	# check md5sums
	my @failedMD5s;

	open(my $log,">",$checksumLTSReport) || dieWithGrace("Couldn't open $checksumLTSReport",$project);
	print $log "Localfile\tRemotefile\tResult\n";

	open(my $LTSReportmd5,">",$checksumLTSReportfile) || dieWithGrace("Couldn't open $checksumLTSReportfile",$project);
	print $LTSReportmd5 "# Bestand ter controle van integriteit van data\n"; #comments possible in md5 files!
	open(my $LTSDatamd5,">",$checksumLTSDatafile) || dieWithGrace("Couldn't open $checksumLTSDatafile",$project);
	print $LTSDatamd5 "# Bestand ter controle van integriteit van data\n";

	foreach my $file ((@samplereports,$txt,$pdf,$cnvreport)){
		my $remotefile = $ltsreportsdir."/".basename($file);

		my $remotemd5 = md5($remotefile);
		print $LTSReportmd5 "$remotemd5 ".basename($file)."\n";

		print $log "$file\t$remotefile\t";
		if(md5($file) eq $remotemd5){
			print $log "OK\n";
		} else {
			print $log "FAIL\n";
			push(@failedMD5s,"File $file has not been copied correctly to $remotefile");
		}
	}

	foreach my $fullpath (glob("$datadir/tmp_files/*/*CBSplot.pdf")){
		my ($file,$dir,undef) = fileparse($fullpath);
		$dir =~ s/\/$//; # remove trailing /

		my @filepath = splitdir($dir);
        my $remotefile = $ltsreportsdir."/CNVs/".$filepath[-1]."/$file";

		my $remotemd5 = md5($remotefile);
		print $LTSReportmd5 "$remotemd5 "."CNVs/".$filepath[-1]."/$file\n";

		print $log "$fullpath\t$remotefile\t";
		if(md5($fullpath) eq $remotemd5){
			print $log "OK\n";
		} else {
			print $log "FAIL\n";
			push(@failedMD5s,"File $file has not been copied correctly to $remotefile");
		}
	}

	foreach my $file (($checksumfile,$checksumresult,$checksumerrorresult,$unusedcountsfile)){
		my $remotefile = $ltsdatadir."/".basename($file);

		my $remotemd5 = md5($remotefile);
		print $LTSDatamd5 "$remotemd5 ".basename($file)."\n";

		print $log "$file\t$remotefile\t";
		if(md5($file) eq $remotemd5){
			print $log "OK\n";
		} else {
			print $log "FAIL\n";
			push(@failedMD5s,"File $file has not been copied correctly to $remotefile");
		}
	}

	my %databufferMD5s;
	open(my $fh,"<",$checksumfile) || dieWithGrace("Couldn't open $checksumfile",$project);
	while(my $line = <$fh>){
		chomp($line);
		my @el = split(" ",$line);
		my $file = basename($el[1]);
		$databufferMD5s{$file} = $el[0];
	}
	close($fh);

	foreach my $file ((@fastqs,$ss_file,$runinfo,$runparam)){
		my $remotefile = $ltsdatadir."/".basename($file);
		my $localmd5 = md5($file);
		# check if still same as before run
		if(defined $databufferMD5s{basename($file)}){
			if($localmd5 eq $databufferMD5s{basename($file)}){
				print $log "$file\tchecksums.databuffer.md5\tOK\n";				
			} else {
				print $log "$file\tchecksums.databuffer.md5\tFAIL\n";
				push(@failedMD5s,"File $file changed during run");
			}
		} else {
			print $log "$file\tchecksums.databuffer.md5\tNOT_PRESENT\n";
			push(@failedMD5s,"File $file should have been present in checksums from databuffer");
		}

		my $remotemd5 = md5($remotefile);
		print $LTSDatamd5 "$remotemd5 ".basename($file)."\n";

		print $log "$file\t$remotefile\t";
		if($localmd5 eq md5($remotefile)){
			print $log "OK\n";
		} else {
			print $log "FAIL\n";
			push(@failedMD5s,"File $file has not been copied correctly to $remotefile");
		}
	}
	
	foreach my $file (glob("$qcfolder/*")){
		my $remotefile = $ltsdatadir."/".basename($qcfolder)."/".basename($file);
		my $localmd5 = md5($file);

		if(defined $databufferMD5s{basename($file)}){
			if($localmd5 eq $databufferMD5s{basename($file)}){
				print $log "$file\tchecksums.databuffer.md5\tOK\n";				
			} else {
				print $log "$file\tchecksums.databuffer.md5\tFAIL\n";
				push(@failedMD5s,"File $file changed during run");
			}
		} else {
			print $log "$file\tchecksums.databuffer.md5\tNOT_PRESENT\n";
			push(@failedMD5s,"File $file should have been present in checksums from databuffer");
		}

		my $remotemd5 = md5($remotefile);
		print $LTSDatamd5 "$remotemd5 ".basename($qcfolder)."/".basename($file)."\n";

		print $log "$file\t$remotefile\t";
		if($localmd5 eq $remotemd5){
			print $log "OK\n";
		} else {
			print $log "FAIL\n";
			push(@failedMD5s,"File $file has not been copied correctly to $remotefile");
		}
	}

	foreach my $file (glob("$UBfolder/*")){
		my $remotefile = $ltsdatadir."/".basename($UBfolder)."/".basename($file);
		my $localmd5 = md5($file);

		if(defined $databufferMD5s{basename($file)}){
			if($localmd5 eq $databufferMD5s{basename($file)}){
				print $log "$file\tchecksums.databuffer.md5\tOK\n";				
			} else {
				print $log "$file\tchecksums.databuffer.md5\tFAIL\n";
				push(@failedMD5s,"File $file changed during run");
			}
		} else {
			print $log "$file\tchecksums.databuffer.md5\tNOT_PRESENT\n";
			push(@failedMD5s,"File $file should have been present in checksums from databuffer");
		}

		my $remotemd5 = md5($remotefile);
		print $LTSDatamd5 "$remotemd5 ".basename($UBfolder)."/".basename($file)."\n";

		print $log "$file\t$remotefile\t";
		if($localmd5 eq $remotemd5){
			print $log "OK\n";
		} else {
			print $log "FAIL\n";
			push(@failedMD5s,"File $file has not been copied correctly to $remotefile");
		}
	}

	foreach my $fullpath (glob("$datadir/tmp_files/*/*CNVData.csv $datadir/tmp_files/*/*CBS.SegmentList.csv")){
		my ($file,$dir,undef) = fileparse($fullpath);
		$dir =~ s/\/$//; # remove trailing /

		my @filepath = splitdir($dir);
        my $remotefile = $ltsdatadir."/CNVs/".$filepath[-1]."/$file";

		my $remotemd5 = md5($remotefile);
		print $LTSDatamd5 "$remotemd5 "."CNVs/".$filepath[-1]."/$file\n";

		print $log "$fullpath\t$remotefile\t";
		if(md5($fullpath) eq $remotemd5){
			print $log "OK\n";
		} else {
			print $log "FAIL\n";
			push(@failedMD5s,"File $file has not been copied correctly to $remotefile");
		}
	}

	foreach my $file (glob("$datadir/Results/*/*.binned.50000.count.filtered.reads.csv")){
		my $remotefile = $ltsdatadir."/BinData/".basename($file);
		print $log "$file\t$remotefile\t";

		my $remotemd5 = md5($remotefile);
		print $LTSDatamd5 "$remotemd5 "."BinData/".basename($file)."\n";

		if(md5($file) eq $remotemd5){
			print $log "OK\n";
		} else {
			print $log "FAIL\n";
			push(@failedMD5s,"File $file has not been copied correctly to $remotefile");
		}
	}

	close($log);
	close($LTSReportmd5);
	close($LTSDatamd5);

	if(@failedMD5s){
		dieWithGrace(join("\n",@failedMD5s),$project);
	}

	# copy LTS DATA md5 & log  to datadir
	$command = "cd $datadir && rsync -av $checksumLTSReport $checksumLTSDatafile $ltsdatadir";
	print "$command\n";
	system($command);

	# copy LTS REPORT md5 to datadir
	$command = "cd $datadir && rsync -av $checksumLTSReportfile $ltsreportsdir";
	print "$command\n";
	system($command);


}


############
# DATABASE #
############

if($mysql){


#---------
# READ QC:
#---------

open(my $QCFH,"<$txt") || dieWithGrace("Couldn't open $txt: $!", "$project");
for(my $i = 0; $i<2; $i++){
	<$QCFH>; # header
}
while(my $QCline = <$QCFH>){
	chomp($QCline);

	my ($samplename,@rest)
		= split("\t",$QCline);

	@{$Samples{$samplename}}{ ("raw_reads",
				   "filtered_reads",
				   "runperc",
				   "procperc",
				   "sdval",
				   "sdresult",
				   "ffx",
				   "ffy",
				   "seqff",
				   "gender",
				   "result") }
			= @rest[2,3,4,5,1,6,7,8,9,10,23];
	# index is one less because Sample is cut off
}
close($QCFH);

#---------
# QUERIES:
#---------

  #BC
my $bcnamequery = 'SELECT bid FROM `Barcodes` WHERE name = ?';
my $bcnamesth = $dbh->prepare($bcnamequery);

  #Results
my $resultquery = "INSERT INTO `Results` ";
$resultquery   .= "(sid,gender,result) ";
$resultquery   .= "VALUES ";
$resultquery   .= "(?,?,?)";
my $resultsth = $dbh->prepare($resultquery);

  #CNVs
my $cnvquery = "INSERT INTO `Maternal_cnvs` ";
$cnvquery   .= "(sid,chr,begin,end,number_of_bins,
		number_of_informative_bins,Zmean,Read_ratio) ";
$cnvquery   .= "VALUES ";
$cnvquery   .= "(".join(",",("?")x8).")";
my $cnvsth = $dbh->prepare($cnvquery);

  #Segments
my $segquery = "INSERT INTO `Chr_segments` ";
$segquery   .= "(sid,chr,begin,end,number_of_bins,
		number_of_informative_bins,Zmean,Read_ratio) ";
$segquery   .= "VALUES ";
$segquery   .= "(".join(",",("?")x8).")";
my $segsth = $dbh->prepare($segquery);

  #QC
my $qcquery = "INSERT INTO `QC` ";
$qcquery   .= "(sid,sd,sd_result,raw_reads,filtered_reads,run_percentage,processing_percentage) ";
$qcquery   .= "VALUES ";
$qcquery   .= "(".join(",",("?")x7).")";
my $qcsth = $dbh->prepare($qcquery);

  #Foetal fraction
my $ffquery = "INSERT INTO `Foetal_fraction` ";
$ffquery   .= "(sid,ffx,ffy,seqff,seqff_enet,seqff_wrsc,automedian,xmedian,ymedian) ";
$ffquery   .= "VALUES ";
$ffquery   .= "(".join(",",("?")x9).")";
my $ffsth = $dbh->prepare($ffquery);

  #CHR values
my $chrquery = "INSERT INTO `Chr_values`";
$chrquery   .= "(sid,chr,count,gc,gc_auto_gr,z,zz,bm,om,gt)";
$chrquery   .= "VALUES";
$chrquery   .= "(".join(",",("?")x10).")";
my $chrsth = $dbh->prepare($chrquery);

  #Y reads
my $yreadsquery = "INSERT INTO `Yreads`";
$yreadsquery   .= "(sid,raw_reads,filtered_reads,".join(",",map { "bin".$_ } 1..20).") ";
$yreadsquery   .= "VALUES";
$yreadsquery   .= "(".join(",",("?")x23).")";
my $yreadssth = $dbh->prepare($yreadsquery);

  #UB counts
my $UBquery = "INSERT INTO `Unused_counts`";
$UBquery   .= "(pid,bc1,bc2,lane,raw_reads) ";
$UBquery   .= "VALUES";
$UBquery   .= "($pid,".join(",",("?")x4).")";
my $UBsth = $dbh->prepare($UBquery);


#---------
# INSERT:
#---------

print "Inserting data in database\n";
#--------- Unused Counts ---------#

if(-e $unusedcountsfile && !-e $rerun_file){
	print "-Unused Counts\n";
	open(my $UBFH,"<$unusedcountsfile") ||
		dieWithGrace("Couldn't open $unusedcountsfile: $!",$project);
	
	my %UBreads;
	while(my $UBline = <$UBFH>){
	        chomp($UBline);
		my ($UBbc,$UBfile,$UBcount) = split(",",$UBline);
		$UBfile =~ m{_L00(\d)_R1_001};
		my $UBlane = $1;
	
		my (undef,$UBbc1,$UBbc2) = split("_",$UBbc);
	
		# BC 1
		$bcnamesth->execute($UBbc1);
		my ($bc1id) = $bcnamesth->fetchrow_array();
		
		# BC 2
		$bcnamesth->execute($UBbc2);
		my ($bc2id) = $bcnamesth->fetchrow_array();	

		$UBsth -> execute( ($bc1id,$bc2id,$UBlane,$UBcount) );

	}
	close($UBFH);
}
$UBsth		-> finish;
$bcnamesth	-> finish;

print "-Sampledata:\n";

foreach $sample (keys(%Samples)) {
	next unless($sample =~ m{^[0-9]{5}$}); # No blanco, UB or validation samples (6 digits)
	
	if($Samples{$sample}{'failed'}){
		print "\t Skipping failed sample: $sample\n";
		next;
	}

	print "\t -$sample\n";

	my $sid = $Samples{$sample}{"sid"};

	#--------- Results ---------#

	$resultsth -> execute( ($sid,
	                        $Samples{$sample}{"gender"},
        	                $Samples{$sample}{"result"})
			      ); 

	#--------- CNVs ---------#

	my $CNVfile = "$datadir/tmp_files/$sample/$sample.CNVs_maternal.txt";
	open(my $CNVFH,"<$CNVfile") || 
		dieWithGrace("Couldn't open $CNVfile: $!", "$project");
        <$CNVFH>; # header
	while(my $CNVline = <$CNVFH>){
	        chomp($CNVline);
#	        my ($samplename,$chr,$position,$Nbins,$NinfBins,$Zmean,$Readratio) 
		my ($samplename,$chr,$position,@rest) = split(" ",$CNVline);
		my ($begin,$end) = split("-",$position);
                $cnvsth -> execute( ($sid,$chr,$begin,$end,@rest) );
	}
	close($QCFH);

	#--------- Segments ---------#

	my $Segfile = "$datadir/tmp_files/$sample/$sample.CBS.SegmentList.csv";
	open(my $SEGFH,"<$Segfile") || 
		dieWithGrace("Couldn't open $Segfile: $!", "$project");
        <$SEGFH>; # header
	while(my $Segline = <$SEGFH>){
	        chomp($Segline);
#	        my ($samplename,$chr,$begin,$end,$Nbins,$NinfBins,$Zmean,$Readratio)
		my ($samplename,@rest) = split("\t",$Segline);
                $segsth -> execute( ($sid,@rest) );
	}
	close($SEGFH);

	#--------- QC ---------#

	$qcsth -> execute( ($sid,
			    @{$Samples{$sample}}{ (
                                		"sdval",
                                		"sdresult",
						"raw_reads",
						"filtered_reads",
                                		"runperc",
		                                "procperc"
						) }
			    ) );

	#--------- Foetal Fraction ---------#

	my $ffx = `cat $datadir/tmp_files/$sample/$sample.fetal.fraction.X`;
	chomp($ffx);
	my $ffy = `cat $datadir/tmp_files/$sample/$sample.fetal.fraction.Y`;
	chomp($ffy);

	my $automedian = `cat $datadir/tmp_files/$sample/$sample.autoMedian`;
	chomp($automedian);
	my $Xmedian = `cat $datadir/tmp_files/$sample/$sample.XMedian`;
	chomp($Xmedian);
	my $Ymedian = `cat $datadir/tmp_files/$sample/$sample.YMedian`;
	chomp($Ymedian);

	my (undef,$seqff,$enet,$wrsc) = `cat $datadir/tmp_files/$sample/$sample.fetal.fraction.SeqFF`; 
	chomp($seqff);chomp($enet); chomp($wrsc);
	my (undef,$seqff_val) = split(",",$seqff);
	my (undef,$enet_val) = split(",",$enet);
	my (undef,$wrsc_val) = split(",",$wrsc);

	$ffsth -> execute( ($sid,$ffx,$ffy,$seqff_val,
			    $enet_val,$wrsc_val,$automedian,$Xmedian,$Ymedian) );

	#--------- CHR Values ---------#

	# Autosomes
	my $detailsfile = "$datadir/Results/$sample/".$Samples{$sample}{"fullid"}.".details.csv";
	open(my $detFH,"<$detailsfile") ||
			dieWithGrace("Couldn't open $detailsfile: $!",$project);
	<$detFH>; # header
	while(my $detline = <$detFH>){
		chomp($detline);

		#my ($chr, $count, $gc, $gc_auto_gr, $z, $zz, $bm, $om, $gt) = 
		my @chrvalues = (split(",",$detline))[0,2,3,9,34,36,37,38,39];

		$chrsth -> execute( ($sid,@chrvalues) );
	}
	close($detFH);

	# X & Y chrom

	my %XYchrdetails = (
			    "F.X" => {"chr" => "XNC", "col" => 10}, #X not corrected
			    "M.X" => {"chr" => "XC",  "col" => 11},  #X corrected
			    "F.Y" => {"chr" => "YNC", "col" => 12}, #Y not corrected
			    "M.Y" => {"chr" => "YC",  "col" => 13}  #Y corrected
			    );
	foreach my $chrGender (keys %XYchrdetails){
	        $detailsfile = "$datadir/Results/$sample/".$Samples{$sample}{"fullid"}.".details.".$chrGender.".csv";
	        open(my $detFH,"<$detailsfile") ||
	                        dieWithGrace("Couldn't open $detailsfile: $!",$project);
	        <$detFH>; # header
	
	        my $detline = <$detFH>;
	        chomp($detline);

	        #my ($count, $gc, $gc_auto_gr, $z, $zz, $bm, $om, $gt) = 
	        my @chrvalues = (split(",",$detline))[2,3,$XYchrdetails{$chrGender}{"col"},34,35,36,37,38];
	
        	$chrsth -> execute( ($sid,$XYchrdetails{$chrGender}{"chr"},@chrvalues) );
        
	       	close($detFH);
	}


	#--------- Yreads ---------#

	my $yreadsfile = "$datadir/tmp_files/$sample/$sample.y_spec_bins.csv";

	my $nrreadsY = `awk -F',' '{if(\$2~/Y/) print \$3}' $datadir/Results/$sample/$sample.chr.count.filtered.reads.csv`; ## these are raw reads.
	chomp($nrreadsY);
	my $nrfilteredY = `tail -n +2 $yreadsfile | cut -f5 -d, | paste -sd+ | bc`;
	chomp($nrfilteredY);

        open(my $yreadsFH,"<$yreadsfile") ||
		dieWithGrace("Couldn't open $yreadsfile: $!",$project);
        <$yreadsFH>; # header

	my @binYreads;
        while(my $yreadsline = <$yreadsFH>){
		chomp($yreadsline);
		push(@binYreads,(split(",",$yreadsline))[-1]);
	}

        close($yreadsFH);

        $yreadssth -> execute( ($sid,$nrreadsY,$nrfilteredY,@binYreads) );


	#--------- Rerun ---------#

	if($Samples{$sample}{'rerunID'}){
		if(-e $rerun_file){
			#this is the composed sample
			 
			my $rerunquery = "UPDATE `Reruns` SET ";
			$rerunquery   .= "composed_sample=$sid ";
			$rerunquery   .= "WHERE id = ?";
			my $rerunsth = $dbh->prepare($rerunquery);
			$rerunsth -> execute($Samples{$sample}{'rerunID'});
			$rerunsth -> finish();
		}
	}

}

$resultsth 	-> finish;
$cnvsth 	-> finish;
$segsth 	-> finish;
$chrsth		-> finish;
$yreadssth	-> finish;
$ffsth		-> finish;

print "-blanco & unused barcodes:\n";
#BLANCO en UNUSED: (Sample en QC)
foreach $sample (keys(%Samples)) {
	next unless($sample =~ m{blanco|unused_barcodes}i);

        if($Samples{$sample}{'failed'}){
                print "\t Skipping failed sample: $sample\n";
                next;
        }

	print "\t -$sample\n";

	#--------- QC ---------#

        $qcsth -> execute( (@{$Samples{$sample}}{ (
						"sid",
                                                "sdval",
                                                "sdresult",
                                                "raw_reads",
                                                "filtered_reads",
                                                "runperc",
                                                "procperc"
                                                ) }
                            ) );

}
$qcsth 		-> finish;


}

############
## FINISH ##
############
## set status.
system("echo 1 > $status_file");

my $endtime = Time::Piece->new();

if($mysql){
	print "Entering last project data\n";
	my $projendquery = "UPDATE `Projects` SET ";
	$projendquery   .= "date_analysis_stop = ?,";
	$projendquery   .= "number_of_samples = ? ";
	$projendquery   .= "WHERE pid = ?";

	my $projendsth = $dbh->prepare($projendquery);


	my @projendvalues = ($endtime->strftime("%y-%m-%d %H:%M:%S"),
			     $sampleCount,
			     $pid);

	$projendsth -> execute(@projendvalues);
	$projendsth -> finish;

}
print "Disconnecting from database\n";
$dbh->disconnect;

&FinishMail($still_failed,$finished_failed,$endtime->epoch());



##############
## ROUTINES ##
##############
sub md5 {
	my $file = shift;
	open(my $fh,"<",$file) || dieWithGrace("Couldn't open $file to make md5sum",$project);
	my $md5 = Digest::MD5->new->addfile($fh)->hexdigest;
	close($fh);
	return $md5;
}

sub latexProofString {
	my $string = shift;
	my %replace =(  '&'     =>      '\\&',
	                '%'     =>      '\\%',
	                '$'     =>      '\\$',
	                '#'     =>      '\\#',           
	                '_'     =>      '\\_',
	                '{'     =>      '\\{',
	                '}'     =>      '\\}',
	                '~'     =>      '\\textasciitilde{}',
	                '^'     =>      '\\textasciicircum{}',
	                '\\'    =>      '\\textbackslash{}'
	);
	my @badchars =('&','%','\$','#','_','{','}','~','\^','\\\\');
	# badchars need to be escaped to match in regex, but only badchar (not escaped) is 
	# returned in $1 => 2 different lists needed
	my $regex = join("|",@badchars);
	$string =~ s/($regex)/$replace{$1}/g;
	return $string;
}


sub CheckErrors {
	# Fetch failed samples and make separate redo-project (copying all necessary files to new dir)
	if (-e "$failedjobs") {
		# if project is a redo, quit, else might keep looping
		if (-e $redo_file){
			my $message = `cat $failedjobs`;
			$message =~ s/\n/\r\n/g;
			$message = "The following jobs(s) reported to be failed, on a second analysis try:\r\n".$message;
			dieWithGrace($message, $project);
		}
		my @failedjobs = split(/\n/, `cat $failedjobs`);
		chomp(@failedjobs);
		foreach (@failedjobs) {
			$sample = (split(/\//, $_))[-2];
			# already processed
			if ($Samples{$sample}{'failed'}) {
				next;
			}
			my $sampler = $sample.'R';
			my $newdir = $basedatadir."/".$project."_ReDo_".$sample;
			# create new project per sample and copy fastq, runinfo.xml and samplesheet
			system("mkdir $newdir");
			system("cp $datadir/$sample*fastq.gz $runinfo '$newdir/'");
			my $headerlinenr = `grep -nri ',Sample_Project,' $ss_file | cut -d: -f1`;
			chomp($headerlinenr);
			system("head -$headerlinenr $ss_file > $newdir/SampleSheet.csv");
			if($sample =~ m{unused_barcodes}){
				system("echo 'unused_barcodes,,,,XXXX,XXXXXXXX,XXXX,XXXXXXXX,NIPT,$ref' >> $newdir/SampleSheet.csv");
			} else {
				system("grep -i $sample $ss_file | sed \"s/$sampler/$sample/gi\" >> $newdir/SampleSheet.csv");
			}
			system("echo $project > $newdir/Redo.txt");
			#md5 checksums
			my @fastqs = split(/\|/, $Samples{$sample}{'fastqs'});
			foreach my $fq (@fastqs){
				system("grep -i $fq $checksumfile >> $newdir/".basename($checksumfile));
			}
			# append new sample to failed samples list and hash
			system("echo '$sample' >> $datadir/tmp_files/Failed_samples.txt");
			$Samples{$sample}{'failed'} = $newdir;

			# transfer.databuffer.done to start run
			system("echo '1' > $newdir/transfer.databuffer.done");
		}
	}
}

sub CheckQueue {
	# do not include completed jobs
	my $command = "qstat | grep NIPT | grep -v ' C ' | wc -l";
	my $output = `$command`;
	chomp($output);
	while ($output > $max_pbs) {
		# sleep for 5 minutes and check again
		sleep(60 * 5);
		$output = `$command`;
		chomp($output);
	}
}

sub commify {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}

sub dieWithGrace {
	my ($message, $proj) = @_;
	my $file;
	if (-e $datadir) {
		$file = "$datadir/dwg.mail.txt";
	}
	else {
		$file = "/tmp/dwg.mail.txt";
	}
	open OUT, ">$file";
	print OUT "To: $admin\n";
	print OUT "subject: NIPT Analysis Died\n";
	print OUT "from: $sender\n\n";
	print OUT "NIPT analysis died on the following project:\r\n";
	print OUT "$proj\r\n\r\n";
	print OUT "The error message was : \r\n\r\n$message\r\n";
	print OUT "\r\n";
	close OUT;

	## send mail
	system("sendmail -t < $file");
	print  "####################\n";
	print  "# CRITICAL PROBLEM #\n";
	print  "####################\n";
	print  "\n";
	print  "NIPT analysis died on project $proj\n";
	print  "Message was:\n";
	print  "$message\n\n";
	print  "\tMail sent to notify admin of crash.\n";
	print  "Monitor will exit now\n";
	## it died, so there was an error => set project to error state.
	system("echo -1 > $status_file");
	exit;
}

sub PBSWriteSubmitNextStep {
	my ($samplename, $jobsub) = @_;
	my $scriptsub = "$datadir/tmp_binaries/$samplename/$samplename.$jobsub.sh";
	my $output    = "JID=`qsub $scriptsub`\n";
	   $output   .= 'while  [[ ! $JID =~ ^[[:digit:]] ]] ; do'."\n";
	   $output   .= "  sleep 5\n";
	   $output   .= '  tmpFile=$(mktemp)'."\n";
	   $output   .= '  TRASH=`qstat -x 2>$tmpFile | grep '."$scriptsub`\n";
	   $output   .= '  # can be non-zero on grep-not-found and qstat-failure. only resumbit on valid qstat check.'."\n";
	   $output   .= '  if  [[ $? -ne 0 ]] ; then'."\n";
	   $output   .= '     if [ ! -s $tmpFile ]; then'."\n";
	   $output   .= "         JID=`qsub $scriptsub`\n";
	   $output   .= "     fi\n";
	   $output   .= "  fi\n";
	   $output   .= '  rm $tmpFile'."\n";
	   $output   .= "done\n";
	return($output);
}

sub FinishMail {
	my $redo_fail = shift();
	my $redo_succes = shift();
	my $endtime = shift();
	## finished (data still needs to be copied)
	my $now = $endtime - $starttime;
	printf("\n\nAnalysis Finished.\nRunning time:%02d:%02d:%02d\n",int($now/3600),int(($now % 3600)/60),int($now % 60));

	## send mail to admin/supervisor.
	my $mail = "$datadir/tmp_files/finished.mail.final.txt";
	open OUT, ">$mail";
	print OUT "To: $addressed\n";
	print OUT "CC: $admin\n";

	if (@$redo_fail || @$redo_succes) {
		print OUT "subject: NIPT analysis finished (failed samples re-analysed): $project\n";
	}
	else {
		print OUT "subject: NIPT analysis finished: $project\n";
	}
	print OUT "from: $sender\n";
	## To interpret HTML in email
	print OUT 'Content-Type: text/html; charset="UTF-8"'."\n";
	print OUT 'MIME-Version: 1.0'."\n";
	print OUT 'Content-Disposition: inline'."\n";
	my $mailcontent;
	{
		local $/;
		open my $mailfh, '<', "$basedir/Mail_Template/NIPTmail.html" 
			|| dieWithGrace("can't open $basedir/Mail_Template/NIPTmail.html: $!",$project);
		$mailcontent = <$mailfh>;
		close $mailfh;
	}
	my $failedprint = "";
	if (@$redo_succes) {
		$failedprint = "The following samples succeeded on re-analysis:<br>\n";
		foreach $sample (@$redo_succes){
			$failedprint .= " - $sample<br>\n";
		}
	}
	if (@$redo_fail) {
		$failedprint .= "<br>" if(@$redo_succes);
		$failedprint .= "The following samples also failed on re-analysis:<br>\n";
		foreach $sample (@$redo_fail) {
			$failedprint .= " - $sample<br>\n";
		}
		$failedprint .= "The NIPT administrator should look into this issue!!<br>\n";
	}
	
	$mailcontent =~ s/#PROJECT#/$project/;
	$mailcontent =~ s*#REPORTS#*fsnf/niptdata2/NIPT/reports/$fullproject*;
	$mailcontent =~ s*#REPORTLOC#*\\\\fsnf\\niptdata2\\NIPT\\reports\\$fullproject*;
	$mailcontent =~ s*#DATA#*fsnf/niptdata2/NIPT/data/$fullproject*;
	$mailcontent =~ s*#DATALOC#*\\\\fsnf\\niptdata2\\NIPT\\data\\$fullproject*;
	$mailcontent =~ s*#WARNINGS#*$failedprint*;
	if($failedprint){
			$mailcontent =~ s*#DOWARN#*block*;
	} else {
			$mailcontent =~ s*#DOWARN#*none*;               
	}
	my $runtime = sprintf("%02d:%02d:%02d<br>\n",int($now/3600),int(($now % 3600)/60),int($now % 60));
	$mailcontent =~ s*#RUNTIME#*$runtime*;
	print OUT $mailcontent;

	close OUT;

	## send mail
	system("sendmail -t < '$mail'");
	#finished
}

sub GetStats {
	## 25-12-2015: Subroutine replaces GetStats.pl
	my ($txtreport, $pdfreport) = @_;

	my $totalrun = 0;

	my %stats;
	my $title = $rundate;
	$title =~ s/_/\\_/g;

	foreach $sample(keys %Samples) {
		## failed samples : skip these
		if ($Samples{$sample}{'failed'}) {
			next;
		}

		my $rgid = $Samples{$sample}{'fullid'};
		print "$sample ($rgid)\n";
		$stats{$sample}{'bc'} = $Samples{$sample}{'bc'};

		#raw fastq counts
		my $count = `cat $datadir/tmp_files/$sample/$sample.fastqlines.txt`;

		chomp($count);
		$count = $count / 4;
		$totalrun += $count;
		$stats{$sample}{'raw'} = $count;

		## get qc values for Results/sample/rgID.summary.csv and Results/sample/rgID.details.csv
		my $details = `cat $datadir/Results/$sample/$rgid.details.csv`;
		my @lines = split(/\n/, $details);
		foreach (@lines) {
			chomp($_);
			my @columns = split(/,/, $_);
			if ($columns[0] == '18') {
				$stats{$sample}{'z18'} = $columns[34];
				$stats{$sample}{'zz18'} = $columns[36];
				$stats{$sample}{'bm18'} = $columns[37];
				$stats{$sample}{'om18'} = $columns[38];
			}
			elsif ($columns[0] == '21') {
				$stats{$sample}{'z21'} = $columns[34];
				$stats{$sample}{'zz21'} = $columns[36];
				$stats{$sample}{'bm21'} = $columns[37];
				$stats{$sample}{'om21'} = $columns[38];
			}
			elsif ($columns[0] == '13') {
				# GT 13 not in summary file, quickest way to get this
				$stats{$sample}{'g13'} = $columns[39];
				$stats{$sample}{'z13'} = $columns[34];
				$stats{$sample}{'zz13'} = $columns[36];
				$stats{$sample}{'bm13'} = $columns[37];
				$stats{$sample}{'om13'} = $columns[38];
			}
		}
		my $summary = `tail -n 1 $datadir/Results/$sample/$rgid.summary.csv`;
		chomp($summary);
		my ($sid,$in_run,$qcval,$testquality,$g21,$z21,$g18,$z18,$otherchr,$gother,$zother,$gender,$wbnumber,$gcnumber) = split(/,/,$summary);
		my $g13 = $stats{$sample}{'g13'};
		$stats{$sample}{'qcval'} = $qcval;
		$stats{$sample}{'gender'} = $gender;
		## more detailed, cfr meeting leuven 2010-10-07, is set in the getReportNIPT.R function.
		$testquality =~ s/ /_/g;
		$stats{$sample}{'qcres'} = $testquality;

		my $res = '';
		if ($g18 eq 'trisomie') {
			$res .= 'T18,';
		}
		elsif ($g18 eq 'monosomie') {
			$res .= 'M18,';
		}
		elsif ($g18 ne 'normaal') {
			$res .= "$g18 18,";
		}
		if ($g21 eq 'trisomie') {
			$res .= 'T21,';
		}
		elsif ($g21 eq 'monosomie') {
			$res .= 'M21,';
		}
		elsif ($g21 ne 'normaal') {  ## e.g. "Geen monosomie"
			$res .= "$g21 21,";
		}
		if ($g13 eq 'trisomie') {
			$res .= 'T13,';
		}
		elsif ($g13 eq 'monosomie') {
			$res .= 'M13,';
		}
		elsif ($g13 ne 'normaal') {  ## e.g. "Geen monosomie"
			$res .= "$g13 13,";
		}
		if ($res ne '') {
			$res =~ s/ /_/g;
			$res = substr($res,0,-1); # remove trailing comma
		}
		if ($testquality ne 'geslaagd') {
			if ($res ne '') {
				$res .= ",(FA)";
			} else {
				$res = "FA";
			}
		}
		if ($res eq '') {
			$res = "NORM";
		}

		$stats{$sample}{'aberr'} = $res;

		## get nr reads
		my $nrreads = `cut -f 3 -d ',' $datadir/Results/$sample/$sample.chr.count.filtered.reads.csv | tail -n +2 | paste -s -d+ | bc`;
		chomp($nrreads);
		$stats{$sample}{'final'} = $nrreads;
		## get nr raw Y reads
		my $nrreadsY = `awk -F',' '{if(\$2~/Y/) print \$3}' $datadir/Results/$sample/$sample.chr.count.filtered.reads.csv`;
		chomp($nrreadsY);
		$stats{$sample}{'rawY'} = $nrreadsY;
		## fetal fraction X
		my $ffx = `cat $datadir/tmp_files/$sample/$sample.fetal.fraction.X`;
		chomp($ffx);
		$stats{$sample}{'ffx'} = $ffx;
		## fetal fraction Y
		my $ffy = `cat $datadir/tmp_files/$sample/$sample.fetal.fraction.Y`;
		chomp($ffy);
		$stats{$sample}{'ffy'} = $ffy;
		## fetal fraction SeqFF
		my $ffSeqFF = `grep 'SeqFF' $datadir/tmp_files/$sample/$sample.fetal.fraction.SeqFF | cut -d, -f2`;
		chomp($ffSeqFF);
		$stats{$sample}{'ffseqff'} = $ffSeqFF;
	}

	open OUT, ">".$txtreport;
	print OUT "Rundate: $rundate\n";
	print OUT "Sample\tBarcode\tSD.Value\tFASTQ.Reads\tReportReads\tRun.Percentage\tProcessing.Percentage\tSD.Result\tFoetal_Fraction_chrX\tFoetal_Fraction_chrY\tFoetal_Fraction_SeqFF\tGender\tZ_chr21\tZZ_chr21\tBM_chr21\tOM_chr21\tZ_chr18\tZZ_chr18\tBM_chr18\tOM_chr18\tZ_chr13\tZZ_chr13\tBM_chr13\tOM_chr13\tResult\n";

	my $unusedTreshold = 144000;
	my $plotUnused = 0;
	my $blancoTreshold = 100000;

	my @samples = sort keys(%stats);
	my $textable = '';
	foreach (@samples) {
		my $s = $_;
		my $bc = $stats{$s}{'bc'};
		my $sd = $stats{$s}{'qcval'};
		my $sdres = $stats{$s}{'qcres'};
		my $sex = $stats{$s}{'gender'};
		my $res = $stats{$s}{'aberr'};
		my $raw = $stats{$s}{'raw'};
		#my $rawY = $stats{$s}{'rawY'};
		my $ffX = sprintf("%.1f",($stats{$s}{'ffx'}*100));
		my $ffY = sprintf("%.1f",($stats{$s}{'ffy'}*100));
		my $ffSeqFF = sprintf("%.1f",($stats{$s}{'ffseqff'}*100));
		my $report = $stats{$s}{'final'};
		my $runperc =  sprintf("%.1f",($raw/$totalrun*100 ));
		my $procperc = sprintf("%.1f",($report/$raw*100 ));
		my $z18 = sprintf("%.2f",$stats{$s}{'z18'});
		my $zz18 = sprintf("%.2f",$stats{$s}{'zz18'});
		my $bm18 = sprintf("%.2f",$stats{$s}{'bm18'});
		my $om18 = sprintf("%.2f",$stats{$s}{'om18'});
		my $z21 = sprintf("%.2f",$stats{$s}{'z21'});
		my $zz21 = sprintf("%.2f",$stats{$s}{'zz21'});
		my $bm21 = sprintf("%.2f",$stats{$s}{'bm21'});
		my $om21 = sprintf("%.2f",$stats{$s}{'om21'});
		my $z13 = sprintf("%.2f",$stats{$s}{'z13'});
		my $zz13 = sprintf("%.2f",$stats{$s}{'zz13'});
		my $bm13 = sprintf("%.2f",$stats{$s}{'bm13'});
		my $om13 = sprintf("%.2f",$stats{$s}{'om13'});
		print OUT "$s\t$bc\t$sd\t$raw\t$report\t$runperc\t$procperc\t$sdres\t$ffX\t$ffY\t$ffSeqFF\t$sex\t$z21\t$zz21\t$bm21\t$om21\t$z18\t$zz18\t$bm18\t$om18\t$z13\t$zz13\t$bm13\t$om13\t$res\n";
		if($s =~ m/unused/i){
			$s = "Unused Barcodes (UB)";
			if($report >= $unusedTreshold){
				$textable .= "\\rowfont{\\color{red}}";
				$plotUnused = 1;
			}
		}
		if($s =~ m/blanco/i){
			$s = ucfirst($s);
			if($report >= $blancoTreshold){
				$textable .= "\\rowfont{\\color{red}}";
			}
		}
		$textable .= latexProofString($s)." & $bc & $sd & ".commify($raw)." & ".commify($report)." & ".($runperc)." & ".($procperc). '\\\\'."\n";
	}
	close OUT;

	open OUT, ">/tmp/nipt.report.$rundate.R";
	print OUT "data <- read.table('$txtreport',skip=1,as.is=TRUE,header=TRUE)\n";
	print OUT "data[data[,1]=='unused_barcodes',1]<-'UB'\n";
	print OUT "data[data[,1]=='blanco',1]<-'Blanco'\n";
	print OUT "runperc <- data[,6]\n";
	print OUT "labels_run <- data[,1]\n";
	print OUT "procperc <- data[,7]\n";
	print OUT "max <- max(runperc) + 10\n if (max > 100) {max <- 100}\n";

	print OUT "options(bitmapType='cairo')\n";

	print OUT "png(file='/tmp/$rundate"."_pool.png',width=240,height=80,units='mm',res=300)\n";
	print OUT 'par(mar=c(3,4.1,0.5,2.1))'."\n";
	print OUT "values <- paste(runperc,'%',sep='')\n";
	print OUT "mp<-barplot(runperc,ylab='Percentage',ylim=c(0,max))\n";
	print OUT "text(mp,runperc +1.8,values,xpd = TRUE,srt=90)\n";
	print OUT "text(mp,par('usr')[3]-0.5,labels=labels_run,srt=45,adj=1,xpd=TRUE)\n";
	print OUT "invisible(dev.off())\n";

	print OUT "png(file='/tmp/$rundate"."_loss.png',width=240,height=80,units='mm',res=300)\n";
	print OUT 'par(mar=c(3,4.1,0.7,2.1))'."\n";
	#print OUT "values <- paste(procperc,'%',sep='')\n";
	print OUT "mp<-barplot(procperc,ylab='Percentage',ylim=c(0,100))\n";
	print OUT "text(mp,par('usr')[3]-1,labels=data[,1],srt=45,adj=1,xpd=TRUE)\n";
	print OUT "invisible(dev.off())\n";
	
	if($plotUnused){
		print OUT "png(file='/tmp/$rundate"."_unused.png',width=480,height=160,units='mm',res=300)\n";
		print OUT 'x<-aggregate(count~name,read.table("'.$unusedcountsfile.'",header=F,sep=",",stringsAsFactors = F,col.names=c("name","lane","count")),sum)'."\n";
		print OUT 'x.sort<-x[order(x$count,decreasing=T),]'."\n";
		print OUT 'x.sort$name <- sapply(x.sort$name,gsub,pattern="unused_",replacement="")'."\n";
		print OUT 'q<-quantile(x.sort$count,names=F,na.rm=T)'."\n";
		print OUT 'x.outliers<-rbind(x.sort[x.sort$count > (q[4]+q[4]-q[2]),],data.frame(name="Rest",count=sum(x.sort[x.sort$count <= (q[4]+q[4]-q[2]),"count"])))'."\n";
		print OUT 'par(mar=c(5,4.1,2.5,2.1))'."\n";
		print OUT 'mp<-barplot(x.outliers$count,ylab="Raw FASTQ Reads")'."\n";
		print OUT 'text(mp,x.outliers$count,adj=c(0.5,-0.3),sapply(x.outliers$count,USE.NAMES = F,FUN = function(x,som) paste(x," (",round(100*x/som,1),"%)",sep=""),som=sum(x.outliers$count)),xpd = TRUE,srt=0)'."\n";
		print OUT 'text(mp,par("usr")[3],labels=x.outliers$name,srt=45,adj=1.2,xpd=TRUE)'."\n";
		print OUT 'invisible(dev.off())'."\n";
	}
	close OUT;
	print "Generating global QC plots\n";
	system("$softwaredir/R/R-3.2.1/bin/Rscript /tmp/nipt.report.$rundate.R");
	system("$softwaredir/R/R-3.2.1/bin/Rscript $basedir/Binaries/WorkFlow_Scripts/GRXYPlot.R $datadir/ $referencedatadir{$ref}/ $datadir/tmp_files/$rundate");

	open OUT, ">/tmp/$rundate.report.tex";
	print OUT '\documentclass[a4paper,10pt]{article}'."\n";
	print OUT '\usepackage[left=2cm,top=1.5cm,right=1.5cm,bottom=2.5cm,nohead]{geometry}'."\n";
	print OUT '\usepackage{longtable}'."\n";
	print OUT '\usepackage[T1]{fontenc}'."\n";
	print OUT '\usepackage{fancyhdr}'."\n";
	print OUT '\usepackage[latin9]{inputenc}'."\n";
	print OUT '\usepackage{color}'."\n";
	print OUT '\usepackage[pdftex]{graphicx}'."\n";
	print OUT '\usepackage{tabu}'."\n"; # to color row
	print OUT '\definecolor{grey}{RGB}{160,160,160}'."\n";
	print OUT '\definecolor{darkgrey}{RGB}{100,100,100}'."\n";
	print OUT '\definecolor{red}{RGB}{255,0,0}'."\n";
	print OUT '\definecolor{orange}{RGB}{238,118,0}'."\n";
	print OUT '\setlength\LTleft{0pt}'."\n";
	print OUT '\setlength\LTright{0pt}'."\n";
	print OUT '\begin{document}'."\n";
	print OUT '\pagestyle{fancy}'."\n";
	print OUT '\fancyhead{}'."\n";
	print OUT '\renewcommand{\footrulewidth}{0.4pt}'."\n";
	print OUT '\renewcommand{\headrulewidth}{0pt}'."\n";
	print OUT '\fancyfoot[R]{\today\hspace{2cm}\thepage\ of \pageref{endofdoc}}'."\n";
	print OUT '\fancyfoot[C]{}'."\n";
	print OUT '\fancyfoot[L]{QC Report for Run ``'.$title.'"}'."\n";
	print OUT '\let\oldsubsubsection=\subsubsection'."\n";
	print OUT '\renewcommand{\subsubsection}{%'."\n";
	print OUT '  \filbreak'."\n";
	print OUT '  \oldsubsubsection'."\n";
	print OUT '}'."\n";
	#table
	print OUT '\subsection*{Run '.substr($rundate,0,2).'-'.substr($rundate,2,2).'-'.substr($rundate,-2).': Quality Details}'."\n";
	print OUT '\hspace{0.2cm}{\scriptsize\begin{tabu}{lllllll}'."\n";
	print OUT '\hline'."\n";
	print OUT '\textbf{Sample} & \textbf{Barcode} & \textbf{SD} & \textbf{Raw.Reads} & \textbf{Final.Reads} & \textbf{Pool.\%} & \textbf{Final.\%} \\\\'."\n";
	print OUT '\hline'."\n";
	print OUT $textable ;
	print OUT '\hline'."\n";
	print OUT '\end{tabu}}'."\n";
	print OUT '\subsection*{Raw Read Distribution By Sample}'."\n";
	print OUT '\includegraphics[width=\textwidth,keepaspectratio=true]{/tmp/'.$rundate.'_pool.png}'."\n";

	print OUT '\subsection*{Reads Retained After Processing}'."\n";
	print OUT '\includegraphics[width=\textwidth,keepaspectratio=true]{/tmp/'.$rundate.'_loss.png}'."\n";

	if($plotUnused){
		print OUT '\subsection*{Raw FASTQ Reads for Unused Barcodes}'."\n";
		print OUT '\includegraphics[width=\textwidth,keepaspectratio=true]{/tmp/'.$rundate.'_unused.png}'."\n";
	}

	print OUT '\subsection*{Genomic Representations of Chromosomes X and Y}'."\n";
	print OUT '\includegraphics[width=\textwidth,keepaspectratio=true]{'.$datadir.'/tmp_files/'.$rundate.'_M_GRXvsSeqFF.pdf}'."\n";
	print OUT '\includegraphics[width=\textwidth,keepaspectratio=true]{'.$datadir.'/tmp_files/'.$rundate.'_F_GRXvsSeqFF.pdf}'."\n";
	print OUT '\includegraphics[width=\textwidth,keepaspectratio=true]{'.$datadir.'/tmp_files/'.$rundate.'_GRYvsSeqFF.pdf}'."\n";
	print OUT '\includegraphics[width=\textwidth,keepaspectratio=true]{'.$datadir.'/tmp_files/'.$rundate.'_GRYvsGRX.pdf}'."\n";

	print OUT '\label{endofdoc}'."\n";
	print OUT '\end{document}'."\n";

	close OUT;
	print "Generating QC.Report.pdf\n";
	system("cd /tmp && (pdflatex $rundate.report.tex && pdflatex $rundate.report.tex) >/dev/null 2>&1 ");
	system("cp /tmp/$rundate.report.pdf $pdfreport");
	system("cp /tmp/$rundate.report.tex $datadir/tmp_files/");
	print "Cleaning up\n";
	system("rm -f /tmp/*$rundate*");

	# Make CNVs.txt file (probably integrate this later in sample report)
	my $header;
	system("touch $cnvreport");
	foreach $sample (@samples) {
		my $samplecnvreport = "$datadir/tmp_files/$sample/$sample.CNVs_maternal.txt";
		if (! $header) {
			$header = `head -1 $samplecnvreport`;
			chomp($header);
			system("echo $header >> $cnvreport");
		}
		system("tail -n +2 $samplecnvreport >> $cnvreport") if($sample !~ m/(?:unused|blanco)/i);

	}

}


sub PBSWriteCheckFailedCommand {
	my ($filehandle, $failed, $failedfile) = @_;
	if (tell($filehandle) == -1) {
		print "Outputfile '$filehandle' is not open for writing! Exit script\n";
		exit;
	}
	print $filehandle 'if [ "$?" -ne "0" ] ; then'."\n";
	print $filehandle "  echo '$failed' >> $failedfile\n";
	print $filehandle '  exit $?'."\n";
	print $filehandle "fi\n\n";

}


sub PBSWriteCommand {
	my ($filehandle, $run_command, $echo) = @_;
	if (tell($filehandle) == -1) {
		print "Outputfile '$filehandle' is not open for writing! Exit script\n";
		exit;
	}
	my $echo_command = $run_command;
	print $filehandle "echo 'Command:'\n";
	print $filehandle "echo '========'\n";
	if ($echo == 1) {
		#$echo_command =~ s/'/\\'/g;
		$echo_command =~ s/"/\\"/g;
		print $filehandle "echo \"$echo_command\"\n";
	}
	print $filehandle "$run_command\n";
}


sub PBSWriteEnd {
	my $filehandle = shift;
	if (tell( $filehandle ) == -1) {
		print "Outputfile '$filehandle' is not open for writing! Exit script\n";
		exit;
	}
	print $filehandle "\n\necho 'End Time : ' `date`\n";
	print $filehandle "printf 'Execution Time = \%dh:\%dm:\%ds\\n' \$((\$SECONDS/3600)) \$((\$SECONDS%3600/60)) \$((\$SECONDS%60))\n";

}


sub PBSWriteHeader {
	my ($filehandle, $samplename, $dir, $jobname, $tomail, $cpu, $mem, $queue, $hpc_account, $additional) = @_;
	if (tell($filehandle) == -1) {
		print "Outputfile '$filehandle' is not open for writing! Exit script\n";
		exit;
	}
	print $filehandle "#!/usr/bin/env bash\n";
	print $filehandle "#PBS -m a\n";
	print $filehandle "#PBS -M $tomail\n";
	print $filehandle "#PBS -d $dir\n";
	print $filehandle "#PBS -l nodes=1:ppn=$cpu,mem=$mem"."g\n";
	print $filehandle "#PBS -N $jobname\n";
	print $filehandle "#PBS -o $dir/Job_Output/$samplename/$jobname.o.txt.\$PBS_JOBID\n";
	print $filehandle "#PBS -e $dir/Job_Output/$samplename/$jobname.e.txt.\$PBS_JOBID\n";
	print $filehandle "#PBS -V\n";
	print $filehandle "#PBS -q $queue\n";
	if ($hpc_account) {
		print $filehandle "#PBS -A $hpc_account\n";
	}
	## additional PBS directives?
	if (ref($additional) eq 'ARRAY') {
		foreach my $directive (@$additional) {
			print $filehandle "$directive\n";
		}
	}
	print $filehandle "\necho 'Running on : ' `hostname`\n";
	print $filehandle "echo 'Start Time : ' `date`\n";
}

