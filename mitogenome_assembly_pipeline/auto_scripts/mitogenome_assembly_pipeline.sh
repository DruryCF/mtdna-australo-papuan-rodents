#!/usr/bin/bash

#############################################################################
#
# mitogenome assembly pipeline
#
# C.F.Drury
#
#############################################################################

echo -e "script executed\n"

## set argument defaults

# arguments requiring option
METHOD="all"
RAW_READS="${PWD}/0_raw"
MAX_MEMORY="15"
OUTPREFIX=""
OUTDIR=${PWD}
SAMPLE_IDS=""
LOG_PATTERN=""
LOG_PATTERN_EXCLUDE=NULL
METHODS_LIST=""
UDOCKER_DIRECTORY=""
MIRA_THREADS=0
mitobim_options=( quick 2step )
kmer_values=(33 21)
concat=( no yes )
ref_guide=( no yes )

# boolean arguments
PREP_READS="no"
DEDUP_READS="yes"
TRIM_READS="yes"
TRIM_ADAPTORS="yes"
OVERRIDE_RAW="no"
NFS_STOP="yes"
UDOCKER="no"
ASSEMBLY_STATS="no"
SKIP_FILE_COMPILATION="no"
RESUME_MITOBIM="no"
TRIM_ADAPTORS="yes"

HELP="no"

# set number of parallel processes for assembly stats
num_processes=50

## extract options and arguments into variables

TEMP=`getopt -o m:r:pdtah --long method:,raw_reads:,prep_reads,dedup_reads,trim_reads,trim_adaptors,remove_dedup,max_memory:,mt_reference:,mt_annotated_ref:,outprefix:,outdir:,sample_id_list:,stop_mira_nfs,udocker,udocker_dir:,assembly_stats,log_pattern:,log_pattern_exclude:,exclude_genes_read_depth:,skip_file_compilation,stats_methods_list:,expected_genes:,resume_mitobim,mitobim_opts:,novoplasty_opts_kmer:,novoplasty_opts_concat:,novoplasty_opts_guide:,help -- "$@"`
eval set -- "$TEMP"

# create empty variable for storing list of supplied arguments
supplied_args=( )

while true ; do
    case "$1" in
		-m|--method)
            case "$2" in
                *) METHOD=$2 ; supplied_args+=( "supplied --method: ${METHOD}\n" ) ; shift 2 ;;
            esac ;;
		-r|--raw_reads)
            case "$2" in
                *) RAW_READS=$2 ; supplied_args+=( "supplied --raw_reads: ${RAW_READS}\n" ) ; shift 2 ;;
            esac ;;
		-p|--prep_reads) PREP_READS="yes" ; supplied_args+=( "supplied --prep_reads: ${PREP_READS}\n" ); shift ;;
		-d|--dedup_reads) DEDUP_READS="no" ; if [ "${PREP_READS}" == "yes" ] ; then supplied_args+=( "supplied --dedup_reads: ${DEDUP_READS}\n" ); fi ; shift ;;
		-t|--trim_reads) TRIM_READS="no" ; if [ "${PREP_READS}" == "yes" ] ; then supplied_args+=( "supplied --trim_reads: ${TRIM_READS}\n" ); fi ; shift ;;
		-a|--trim_adaptors) TRIM_ADAPTORS="no" ; if [ "${PREP_READS}" == "yes" ] ; then supplied_args+=( "supplied --trim_adaptors: ${TRIM_ADAPTORS}\n" ); fi ; shift ;;
		--remove_dedup) REMOVE_DEDUP="yes"; if [ "${PREP_READS}" == "yes" ] ; then supplied_args+=( "supplied --remove_dedup: ${REMOVE_DEDUP}\n" ); fi ; shift ;;
		--max_memory)
            case "$2" in
                *) MAX_MEMORY=$2; supplied_args+=( "supplied --max_memory: ${MAX_MEMORY}\n" ) ; shift 2 ;;
            esac ;;
		--mt_annotated_ref)
            case "$2" in
                *) MT_REFERENCE_ANNOTATED=$2 ; supplied_args+=( "supplied --mt_annotated_ref: ${MT_REFERENCE_ANNOTATED}\n" ) ; shift 2 ;;
            esac ;;
		--mt_reference)
            case "$2" in
                *) MT_REFERENCE=$2 ; supplied_args+=( "supplied --mt_reference: ${MT_REFERENCE}\n" ) ; shift 2 ;;
            esac ;;
		--outprefix)
            case "$2" in
                *) OUTPREFIX="$2"_ ; supplied_args+=( "supplied --outprefix: ${OUTPREFIX}\n" ) ; shift 2 ;;
            esac ;;
		--outdir)
            case "$2" in
                *) OUTDIR="$2" ; supplied_args+=( "supplied --outdir: ${OUTDIR}\n" ); shift 2 ;;
            esac ;;
		--sample_id_list)
            case "$2" in
                *) SAMPLE_IDS="$2" ; OVERRIDE_RAW="yes" ; supplied_args+=( "supplied --sample_id_list: ${SAMPLE_IDS}\n" ) ; shift 2 ;;
            esac ;;
		--stop_mira_nfs) NFS_STOP="no" ; supplied_args+=( "supplied --stop_mira_nfs: ${NFS_STOP}\n" ); shift ;;
                --resume_mitobim) RESUME_MITOBIM="yes" ; shift ;;
		--mitobim_opts)
            case "$2" in
			"all") mitobim_options=( quick 2step ) ; supplied_args+=( "supplied --mitobim_opts: ${mitobim_options[@]}\n" ); shift 2 ;;
			*) mitobim_options=("$2") ; supplied_args+=( "supplied --mitobim_opts: ${mitobim_options[@]}\n" ); shift 2 ;;
            esac ;;
		--novoplasty_opts_kmer)
            case "$2" in
                *) kmer_values=($2) ; CHECK_UNIQUE="TRUE" ; supplied_args+=( "supplied --novoplasty_opts_kmer: ${kmer_values[@]}\n" ); shift 2 ;;
            esac ;;
		--novoplasty_opts_concat)
            case "$2" in
			"all") concat=( no yes ) ; supplied_args+=( "supplied --novoplasty_opts_concat: ${concat[@]}\n" ); shift 2 ;;
                *) concat=("$2") ; supplied_args+=( "supplied --novoplasty_opts_concat: ${concat[@]}\n" ); shift 2 ;;
            esac ;;
		--novoplasty_opts_guide)
            case "$2" in
			"all") ref_guide=( no yes ) ; supplied_args+=( "supplied --novoplasty_opts_guide: ${ref_guide}\n" ); shift 2 ;;
                *) ref_guide=("$2") ; supplied_args+=( "supplied --novoplasty_opts_guide: ${ref_guide}\n" ); shift 2 ;;
            esac ;;
		--udocker) UDOCKER="yes" ; supplied_args+=( "supplied --udocker: ${UDOCKER}\n" ) ; shift ;;
		--udocker_dir)
            case "$2" in
                *) UDOCKER_DIRECTORY=$2 ; supplied_args+=( "supplied --udocker_dir: ${UDOCKER_DIRECTORY}\n" ) ; shift 2 ;;
            esac ;;
		--assembly_stats) ASSEMBLY_STATS="yes" ; supplied_args+=( "supplied --assembly_stats: ${ASSEMBLY_STATS}\n" ) ; shift ;;
		--log_pattern)
            case "$2" in
                *) LOG_PATTERN="$2" ; if [ "${ASSEMBLY_STATS}" == "yes" ] ; then supplied_args+=( "supplied --log_pattern: ${LOG_PATTERN}\n" ); fi ; shift 2 ;;
            esac ;;
		--log_pattern_exclude)
            case "$2" in
                *) LOG_PATTERN_EXCLUDE="$2" ; if [ "${ASSEMBLY_STATS}" == "yes" ] ; then supplied_args+=( "supplied --log_pattern_exclude: ${LOG_PATTERN_EXCLUDE}\n" ); fi ; shift 2 ;;
            esac ;;
		--exclude_genes_read_depth)
            case "$2" in
                *) EXCLUDE_GENES_READ_DEPTH="$2" ; if [ "${ASSEMBLY_STATS}" == "yes" ] ; then supplied_args+=( "supplied --exclude_genes_read_depth: ${EXCLUDE_GENES_READ_DEPTH}\n" ); fi ; shift 2 ;;
            esac ;;
		--skip_file_compilation) SKIP_FILE_COMPILATION="yes" ; if [ "${ASSEMBLY_STATS}" == "yes" ] ; then supplied_args+=( "supplied --skip_file_compilation: ${SKIP_FILE_COMPILATION}\n" ); fi ; shift ;;
		--stats_methods_list)
            case "$2" in
                 *) METHODS_LIST="$2" ; if [ "${ASSEMBLY_STATS}" == "yes" ] ; then supplied_args+=( "supplied --stats_methods_list: ${METHODS_LIST}\n" ); fi ; shift 2 ;;
            esac ;;
	    	--expected_genes)
		case "$2" in
			*) EXPECTED_GENES="$2" ; if [ "${ASSEMBLY_STATS}" == "yes" ] ; then supplied_args+=( "supplied --expected_genes: ${EXPECTED_GENES}\n" ); fi ; shift 2 ;;
		esac ;;
		-h|--help) HELP="yes" ; shift ;;
		--) shift ; break ;;
		*) echo "Internal error!" ; exit 1 ;;
    esac
done


# If help is requested print help and exit

if [ "${HELP}" == "yes" ]; then

echo -e "This script assembles mitochondrial genomes from NGS bycatch using MitoZ (Megahit), MITObim, and NOVOPlasty and compiles various statistics.

	-m|--method
		Determines assembly method. Can be set to: all (default), mitoz, mitobim, novoplasty, not_mitoz, not_mitobim, not_novoplasty, prep_reads, or none
        	    prep_reads: no assembly will occur but raw reads will be prepared for all methods
                    none: no assembly will occur, use this when compiling statistics on previous assemblies.
        -r|--raw_reads
		Path to folder containing raw reads. Defaults to ${PWD}/0_raw.
        -p|--prep_reads
		Boolean flag. Will prepare reads as required by selected assembly method/s.
        -d|--dedup_reads
		Boolean flag. Will TURN OFF deduplication of reads during read preparation. Use when retrimming is required but not rededuplication.
        -t|--trim_reads
		Boolean flag. Will TURN OFF trimming of reads during read preparation. Use when only deduplicating reads at this time OR when concatenation of previously trimmed reads is required for MITObim.
	-a|--trim_adaptors
		Boolean flag. Will TURN OFF trimming of only adaptors from reads during read preparation. Use when not planning NOVOPlasty assembly.
        --remove_dedup
		Boolean flag. Will REMOVE deduplicated read files in the dedup/ directory following trimming.
        --max_memory
		Set the maximum memory in GB. Supply only the number. Defaults to 15.
        --mt_annotated_ref
		Path to the annotated mitochondrial reference file. Required argument.
        --mt_reference
		Path to the unannotated mitochondrial reference file. Required argument.
        --outprefix
		Adds prefix to log files ONLY. Underscore is added to supplied text. Use to separate logs of reruns after step failures.
        --outdir
		Adds prefix and/or path to sequence output (*_bm) directories ONLY.
        --sample_id_list
		Optionally supply path to a list of sample IDs. Default will extract sample IDs from raw file names.
        --stop_mira_nfs
		Boolean flag. Will TURN ON Mira option to terminate if writing files over NFS mount.
        --resume_mitobim
		Boolean flag. Will resume previous MITObim assemblies from the latest successful iteration.
	--mitobim_opts
		Set the MITObim assembly options. Can be set to all (default), quick, or 2step.
	--novoplasty_opts_*
		Set the NOVOPlasty assembly options.	
		--novoplasty_opts_kmer
			Sets the kmer value. Can be set to any integer or list of intergers (e.g. "20 30 55"), defaults to 33 and 21.
		--novoplasty_opts_concat
			Determines inclusion of unpaired reads. Can be set to all (default), yes, or no.
		--novoplasty_opts_guide
			Determines the use of a reference guide. Can be set to all (default), yes, or no.
	--udocker
		Boolean flag. Will run MitoZ through udocker over docker.
        --udocker_dir
		Optionally specify an alternative directory to home for udocker files (e.g. a tmp directory).
        --assembly_stats
		Boolean flag. Will compile assembly statistics if supplied.
        --log_pattern
		For use with assembly stats. Will only count time elapsed log files matching supplied pattern.
        --log_pattern_exclude
		For use with assembly stats. Will NOT count time elapsed log files matching supplied pattern.
        --exclude_genes_read_depth
		For use with assembly stats. Will additionally calculate average read depth excluding suplied genes (e.g. "COX1 CYTB").
        --skip_file_compilation
		Boolean flag. For use with assembly stats. Assumes files are previously compiled and skips file compilation. Speeds up assembly stats reruns when underlying files are unchanged.
        --stats_methods_list
		For use with assembly stats. Optionally supply path to a list of folders from which to calculate statistics. Defaults to all folders with _bm suffix.
	--expected_genes
		For use with assembly stats. Optionally supply path to a file containing expected genes and their expected copy numbers. Use when expectation differs from the 37 conserved loci.
        -h|--help
		Displays this helpful text.
        " | while IFS= read -r line; do
    			if [ -z "$line" ]; then
        			echo "" # Preserves blank lines safely
    			else
				# Extract leading whitespace
        			indent="${line%%[![:space:]]*}"
        			echo "$line" | fold -s -w 120 | sed "2,\$s/^/$indent/"
    			fi
		done

	exit 0

fi

# print supplied arguments
echo -e " ${supplied_args[@]}"

# set timeformat for logging assembly times
TIMEFORMAT=`echo -e "%R\n%U\n%S"` #writes real, user, and system times on separate lines

## check if arguments without defaults are given when required

if [ ! -n "${MT_REFERENCE_ANNOTATED}" ] && [ "${ASSEMBLY_STATS}" = "yes" ] && [ ! "${METHOD}" = "prep_reads" ] ; then

	echo "ERROR: path to annotated mitochondrial reference is required"
	exit 1

fi

if [ ! -n "${MT_REFERENCE}" ] && [ ! "${METHOD}" = "prep_reads" ]; then

	echo "ERROR: path to mitochondrial reference is required"
	exit 1

fi

## check if supplied arguments are sensible

# check paths exist
if [ ! -d "${RAW_READS}" ] && [ ! "${METHOD}" = "none" ] ; then

	echo "ERROR: path to raw reads does not exist"
	exit 1

fi
if [ -n "${MT_REFERENCE_ANNOTATED}" ] && [ ! -e "${MT_REFERENCE_ANNOTATED}" ]; then

	echo "ERROR: path to annotated reference does not exist"
	exit 1

fi
if [ -n "${MT_REFERENCE}" ] && [ ! -e "${MT_REFERENCE}" ]; then

	echo "ERROR: path to reference does not exist"
	exit 1

fi
if [ -n "${SAMPLE_IDS}" ] && [ ! -f "${SAMPLE_IDS}" ]; then

	echo "ERROR: path to sample id list text file does not exist or is not a regular file"
	exit 1

fi

# check max memory is numeric
if ! [[ "${MAX_MEMORY}" =~ ^[0-9]+$ ]]; then

	echo "ERROR: supplied max memory is not an integer"
	exit 1

fi

# check mira threads is numeric
if ! [[ "${MIRA_THREADS}" =~ ^[0-9]+$ ]]; then

        echo "ERROR: supplied mira threads is not an integer"
        exit 1

fi

# translate supplied method and sanity check
case "$METHOD" in
	"all")
		NOVOPLASTY="TRUE" ; MITOZ="TRUE" ; MITOBIM="TRUE" ;;
	"not_mitobim")
		NOVOPLASTY="TRUE" ; MITOZ="TRUE" ; MITOBIM="FALSE" ;;
	"not_mitoz")
		NOVOPLASTY="TRUE" ; MITOZ="FALSE" ; MITOBIM="FALSE" ;;
	"not_novoplasty")
                NOVOPLASTY="FALSE" ; MITOZ="TRUE" ; MITOBIM="TRUE" ;;
	"novoplasty")
		NOVOPLASTY="TRUE" ; MITOZ="FALSE" ; MITOBIM="FALSE" ;;
	"mitoz")
		NOVOPLASTY="FALSE" ; MITOZ="TRUE" ; MITOBIM="FALSE" ;;
	"mitobim")
		NOVOPLASTY="FALSE" ; MITOZ="FALSE" ; MITOBIM="TRUE" ;;
	"prep_reads"|"none")
		NOVOPLASTY="FALSE" ; MITOZ="FALSE" ; MITOBIM="FALSE" ; UDOCKER="no" ;;
	*)
		echo "ERROR: inappropriate method selected" ; exit 1 ;;
esac

# check mitobim options are sensible
if [ "$MITOBIM" = "TRUE" ]; then

	# check options validity
	# store array as string
	mitobim_opts=`echo "${mitobim_options[@]}"`
	# check valid
	if [ ! "${mitobim_opts}" = "quick 2step" ] && [ ! "${mitobim_opts}" = "2step quick" ] && [ ! "${mitobim_opts}" = "2step" ] && [ ! "${mitobim_opts}" = "quick" ]; then
		echo -e "\nERROR: inappropriate mitobim options selected"
		echo -e "\n please select valid mitobim options: quick, 2step, or 'all'"
		exit 1
	fi

	# confirm array
	mitobim_options=( ${mitobim_options[@]} )

fi

# check novoplasty options are sensible
if [ "$NOVOPLASTY" = "TRUE" ]; then
	# check options validity
	# store array as string
	guide_opts=`echo "${ref_guide[@]}"`
	# check valid
	if [ ! "${guide_opts}" = "no yes" ] && [ ! "${guide_opts}" = "yes no" ] && [ ! "${guide_opts}" = "no" ] && [ ! "${guide_opts}" = "yes" ]; then
		echo -e "\nERROR: inappropriate novoplasty reference guide options selected"
		echo -e "\n please select valid novoplasty reference guide options: yes, no, or both"
		exit 1
	fi

	# check options validity
	# store array as string
	concat_opts=`echo "${concat[@]}"`
	# check valid
	if [ ! "${concat_opts}" = "no yes" ] && [ ! "${concat_opts}" = "yes no" ] && [ ! "${concat_opts}" = "no" ] && [ ! "${concat_opts}" = "yes" ]; then
		echo -e "\nERROR: inappropriate novoplasty concatenated reads options selected"
		echo -e "\n please select valid novoplasty concatenated reads options: yes, no, or both"
		exit 1
	fi

	# check that kmer_values only contain numbers
	for k in "${kmer_values[@]}"; do
    		if ! [[ $k =~ ^[0-9]+$ ]]; then
        		echo -e "\nERROR: inappropriate novoplasty kmer options selected"
			echo -e "\nplease only provide integers"
       			exit 1
    		fi
	done

	if [ "$CHECK_UNIQUE" = "TRUE" ]; then

		# confirm kmer_values are unique
		# send kmer values as string to file
		kmer_tmp=$(mktemp)
		echo "${kmer_values[@]}" > $kmer_tmp
		# separate values into lines
		sed -i 's/ /\n/g' $kmer_tmp
		kmer_tmp2=$(mktemp)
		# get unique lines
		uniq $kmer_tmp > $kmer_tmp2
		# read lines back in as array
		readarray -t kmer_values < $kmer_tmp2
		# remove tmp files
		rm ${kmer_tmp} ${kmer_tmp2}
	fi

	# confirm arrays
	kmer_values=( ${kmer_values[@]} )
	concat=( ${concat[@]} )
	ref_guide=( ${ref_guide[@]} )

# close if: NOVOPLASTY is TRUE
fi

## Update arguments as needed based on other arguments

# if method equals prep reads, set prep reads to yes

if [ "${METHOD}" = "prep_reads" ]; then

	PREP_READS="yes"

fi

# if udocker then add udocker to path and load image
if [ "${UDOCKER}" = "yes" ]; then
	export PATH="${PWD}/external_software/udocker-1.3.16/udocker:$PATH"
	
	# if udocker directory is given, write udocker files to other directory
	if [ ! "${UDOCKER_DIRECTORY}" = "" ]; then export UDOCKER_DIR="${UDOCKER_DIRECTORY}"; fi
	
	# load image from file if not existing
	if ! udocker images | grep -i -q mitoz_v3.6u ; then # check udocker image list for mitoz image
		echo -e "\nloading udocker image from ${PWD}/external_software/mitoz_v3.6u.tar\n"
		
		#load image
		udocker load -i ${PWD}/external_software/mitoz_v3.6u.tar mitoz_v3.6u
		
		# confirm image load success
		if udocker images | grep -i -q mitoz_v3.6u ; then
			echo -e "\nimage loaded successfully\n"
		else
			echo -e "\nERROR: image failed to load\n"
			exit 1
		fi
	else
		echo -e "\n udocker image exists\n"
	fi
fi

if [ ! "$METHOD" = "prep_reads" ] && [ ! "$METHOD" = "none" ] && [ "$UDOCKER" = "no" ]; then
	echo -e "\ndocker not currently supported, please enable udocker\n"
	exit 1
fi

if [ ! "$OUTDIR" = $PWD ]; then
	echo -e "\noutdirectory is : ${OUTDIR}, will create *_bm folders in working directory and move to outdirectory upon completion\n"
fi

# make subdirectory for log files
mkdir -p $PWD/logs

## create csv for logging file times ##

if ! ls $PWD/logs | grep -q ${OUTPREFIX}assembly_times_log.csv && [ ! $METHOD == "prep_reads" ]; then

	printf "sample_ID,assembly_method,start_date,start_time,end_date,end_time\n" > logs/${OUTPREFIX}assembly_times_log.csv

fi

## create csv for logging time elapsed ##

if ! ls $PWD/logs | grep -q ${OUTPREFIX}assembly_time_elapsed_log.csv && [ ! $METHOD == "prep_reads" ]; then

	printf "sample_ID,assembly_method,real,user,system\n" > logs/${OUTPREFIX}assembly_time_elapsed_log.csv

fi

## create csv for logging failed steps ##

if ! ls $PWD/logs | grep -q ${OUTPREFIX}failed_samples_log.csv ; then

	printf "sample_ID,failed_step\n" > logs/${OUTPREFIX}failed_samples_log.csv

fi

#store mira path
mirapath="$PWD/external_software/MITObim/docker/external_software/mira_4.0.2"

# add mira to path
export PATH="${mirapath}:$PATH"

# if method is not none, get sample ids
if [ ! "$METHOD" = "none" ]; then

## get input read file names ##

	if [ $OVERRIDE_RAW = "no" ]; then

		echo -e "\nretrieving sample_IDs from raw file names in ${RAW_READS}\n"

		# store raw_file names
		file_names=( $(ls "${RAW_READS}") )

		# get R1 file names and store as string; add .[IDENTIFIER] to _R1 to access certain type of reads e.g. trimmed vs untrimmed; paired vs unpaired
		R1_files=$(printf -- '%s\n' "${file_names[@]}" | grep _R1 | tr '\n' ' ')

		# convert string to array
		sample_array=($R1_files)

	else # if override_raw is yes

		echo -e "\nretrieving sample_IDs from supplied text file\n"

		# read sample_IDs from supplied text file
		readarray -t sample_array < "$SAMPLE_IDS"

	fi

	# loop over sample array to perform read prep and assembly as directed
	for i in "${sample_array[@]}";
	do
		if [ $OVERRIDE_RAW = "no" ]; then
			sample_ID="${i%_R1.*}" # extract sample id from stored file names
		else
			sample_ID="$i" # sample id supplied directly via text file
		fi

#############################################################################
#
# PREPARATION OF RAW READS
#
# deduplication and trimming as appropriate for mitoz, mitobim, and/or novoplasty
#
#############################################################################

		if [ "$PREP_READS" = "yes" ]; then

			echo -e "\npreparing reads for assembly and annotation as required by selected methods\n"

			# create csv for logging read prep times #

			if ! ls $PWD/logs | grep -q ${OUTPREFIX}read_prep_times_log.csv; then

				printf "sample_ID,preparation_step,start_date,start_time,end_date,end_time\n" > logs/${OUTPREFIX}read_prep_times_log.csv

			fi

			## create csv for logging time elapsed ##

			if ! ls $PWD/logs | grep -q ${OUTPREFIX}read_prep_time_elapsed_log.csv ; then

				printf "sample_ID,preparation_step,real,user,system\n" > logs/${OUTPREFIX}read_prep_time_elapsed_log.csv

			fi

###########     Deduplication of raw reads     ####################################################	

			if [ "$DEDUP_READS" = "yes" ]; then

				preparation_step="dedup"
				echo -e "\nread preparation step: $preparation_step\n"

				# make directory #
				mkdir -p $PWD/dedup

				# store start date and time #
				start_date=$(date +"%Y-%m-%d")
				start_time=$(date +"%H:%M:%S")

				# begin tracking of elapsed time #
				time (
				{

				# get appropriate file path #
				file_R1=$(find ${RAW_READS} | grep "$sample_ID"_R1);
				file_R2=$(find ${RAW_READS} | grep "$sample_ID"_R2);

				echo -e "\ndeduplicate $sample_ID reads"
				echo "output stored in ${PWD}/dedup"
				echo "R1 file: $file_R1"
				echo -e "R2 file: $file_R2\n"
				
				# dedup reads
				external_software/bbmap/dedupe.sh -Xmx"$MAX_MEMORY"000m in1=$file_R1 in2=$file_R2 out=$PWD/dedup/"$sample_ID".deduped.fastq.gz pigz=t unpigz=t overwrite=true
				#bash dedupe.sh -Xmx"$MAX_MEMORY"000m in1=$file_R1 in2=$file_R2 out=$PWD/dedup/"$sample_ID".deduped.fastq.gz pigz=t unpigz=t overwrite=true

				# check for success #
			
				if ! ls $PWD/dedup | grep -q "$sample_ID".deduped.fastq.gz ; then
				
					echo -e "\nWARNING: deduplication of $sample_ID failed, trying 'repair.sh' with 'nullifybrokenquality'...\n"
					printf "$sample_ID,WARNING: deduplication of $sample_ID failed trying 'repair.sh' with 'nullifybrokenquality'...\n" >> logs/${OUTPREFIX}check_warnings.csv

					external_software/bbmap/repair.sh -Xmx"$MAX_MEMORY"000m in1=$file_R1 in2=$file_R2 out=$PWD/dedup/"$sample_ID".repaired.fastq.gz nullifybrokenquality=true pigz=t unpigz=t overwrite=true
					#bash repair.sh -Xmx"$MAX_MEMORY"000m in1=$file_R1 in2=$file_R2 out=$PWD/dedup/"$sample_ID".repaired.fastq.gz nullifybrokenquality=true pigz=t unpigz=t overwrite=true

					external_software/bbmap/dedupe.sh -Xmx"$MAX_MEMORY"000m in=$PWD/dedup/"$sample_ID".repaired.fastq.gz out=$PWD/dedup/"$sample_ID".deduped.fastq.gz pigz=t unpigz=t overwrite=true
					#bash dedupe.sh -Xmx"$MAX_MEMORY"000m in=$PWD/dedup/"$sample_ID".repaired.fastq.gz out=$PWD/dedup/"$sample_ID".deduped.fastq.gz pigz=t unpigz=t overwrite=true

					rm  $PWD/dedup/"$sample_ID".repaired.fastq.gz

					if ! ls $PWD/dedup | grep -q "$sample_ID".deduped.fastq.gz ; then

						printf "$sample_ID,$preparation_step\n" >> logs/${OUTPREFIX}failed_samples_log.csv
						printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
						echo -e "\nERROR: deduplication of $sample_ID still failed, starting next sample\n"
						continue
					fi

				fi

				if ls $PWD/dedup | grep -q "$sample_ID".deduped.fastq.gz ; then

					# reformat deduped reads as separate R1 and R2 files
					external_software/bbmap/reformat.sh in=$PWD/dedup/"$sample_ID".deduped.fastq.gz out1=$PWD/dedup/"$sample_ID"_R1.DR.fastq out2=$PWD/dedup/"$sample_ID"_R2.DR.fastq pigz=t unpigz=t overwrite=true
					#bash reformat.sh in=$PWD/dedup/"$sample_ID".deduped.fastq.gz out1=$PWD/dedup/"$sample_ID"_R1.DR.fastq out2=$PWD/dedup/"$sample_ID"_R2.DR.fastq pigz=t unpigz=t overwrite=true

					# remove combined deduped reads file
					rm  $PWD/dedup/"$sample_ID".deduped.fastq.gz

					echo -e "\ndeduplication concluded\n"

					# store end date and time #

					end_date=$(date +"%Y-%m-%d")
					end_time=$(date +"%H:%M:%S")

					# log start and end times #

					printf "$sample_ID,$preparation_step,$start_date,$start_time,$end_date,$end_time\n" >> logs/${OUTPREFIX}read_prep_times_log.csv
				fi

				# end tracking of elapsed time #
				} 2>&1
				) 2>${OUTPREFIX}time_elapsed.log

				# store elapsed time
				readarray -t time_elapsed < ${OUTPREFIX}time_elapsed.log
				printf "$sample_ID,$preparation_step,${time_elapsed[0]},${time_elapsed[1]},${time_elapsed[2]}\n" >> logs/${OUTPREFIX}read_prep_time_elapsed_log.csv
				rm ${OUTPREFIX}time_elapsed.log
			
			# close if: dedup reads is yes
			fi

			if [ "$DEDUP_READS" = "no" ]; then

				echo -e "\nskipping deduplication of reads and assuming previously deduped reads are in /dedup\n"

			fi

###########     Trimming of deduplicated reads     ###########################################

			if [ "$TRIM_READS" = "yes" ]; then

				preparation_step="trim"
				echo -e "\nread preparation step: $preparation_step\n"

				# make directory #
				mkdir -p $PWD/deduptrim
				mkdir -p $PWD/stats

				# store start date and time #
				start_date=$(date +"%Y-%m-%d")
				start_time=$(date +"%H:%M:%S")

				# begin tracking of elapsed time #
				time (
				{

				# get appropriate file path #
				file_R1=$(find $PWD/dedup/ | grep "$sample_ID"_R1);
				file_R2=$(find $PWD/dedup/ | grep "$sample_ID"_R2);

				echo -e "\ntrim $sample_ID reads"
				echo "output stored in $PWD/deduptrim"
				echo "R1 file: $file_R1"
				echo -e "R2 file: $file_R2\n"

				# trim reads
				java -Xmx"$MAX_MEMORY"000m -jar $PWD/external_software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 "$file_R1" "$file_R2" $PWD/deduptrim/"$sample_ID"_R1.DRTRIMP.fastq.gz $PWD/deduptrim/"$sample_ID"_R1.DRTRIMU.fastq.gz $PWD/deduptrim/"$sample_ID"_R2.DRTRIMP.fastq.gz $PWD/deduptrim/"$sample_ID"_R2.DRTRIMU.fastq.gz ILLUMINACLIP:$PWD/external_software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:1:35:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:8:20 MINLEN:40 2>&1 >/dev/null | awk '/Input Read Pairs:/ {print $4,$7,$8,$12,$13,$17,$18,$20,$21}' OFS="\t" > $PWD/stats/"$sample_ID".subtrimsumPE.log

				# adapter trimming only
				#java -Xmx"$MAX_MEMORY"000m -jar $PWD/external_software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 "$file_R1" "$file_R2" $PWD/deduptrim/"$sample_ID"_R1.DRTRIMP.fastq.gz $PWD/deduptrim/"$sample_ID"_R1.DRTRIMU.fastq.gz $PWD/deduptrim/"$sample_ID"_R2.DRTRIMP.fastq.gz $PWD/deduptrim/"$sample_ID"_R2.DRTRIMU.fastq.gz ILLUMINACLIP:$PWD/external_software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:1:35:15 2>&1 >/dev/null | awk '/Input Read Pairs:/ {print $4,$7,$8,$12,$13,$17,$18,$20,$21}' OFS="\t" > $PWD/stats/"$sample_ID".subtrimsumPE.log
				
				# check for success				

				if ! ls $PWD/deduptrim | grep -q "$sample_ID"_R1.DRTRIMP.fastq.gz ; then

					printf "$sample_ID,$preparation_step\n" >> logs/${OUTPREFIX}failed_samples_log.csv
					printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
					deduptrim_failed="true"
					echo -e "\nERROR: trimming of $sample_ID failed, skipping to next independent step\n"

				fi

				# end tracking of elapsed time #
				} 2>&1
				) 2>${OUTPREFIX}time_elapsed.log

				# store elapsed time
				readarray -t time_elapsed < ${OUTPREFIX}time_elapsed.log
				printf "$sample_ID,$preparation_step,${time_elapsed[0]},${time_elapsed[1]},${time_elapsed[2]}\n" >> logs/${OUTPREFIX}read_prep_time_elapsed_log.csv
				rm ${OUTPREFIX}time_elapsed.log

				if ls $PWD/deduptrim | grep -q "$sample_ID"_R1.DRTRIMP.fastq.gz ; then
				
					echo "trimming concluded"
	
					# store end date and time #
	
					end_date=$(date +"%Y-%m-%d")
					end_time=$(date +"%H:%M:%S")

					# log start and end times #

					printf "$sample_ID,$preparation_step,$start_date,$start_time,$end_date,$end_time\n" >> logs/${OUTPREFIX}read_prep_times_log.csv
						
				fi
				
				# if directed to remove dedup files, remove here
                if [ "${REMOVE_DEDUP}" = "yes" ]; then
					rm ${PWD}/dedup/${sample_ID}_R1.DR.fastq
					rm ${PWD}/dedup/${sample_ID}_R2.DR.fastq
                fi
			
			#close if: trim reads is yes
			fi
			
			if [ "$TRIM_READS" = "no" ]; then

				echo -e "\nskipping trimming of reads and assuming previously trimmed reads are in /deduptrim. Reads will still be concatenated into /deduptrimCat for MITObim\n"

			fi
			
###################################     Concatenation of trimmed reads     ########################

			# concatenation will still occur if prep reads is true but dedup reads and trim reads are false

			# if deduptrim files exist
			if ls $PWD/deduptrim | grep -q "$sample_ID"_R1.DRTRIMP.fastq.gz ; then
			
				# and if mitobim is true or method is prep reads
				if [ "$MITOBIM" = "TRUE" ] || [ "$METHOD" = "prep_reads" ]; then

					preparation_step="cat_trim"
					echo -e "\nread preparation step: $preparation_step\n"

					# make directory #
					mkdir -p $PWD/deduptrimCat

					# store start date and time #
					start_date=$(date +"%Y-%m-%d")
					start_time=$(date +"%H:%M:%S")

					# begin tracking of elapsed time #
					time (
					{
				
					# get appropriate file path #
					file_R1_P=$(find $PWD/deduptrim/ | grep "$sample_ID"_R1.DRTRIMP.fastq);
					file_R2_P=$(find $PWD/deduptrim/ | grep "$sample_ID"_R2.DRTRIMP.fastq);
					file_R1_U=$(find $PWD/deduptrim/ | grep "$sample_ID"_R1.DRTRIMU.fastq);
					file_R2_U=$(find $PWD/deduptrim/ | grep "$sample_ID"_R2.DRTRIMU.fastq);


					echo -e "\nconcatenate deduptrim $sample_ID reads"
					echo -e "output stored in $PWD/deduptrimCat\n"
	
					echo -e "\nconcatenating: \n$file_R1_P\n$file_R1_U\n$file_R2_P\n$file_R2_U\n"
	
					cat "$file_R1_P" "$file_R2_P" "$file_R1_U" "$file_R2_U" > $PWD/deduptrimCat/"$sample_ID".DRTRIMC.fastq.gz
	
					echo -e "\nconcatenation concluded\n"

					# end tracking of elapsed time #
					} 2>&1
					) 2>${OUTPREFIX}time_elapsed.log

					# store elapsed time
					readarray -t time_elapsed < ${OUTPREFIX}time_elapsed.log
					printf "$sample_ID,$preparation_step,${time_elapsed[0]},${time_elapsed[1]},${time_elapsed[2]}\n" >> logs/${OUTPREFIX}read_prep_time_elapsed_log.csv
					rm ${OUTPREFIX}time_elapsed.log				

					# store end date and time #
	
					end_date=$(date +"%Y-%m-%d")
					end_time=$(date +"%H:%M:%S")
	
					# log start and end times #
	
					printf "$sample_ID,$preparation_step,$start_date,$start_time,$end_date,$end_time\n" >> logs/${OUTPREFIX}read_prep_times_log.csv
				# close if: method is mitobim or prep_reads
				fi
				
			# close if: deduptrim files exists
			fi
			
			if [ "$TRIM_ADAPTORS" = "yes" ]; then

					preparation_step="trim_adaptors"
					echo -e "\nread preparation step: $preparation_step\n"
	
					# make directory #
					mkdir -p $PWD/deduptrimA
					mkdir -p $PWD/stats

					# store start date and time #
					start_date=$(date +"%Y-%m-%d")
					start_time=$(date +"%H:%M:%S")

					# begin tracking of elapsed time #
					time (
					{
	
					# get appropriate file path #
					file_R1=$(find $PWD/dedup/ | grep "$sample_ID"_R1);
					file_R2=$(find $PWD/dedup/ | grep "$sample_ID"_R2);
	
					echo -e "\ntrim $sample_ID reads"
					echo "output stored in $PWD/deduptrimA"
					echo "R1 file: $file_R1"
					echo -e "R2 file: $file_R2\n"

					java -Xmx"$MAX_MEMORY"000m -jar $PWD/external_software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 "$file_R1" "$file_R2" $PWD/deduptrimA/"$sample_ID"_R1.DRTRIMAP.fastq.gz $PWD/deduptrimA/"$sample_ID"_R1.DRTRIMAU.fastq.gz $PWD/deduptrimA/"$sample_ID"_R2.DRTRIMAP.fastq.gz $PWD/deduptrimA/"$sample_ID"_R2.DRTRIMAU.fastq.gz ILLUMINACLIP:$PWD/external_software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:1:35:15 2>&1 >/dev/null | awk '/Input Read Pairs:/ {print $4,$7,$8,$12,$13,$17,$18,$20,$21}' OFS="\t" > $PWD/stats/"$sample_ID".subtrimAsumPE.log

					# check for success				

					if ! ls $PWD/deduptrimA | grep -q "$sample_ID"_R1.DRTRIMAP.fastq.gz ; then

						printf "$sample_ID,$preparation_step\n" >> logs/${OUTPREFIX}failed_samples_log.csv
						printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
						deduptrimA_failed="true"
						echo -e "\nERROR: adaptor only trimming of $sample_ID failed, skipping to next independent step\n"

					fi

					# end tracking of elapsed time #
					} 2>&1
					) 2>${OUTPREFIX}time_elapsed.log

					# store elapsed time
					readarray -t time_elapsed < ${OUTPREFIX}time_elapsed.log
					printf "$sample_ID,$preparation_step,${time_elapsed[0]},${time_elapsed[1]},${time_elapsed[2]}\n" >> logs/${OUTPREFIX}read_prep_time_elapsed_log.csv
					rm ${OUTPREFIX}time_elapsed.log

					if ls $PWD/deduptrimA | grep -q "$sample_ID"_R1.DRTRIMAP.fastq.gz ; then

						echo -e "\nadaptor trimming concluded\n"
	
						# store end date and time #
						end_date=$(date +"%Y-%m-%d")
						end_time=$(date +"%H:%M:%S")

						# log start and end times #
						printf "$sample_ID,$preparation_step,$start_date,$start_time,$end_date,$end_time\n" >> logs/${OUTPREFIX}read_prep_times_log.csv

###################################     Concatenation of adaptor trimmed reads     ################

						preparation_step="cat_trimA"
						echo -e "\nread preparation step: $preparation_step\n"

						# make directory #
						mkdir -p $PWD/deduptrimACat

						# store start date and time #
						start_date=$(date +"%Y-%m-%d")
						start_time=$(date +"%H:%M:%S")

						# begin tracking of elapsed time #
						time (
						{
				
						# get appropriate file path #
						file_R1_P=$(find $PWD/deduptrimA/ | grep "$sample_ID"_R1.DRTRIMAP.fastq);
						file_R2_P=$(find $PWD/deduptrimA/ | grep "$sample_ID"_R2.DRTRIMAP.fastq);
						file_R1_U=$(find $PWD/deduptrimA/ | grep "$sample_ID"_R1.DRTRIMAU.fastq);
						file_R2_U=$(find $PWD/deduptrimA/ | grep "$sample_ID"_R2.DRTRIMAU.fastq);

						echo -e "\nconcatenate paired and unpaired deduptrimA $sample_ID reads"
						echo -e "output stored in $PWD/deduptrimACat\n"
	
						echo -e "\nconcatenating: \n$file_R1_P\n$file_R1_U\n"

						cat "$file_R1_P" "$file_R1_U" > $PWD/deduptrimACat/"$sample_ID"_R1.DRTRIMAC.fastq.gz

						echo -e "\nconcatenating: \n$file_R2_P\n$file_R2_U\n"

						cat "$file_R2_P" "$file_R2_U" > $PWD/deduptrimACat/"$sample_ID"_R2.DRTRIMAC.fastq.gz
	
						echo -e "\nconcatenation concluded\n"

						# end tracking of elapsed time #
						} 2>&1
						) 2>${OUTPREFIX}time_elapsed.log

						# store elapsed time
						readarray -t time_elapsed < ${OUTPREFIX}time_elapsed.log
						printf "$sample_ID,$preparation_step,${time_elapsed[0]},${time_elapsed[1]},${time_elapsed[2]}\n" >> logs/${OUTPREFIX}read_prep_time_elapsed_log.csv
						rm ${OUTPREFIX}time_elapsed.log				

						# store end date and time #
	
						end_date=$(date +"%Y-%m-%d")
						end_time=$(date +"%H:%M:%S")
	
						# log start and end times #
	
						printf "$sample_ID,$preparation_step,$start_date,$start_time,$end_date,$end_time\n" >> logs/${OUTPREFIX}read_prep_times_log.csv
					# close if: deduptrimA reads exist
					fi
				# close if: trim adaptors is yes
				fi
				
				if [ "$TRIM_ADAPTORS" = "no" ]; then

					echo -e "\nskipping trimming of only adaptors from reads and assuming previously trimmed reads are in /deduptrimA - required for novoplasty\n"

				fi
			
		# close if: prep reads is yes
		fi

# RUN TEST: read preparation		
#done # close loop: sample ids
#fi # close id: method is NOT none
#exit 0
# END TEST: read preparation
		
#############################################################################
#
# MITOGENOME ASSEMBLY AND ANNOTATION
#
# assembly via mitoz, mitobim, and/or novoplasty; annotation via mitoz
#
#############################################################################

		if [ ! "${METHOD}" = "prep_reads" ] ; then
			echo -e "\nAssembling: ${sample_ID}\n"
			
###     MitoZ assembly and annotation     #########################################################

			if [ "$MITOZ" = "TRUE" ]; then
			
				# check existence of dependent files, if previous dependent step not failed log step failure
				
				if ! ls $PWD/deduptrim | grep -q "$sample_ID"_R1.DRTRIMP.f && [ ! -n "$deduptrim_failed" ]; then

					printf "$sample_ID,mitoz_assembly\n" >> logs/${OUTPREFIX}failed_samples_log.csv
					printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
					echo -e "\nERROR: mitoz assembly failed as deduplicated and trimmed $sample_ID files not found in /deduptrim, skipping to next independent step\n"

				fi
				
				if ls $PWD/deduptrim | grep -q "$sample_ID"_R1.DRTRIMP.f ; then

					echo -e "\nAssembly method: mitoz\n"

					# store start date and time #
					start_date=$(date +"%Y-%m-%d")
					start_time=$(date +"%H:%M:%S")

					# begin tracking of elapsed time #
					time (
					{

					# get appropriate file path #
					file_R1=$(find $PWD/deduptrim/ | grep "$sample_ID"_R1.DRTRIMP.f);
					file_R2=$(find $PWD/deduptrim/ | grep "$sample_ID"_R2.DRTRIMP.f);
	
					# make directory #
					mkdir -p $PWD/mitoz_bm

					# decompress .gz files #

					# check file types
					file_type_R1="${file_R1##*.}"
					file_type_R2="${file_R2##*.}"

					if [ "$file_type_R1" = "gz" ]; then

						# decompress .gz file
						gzip -dk $file_R1
						# get decompressed file name
						file_R1="${file_R1%.gz}"

					fi

					if [ "$file_type_R2" = "gz" ]; then

						# decompress .gz file
						gzip -dk $file_R2
						# get decompressed file name
						file_R2="${file_R2%.gz}"

					fi

					# make subdirectory #

					mkdir -p $PWD/mitoz_bm/"$sample_ID"

					# run mitoz assembly module

					if [ $UDOCKER = "yes" ]; then
						# run mitoz #
						udocker run --volume=$PWD:$PWD --workdir=$PWD --rm mitoz_v3.6u:latest  mitoz assemble --workdir $PWD/mitoz_bm/"$sample_ID" --fq1 "$file_R1" --fq2 "$file_R2" --clade Chordata --outprefix $sample_ID --insert_size 300 --requiring_taxa Chordata --assembler megahit --kmers_megahit 59 79 99 119 141
					fi

					# check success

					mitogenome_assembly_file_name=`ls $PWD/mitoz_bm/${sample_ID}/${sample_ID}*.result | grep "mitogenome.fa"`
					if [ ! "${mitogenome_assembly_file_name}" = "" ]; then mitogenome_assembly_path=`ls $PWD/mitoz_bm/${sample_ID}/${sample_ID}*.result/${mitogenome_assembly_file_name}` ; else mitogenome_assembly_path=""; fi

					if [ ! -s "${mitogenome_assembly_path}" ]; then

						echo -e "\nWARNING: mitoz assembly of $sample_ID failed, trying --min_abundance 0 AND --slow_search\n"
						printf "${sample_ID},WARNING: mitoz assembly of ${sample_ID} failed trying --min_abundance 0 AND --slow_search\n" >> logs/${OUTPREFIX}check_warnings.csv

						# remove previous attempt as temporary files can cause rerun issues
						if ls $PWD/mitoz_bm | grep -q "$sample_ID"; then
							rm -r $PWD/mitoz_bm/$sample_ID
							echo -e "\nremoved $PWD/mitoz_bm/$sample_ID directory in preparation for new attempt\n"
						fi

						# recreate subdirectory
						mkdir -p $PWD/mitoz_bm/"$sample_ID"

						if [ $UDOCKER = "yes" ]; then
							# run mitoz #
							udocker run --volume=$PWD:$PWD --workdir=$PWD --rm mitoz_v3.6u:latest  mitoz assemble --workdir $PWD/mitoz_bm/"$sample_ID" --fq1 "$file_R1" --fq2 "$file_R2" --clade Chordata --outprefix $sample_ID --insert_size 300 --requiring_taxa Chordata --assembler megahit --kmers_megahit 59 79 99 119 141 --min_abundance 0 --slow_search
						fi
					
						# retrieve file name and path as previous failure means variables are blank
						mitogenome_assembly_file_name=`ls $PWD/mitoz_bm/${sample_ID}/${sample_ID}*.result | grep "mitogenome.fa"`
						if [ ! "${mitogenome_assembly_file_name}" = "" ]; then mitogenome_assembly_path=`ls $PWD/mitoz_bm/${sample_ID}/${sample_ID}*.result/${mitogenome_assembly_file_name}` ; else mitogenome_assembly_path=""; fi

						# check success
						if [ ! -s "${mitogenome_assembly_path}" ]; then
							printf "$sample_ID,mitoz_assembly\n" >> logs/${OUTPREFIX}failed_samples_log.csv
							printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
							echo -e "\nERROR: mitoz assembly failed\n"
						fi

					fi
					
					if [ -s "${mitogenome_assembly_path}" ]; then

						# Check if high abundance contigs have been dropped due to low score and add to fasta file
						# get path to 'not-picked' contigs
						dropped_contigs_path=`ls $PWD/mitoz_bm/${sample_ID}/${sample_ID}*.result/*like.sorted.Not-picked.fa`
						if [ -s ${dropped_contig_path} ] ; then 
							sed 's/abun=.*/topology=linear/g' ${dropped_contig_path} >> ${mitogenome_assembly_path}
							printf "${sample_ID},WARNING: restored low scoring high abundance contigs dropped by mitoz\n" >> logs/${OUTPREFIX}check_warnings.csv
						fi


						echo -e "\nmitoz assembly concluded\n"

						# REMOVE HIGH SPACE INTERMEDIARAY FILES
						if [ -d $PWD/mitoz_bm/"$sample_ID"/metaspades ]; then
						       rm -r $PWD/mitoz_bm/"$sample_ID"/metaspades/corrected $PWD/mitoz_bm/"$sample_ID"/metaspades/K21 $PWD/mitoz_bm/"$sample_ID"/metaspades/K33 $PWD/mitoz_bm/"$sample_ID"/metaspades/K55
						fi

						# store fasta_file path
						fasta_file="${mitogenome_assembly_path}"

						# replace space before topology with _; prevents header and topology being read as separate array objects
						sed -i 's/ topology/_topology/' $fasta_file

						# store contig header names
						header_lines=$(printf -- '%s\n' | grep -i ">" $fasta_file | tr '\n' ' ')

						# convert string to array
						header_lines=($header_lines)

						# check if only one header and topology is circular
						if [ "${#header_lines[@]}" = 1 ] && echo "${header_lines[@]}" | grep -q "topology=circular" ; then

							echo -e "\nassembly is circular, checking alignment to reference and reordering nucleotides if needed\n"

							# rename contig header to remove space
							sed -i 's/ topology=circular/_topology=circular/' ${fasta_file}
															
							# run reorder nucleotides -> uses mafft to ensure nucleotides start at origin where possible
							bash auto_scripts/reorder_nucleotides.sh -f $fasta_file -r $MT_REFERENCE -o "$PWD/mitoz_bm/${sample_ID}" -c "yes"

							# replace fasta_file with reordered file if existing -> will not exist if reordering is not needed
							if [ -s "$PWD/mitoz_bm/${sample_ID}/realigned.fasta" ]; then
								# rename file
								mv $PWD/mitoz_bm/${sample_ID}/realigned.fasta $PWD/mitoz_bm/${sample_ID}/${sample_ID}_realigned.fasta
	
								# restore space in contig header in original fasta
								sed -i 's/_topology=circular/ topology=circular/' ${fasta_file}

								# reassign fasta_file
								fasta_file=$PWD/mitoz_bm/${sample_ID}/${sample_ID}_realigned.fasta

								# restore space in contig header in realigned fasta
								sed -i 's/_topology=circular/ topology=circular/' ${fasta_file}

							fi

							# if realigned.fasta does not exist, restore space in original contig header
							if [ ! -s "$PWD/mitoz_bm/${sample_ID}/realigned.fasta" ]; then
								# restore space in contig header in original fasta
								sed -i 's/_topology=circular/ topology=circular/' ${fasta_file}
							fi

							# check if max iteration reached warning file exists and write to check warnings.csv
							if [ -s "$PWD/mitoz_bm/${sample_ID}/WARNING_reorder_nucleotides.txt" ]; then
								echo -e "\nWARNING: reordering of nucleotides for ${fasta_file} reached maximum number of iterations, reordering may be incomplete\n"
								printf "$sample_ID,WARNING: reordering of nucleotides for ${fasta_file} reached maximum number of iterations reordering may be incomplete\n" >> logs/${OUTPREFIX}check_warnings.csv
							fi
									
						fi

###################### START circular but multiple headers

						# check if circular topology but multiple headers
						if [ "${#header_lines[@]}" > 1 ] && echo "${header_lines[@]}" | grep -q "topology=circular" ; then

							echo -e "\nassembly contains circular contigs, checking alignment of circular contigs to reference and reordering nucleotides if needed\n"

							# remove any existing circular contigs file
							if [ -s mitoz_bm/${sample_ID}/circular_contigs.fasta ] ; then
								rm mitoz_bm/${sample_ID}/circular_contigs.fasta
							fi

							# get line numbers of headers with circular topology
							circular_header_line_numbers=$(printf -- '%s\n' | grep -n -i "topology=circular" $fasta_file | cut -d : -f 1 | tr '\n' ' ')

							# get line numbers with linear topology
							linear_header_line_numbers=$(printf -- '%s\n' | grep -n -i "topology=linear" $fasta_file | cut -d : -f 1 | tr '\n' ' ')

							# loop through circular headers
							for c in ${circular_header_line_numbers[@]};
							do

								# get tail of file from header
								tail -n "+$c" ${fasta_file} > mitoz_bm/${sample_ID}/tmp_circ.fasta

								# find line number of next header
								scndheader=`grep -n -i -m2 ">" mitoz_bm/${sample_ID}/tmp_circ.fasta | cut -d : -f 1 | tail -n1`

								# check if line number of second header equals first line
								# if no, then take head lines of file up to before second header and append to circular contigs file
								if [ ! $scndheader = 1 ]; then
									head -n $((${scndheader}-1)) mitoz_bm/${sample_ID}/tmp_circ.fasta >> mitoz_bm/${sample_ID}/circular_contigs.fasta
								fi

								# if yes, then append tmp file to circular contigs file
								if [ $scndheader = 1 ]; then
									cat mitoz_bm/${sample_ID}/tmp_circ.fasta >> mitoz_bm/${sample_ID}/circular_contigs.fasta
								fi

								# remove tmp circ file
								rm mitoz_bm/${sample_ID}/tmp_circ.fasta

							done

							# store circular_contigs
							circular_file="mitoz_bm/${sample_ID}/circular_contigs.fasta"

							# rename contig headers in circular_contigs file to remove space
							sed -i 's/ topology=circular/_topology=circular/' ${circular_file}

							# run reorder nucleotides with circular enabled -> uses mafft to ensure nucleotides start at origin where possible
							bash auto_scripts/reorder_nucleotides.sh -f $circular_file -r $MT_REFERENCE -o "$PWD/mitoz_bm/${sample_ID}" -c "yes"

							# if realigned fasta file exists
							if [ -s "$PWD/mitoz_bm/${sample_ID}/realigned.fasta" ]; then

								# rename file
								mv $PWD/mitoz_bm/${sample_ID}/realigned.fasta $PWD/mitoz_bm/${sample_ID}/${sample_ID}_realigned.fasta

								# restore space in contig header in circular fasta
								sed -i 's/_topology=circular/ topology=circular/' ${circular_file}

								# loop through linear headers and add to realigned file
								for l in ${linear_header_line_numbers[@]};
								do

									# get tail of file from header
									tail -n "+$l" ${fasta_file} > mitoz_bm/${sample_ID}/tmp_linear.fasta

									# find line number of next header
									scndheader=`grep -n -i -m2 ">" mitoz_bm/${sample_ID}/tmp_linear.fasta | cut -d : -f 1 | tail -n1`

									# check if line number of second header equals first line
									# if no, then take head lines of file up to before second header and append to realigned_fasta file
									if [ ! $scndheader = 1 ]; then
										head -n $((${scndheader}-1)) mitoz_bm/${sample_ID}/tmp_linear.fasta >> $PWD/mitoz_bm/${sample_ID}/${sample_ID}_realigned.fasta
									fi

									# if yes, then append tmp file to realigned_fasta file
									if [ $scndheader = 1 ]; then
										cat mitoz_bm/${sample_ID}/tmp_linear.fasta >> $PWD/mitoz_bm/${sample_ID}/${sample_ID}_realigned.fasta
									fi

									# remove tmp linear file
									rm mitoz_bm/${sample_ID}/tmp_linear.fasta

								done

								# reassign fasta_file to realigned file
								fasta_file=$PWD/mitoz_bm/${sample_ID}/${sample_ID}_realigned.fasta

								# restore space in contig headers in realigned fasta
								sed -i 's/_topology=circular/ topology=circular/' ${fasta_file}

							fi

							# if realigned.fasta does not exist, remove circular fasta file
							if [ ! -s "$PWD/mitoz_bm/${sample_ID}/realigned.fasta" ]; then
								rm ${circular_file}
							fi

							# check if max iteration reached warning file exists and write to check warnings.csv
							if [ -s "$PWD/mitoz_bm/${sample_ID}/WARNING_reorder_nucleotides.txt" ]; then
								echo -e "\nWARNING: reordering of nucleotides for ${fasta_file} reached maximum number of iterations, reordering may be incomplete\n"
								printf "$sample_ID,WARNING: reordering of nucleotides for ${fasta_file} reached maximum number of iterations reordering may be incomplete\n" >> logs/${OUTPREFIX}check_warnings.csv
							fi
						fi


###################### END circular but multiple headers

						# check if contigs need renaming

						# if contigs contain _R_ prefix from alignment, remove suffix
						if grep -q ">_R_" ${fasta_file} ; then
							sed -i 's/>_R_/>/' $fasta_file;
						fi

						# if any contigs contain _FivePCGs suffix, rename contigs
						if grep -q "_FivePCGs" ${fasta_file} ; then

							# confirm space between topology and header is marked with _ -> will not be if realigned fasta
							sed -i 's/ topology=circular/_topology=circular/' ${fasta_file}

							# store header line text
							header_lines=$(printf -- '%s\n' | grep -i ">" $fasta_file | tr '\n' ' ')
							header_lines=( $header_lines )

							# store header line numbers
							header_line_numbers=$(printf -- '%s\n' | grep -n -i ">" $fasta_file | cut -d : -f 1 | tr '\n' ' ')
							header_line_numbers=( $header_line_numbers )

							# loop over header line array indices
							for l in ${!header_lines[@]};
							do
								# store topology information
								topology=`echo ${header_lines[$l]} | rev | cut -d _ -f 1 | rev`

								# if header contains _FivePCGs suffix, rename as >k_*_5PCGs
								if echo ${header_lines[$l]} | grep -q "_FivePCGs" ; then
										sed -i "${header_line_numbers[$l]} s/>.*/>k_${l}_5PCGs ${topology}/" $fasta_file;
								# if header does not contain _FivePCGs suffix, rename as >k_*
								else							
										sed -i "${header_line_numbers[$l]} s/>.*/>k_${l} ${topology}/" $fasta_file;
								fi
							done
						fi

						# restore space before topology if existing
						sed -i 's/_topology/ topology/' $fasta_file


						echo -e "\nannotating ${fasta_file}\n"

						# make annotation subdirectory
						mkdir -p "$PWD"/mitoz_bm/"${sample_ID}"/mt_annotation
					
						# run mitoz annotation module

						if [ $UDOCKER = "yes" ]; then
							# run mitoz #
							udocker run --volume=$PWD:$PWD --workdir=$PWD --rm mitoz_v3.6u:latest  mitoz annotate --workdir ./mitoz_bm/"$sample_ID"/mt_annotation --fastafiles "$fasta_file" --fq1 "$file_R1" --fq2 "$file_R2" --clade Chordata --outprefix "$sample_ID" --species_name "$sample_ID"
						fi

						# check annotation success
						if ! ls $PWD/mitoz_bm/"$sample_ID"/mt_annotation/"$sample_ID"."${fasta_file##*/}".result | grep -q summary.txt; then

							printf "$sample_ID,mitoz_annotation\n" >> logs/${OUTPREFIX}failed_samples_log.csv
							printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
							echo -e "\nERROR: mitoz annotation failed\n"

						fi
						if ls $PWD/mitoz_bm/"$sample_ID"/mt_annotation/"$sample_ID"."${fasta_file##*/}".result | grep -q summary.txt; then
							echo -e "\nmitoz annotation concluded\n"
						fi
					fi

					# remove decompressed files #

					if [ "$file_type_R2" = "gz" ]; then

						rm "$file_R1"

					fi

					if [ "$file_type_R2" = "gz" ]; then

						rm "$file_R2"

					fi

					# end tracking of elapsed time #
					} 2>&1
					) 2>${OUTPREFIX}time_elapsed.log

					# store elapsed time
					readarray -t time_elapsed < ${OUTPREFIX}time_elapsed.log
					printf "$sample_ID,mitoz_bm,${time_elapsed[0]},${time_elapsed[1]},${time_elapsed[2]}\n" >> logs/${OUTPREFIX}assembly_time_elapsed_log.csv
					rm ${OUTPREFIX}time_elapsed.log

					# store end date and time #

					end_date=$(date +"%Y-%m-%d")
					end_time=$(date +"%H:%M:%S")

					# log start and end times #

					printf "$sample_ID,mitoz_bm,$start_date,$start_time,$end_date,$end_time\n" >> logs/${OUTPREFIX}assembly_times_log.csv
				# close if: deduptrim reads exist
				fi
			# close if: mitoz is true
			fi
			
###     NOVOPlasty assembly     ###################################################################

			if [ "$NOVOPLASTY" = "TRUE" ]; then

				## novoplasty ##

				# check existence of dependent files, if previous dependent step not failed log step failure
				
				if ! ls $PWD/deduptrimA | grep -q "$sample_ID"_R1.DRTRIMAP.fastq.gz && [ ! -n "$deduptrimA_failed" ]; then

					printf "$sample_ID,novoplasty_assembly\n" >> logs/${OUTPREFIX}failed_samples_log.csv
					printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
					echo -e "\nERROR: novoplasty assembly failed as deduplicated and adaptor trimmed $sample_ID files not found in /deduptrimA, skipping to next independent step\n"

				fi
		
				if ls $PWD/deduptrimA | grep -q "$sample_ID"_R1.DRTRIMAP.fastq.gz; then

					echo -e "\nAssembly method: novoplasty\n"


					# loop through concatenated vs non-concatenated reads #
					# default is concat=("no" "yes")
					for c in "${concat[@]}";
					do

						if [ $c = "no" ]; then

							# get appropriate file path #
							file_R1=$(find $PWD/deduptrimA/ | grep "$sample_ID"_R1.DRTRIMAP);
							file_R2=$(find $PWD/deduptrimA/ | grep "$sample_ID"_R2.DRTRIMAP);

						fi
						if [ $c = "yes" ]; then

							# get appropriate file path #
							file_R1=$(find $PWD/deduptrimACat/ | grep "$sample_ID"_R1.DRTRIMAC);
							file_R2=$(find $PWD/deduptrimACat/ | grep "$sample_ID"_R2.DRTRIMAC);

						fi
						# deduptrimA files do not need decompressing

						# loop through kmers #
						# default is kmer_values=(33 21);
						for k in "${kmer_values[@]}";
						do

							# loop through reference guide yes/no #
							# default is ref_guide=("no" "yes");
							for g in "${ref_guide[@]}";
							do

								if [ $g = "no" ]; then
									guide=""
									var_det="no"
								fi
								if [ $g = "yes" ]; then
									guide="$MT_REFERENCE"
									var_det="yes"
								fi

								echo -e "\nAssembly method options: $c concatenated reads\n"
								echo -e "\nAssembly method options: $k kmer\n"
								echo -e "\nAssembly method options: $g reference guide\n"


								# store start date and time #
								start_date=$(date +"%Y-%m-%d")
								start_time=$(date +"%H:%M:%S")

								# begin tracking of elapsed time #
								time (
								{

								# make directories #

								mkdir -p novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm;
								mkdir -p novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/config;
								mkdir -p novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID";

								# create config file #

								echo -e "Project:
-----------------------
Project name          = "$sample_ID"_assembly
Type                  = mito
Genome Range          = 16000-17000
K-mer                 = "$k"
Max memory            = "$MAX_MEMORY"
Extended log          = 0
Save assembled reads  = no
Seed Input            = "$MT_REFERENCE"
Extend seed directly  = no
Reference sequence    = "$guide"
Variance detection    = "$var_det"
Chloroplast sequence  = 

Dataset 1:
-----------------------
Read Length           = 101
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        = 
Forward reads         = "$file_R1"
Reverse reads         = "$file_R2"
Store Hash            =

Heteroplasmy:
-----------------------
MAF                   = 
HP exclude list       = 
PCR-free              = 

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = yes
Reduce ambigious N's  = 
Output path           = novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID"/" > $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/config/"$sample_ID"_config.txt

								# run novoplasty #

								perl $PWD/external_software/NOVOPlasty4.3.5.pl -c $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/config/"$sample_ID"_config.txt
			
								# check success
				
								if ! ls $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID" | grep -q "$sample_ID"_assembly.fasta$; then

									printf "$sample_ID,novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_assembly\n" >> logs/${OUTPREFIX}failed_samples_log.csv
									printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
									echo -e "\nERROR: novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat assembly failed, skipping to next independent step\n"
		
									# set fasta file name as failed for log
					
									fasta_file_name="assembly_failed"

									continue

								fi
				
								if ls $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID" | grep -q "$sample_ID"_assembly.fasta$; then

									echo -e "\nnovoplasty assembly concluded\n"

###################     NOVOPlasty assembly annotation via MitoZ     ##############################

									# annotation #

									# make subdirectory

									mkdir -p $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID"/annotation

									# get novoplasty output fasta file

									fasta_file=$(find $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID"/ | grep "$sample_ID"_assembly.fasta$)

									# store fasta file name for log
					
									fasta_file_name="${fasta_file##*/}"

									# rename contigs in file -> too long for mitoz function

									# get header line numbers and store as string
									header_lines=$(printf -- '%s\n' | grep -n -i ">" $fasta_file | cut -d : -f 1 | tr '\n' ' ')

									# convert string to array
									header_lines=($header_lines)

									# loop over contig header line array and replace contig names; sed may need sudo permissions
									for l in "${!header_lines[@]}" ;
									do
									
										sed -i "${header_lines[$l]} s/>.*/>contig_$l/" $fasta_file;
									
									done

									# check if original fasta_file is named circular and has only one contig to run mafft
									if echo $fasta_file_name | grep -q Circular && [ "${#header_lines[@]}" = 1 ]; then
										
										# run reorder nucleotides -> uses mafft to ensure nucleotides start at origin where possible
										bash auto_scripts/reorder_nucleotides.sh -f $fasta_file -r $MT_REFERENCE -o "$PWD/novoplasty_${k}kmer_${g}RefGuide_${c}Concat_bm/${sample_ID}" -c "yes"

										# replace fasta_file with reordered file if existing -> will not exist if reordering is not needed
										if [ -s "$PWD/novoplasty_${k}kmer_${g}RefGuide_${c}Concat_bm/${sample_ID}/realigned.fasta" ]; then
											# rename file
											mv $PWD/novoplasty_${k}kmer_${g}RefGuide_${c}Concat_bm/${sample_ID}/realigned.fasta $PWD/novoplasty_${k}kmer_${g}RefGuide_${c}Concat_bm/${sample_ID}/${sample_ID}_realigned.fasta
	
											# reassign fasta_file
											fasta_file=$PWD/novoplasty_${k}kmer_${g}RefGuide_${c}Concat_bm/${sample_ID}/${sample_ID}_realigned.fasta
										fi

										# check if max iterations reaches warning file exists and write to check warnings.csv
										if [ -s "$PWD/novoplasty_${k}kmer_${g}RefGuide_${c}Concat_bm/${sample_ID}/WARNING_reorder_nucleotides.txt" ]; then
											printf "$sample_ID,WARNING: reordering of nucleotides for ${fasta_file} reached maximum number of iterations, reordering may be incomplete\n" >> logs/${OUTPREFIX}check_warnings.csv
										fi
									
									fi

									#check if fasta file contains asterisk
									if grep -i -q "[*]" $fasta_file ; then
										echo -e "\nWARNING: ${fasta_file} contains asterisks, removing...\n"
										cp $fasta_file $PWD/novoplasty_${k}kmer_${g}RefGuide_${c}Concat_bm/${sample_ID}/${sample_ID}_contains_asterisks.fasta
										if [ $IGNORE_SUDO = "no" ]; then
											echo $password | sudo -S sed -i 's/[*]//g' $fasta_file;
										fi
										if [ $IGNORE_SUDO = "yes" ]; then
											sed -i 's/[*]//g' $fasta_file;
										fi
										printf "$sample_ID,WARNING: ${fasta_file} contained asterisks that were removed\n" >> logs/${OUTPREFIX}check_warnings.csv
									fi

									# check if original fasta_file is named circular and has only one contig to add circular information -> topology info must be added after mafft due to space character in header
									# check named circular
									if echo $fasta_file_name | grep -q Circular ; then
										# check only one contig
										if [ "${#header_lines[@]}" = 1 ]; then

											if [ $IGNORE_SUDO = "no" ]; then
												echo $password | sudo -S sed -i "${!header_lines[0]} s/>.*/>contig_0 topology=circular/" $fasta_file;
											fi
											if [ $IGNORE_SUDO = "yes" ]; then
												sed -i "${!header_lines[0]} s/>.*/>contig_0 topology=circular/" $fasta_file;
											fi									

										fi
									fi

									# confirm fasta_file name -> if fasta_file has been reassigned to realigned.fasta then name will be reassigned
									fasta_file_name="${fasta_file##*/}"

									# run mitoz
									if [ $UDOCKER = "no" ]; then
										docker run -v $PWD:$PWD -w $PWD --rm guanliangmeng/mitoz:3.6  mitoz annotate --workdir $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID"/annotation --fastafiles "$fasta_file" --fq1 "$file_R1" --fq2 "$file_R2" --clade Chordata --outprefix "$sample_ID" --species_name "$sample_ID"
									fi
									if [ $UDOCKER = "yes" ]; then
										udocker run  --volume=$PWD:$PWD --workdir=$PWD --rm mitoz_v3.6u:latest  mitoz annotate --workdir $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID"/annotation --fastafiles "$fasta_file" --fq1 "$file_R1" --fq2 "$file_R2" --clade Chordata --outprefix "$sample_ID" --species_name "$sample_ID"
									fi

									# check success

									if ! ls $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID"/annotation/"$sample_ID"."${fasta_file##*/}".result | grep -q summary.txt; then
										printf "$sample_ID,novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_mitoz_annotation\n" >> logs/${OUTPREFIX}failed_samples_log.csv
										printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
										echo -e "\nERROR: mitoz annotation of novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat failed\n"
										continue
									fi

									if ls $PWD/novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm/"$sample_ID"/annotation/"$sample_ID"."${fasta_file##*/}".result | grep -q summary.txt; then				

										echo -e "\nnovoplasty assembly mitoz annotation concluded\n"
		
										# store end date and time #
										end_date=$(date +"%Y-%m-%d")
										end_time=$(date +"%H:%M:%S")

										# log start and end times #
										printf "$sample_ID,novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm,$start_date,$start_time,$end_date,$end_time\n" >> logs/${OUTPREFIX}assembly_times_log.csv
					
									fi
								# close if: fasta exists
								fi
					
								# end tracking of elapsed time #
								} 2>&1
								) 2>${OUTPREFIX}time_elapsed.log

								# store elapsed time
								readarray -t time_elapsed < ${OUTPREFIX}time_elapsed.log
								printf "$sample_ID,novoplasty_"$k"kmer_"$g"RefGuide_"$c"Concat_bm,${time_elapsed[0]},${time_elapsed[1]},${time_elapsed[2]}\n" >> logs/${OUTPREFIX}assembly_time_elapsed_log.csv
								rm ${OUTPREFIX}time_elapsed.log

								# append annotated fasta file information to log
								#printf "${sample_ID},${k}kmer_${g}RefGuide_${c}Concat,${fasta_file_name}\n" >> novoplasty_fasta_file.csv
				
							# end ref_guide loop #
							done

						# end kmer loop #
						done

					# end concat_reads loop #
					done
				# close if: deduptrimA file exists
				fi
			# close if: novoplasty is true
			fi		
			
if [ "$MITOBIM" = "TRUE" ]; then

				## mitobim ##

				# check existence of dependent files, if previous dependent step not failed log step failure
				
				if ! ls $PWD/deduptrimCat | grep -q "$sample_ID".DRTRIMC.f && [ ! -n "$deduptrim_failed" ]; then

					printf "$sample_ID,mitobim_assembly\n" >> logs/${OUTPREFIX}failed_samples_log.csv
					printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
					echo -e "\nERROR: mitobim assembly failed as concatenated deduplicated and trimmed $sample_ID files not found in /deduptrimCat, skipping to next independent step\n"

				fi
		
				if ls $PWD/deduptrimCat | grep -q "$sample_ID".DRTRIMC.f ; then

					echo -e "\nAssembly method: mitobim\n"
			
					echo -e "\npath to mira binaries: $mirapath"

					export LC_ALL=C

					# print NFS warning status
			
					echo -e "stop mira if tmp directory is on NFS mount: $NFS_STOP"
			
					# get appropriate file path #
					reads=$(find $PWD/deduptrimCat/ | grep "$sample_ID".DRTRIMC.f)

					# store present working directory; needed to restore directory after changing to output directory
					if [ ! -n "$master_dir" ]; then
						master_dir=$PWD
					fi

					# loop through selected options; quick and two step
					# default is mitobim_options=( quick 2step )
					for o in "${mitobim_options[@]}";
					do

						echo -e "\nAssembly method options: $o\n"

						# check directory as changes occur in loop
						if [ -n "$master_dir" ]; then
							if [ ! "$master_dir" = "$PWD" ]; then
								cd "${master_dir}"
								echo -e "\n working directory is ${PWD}\n"
							fi
						fi

						# make directories
						mkdir -p "$PWD/mitobim_${o}_bm"

						# store start date and time #
						start_date=$(date +"%Y-%m-%d")
						start_time=$(date +"%H:%M:%S")

						# begin tracking of elapsed time #
						time (
						{

						if [ "$o" = quick ]; then
					
							## quick ##

							# create and change to output directory
							mkdir -p $PWD/mitobim_${o}_bm/$sample_ID
							cd $PWD/mitobim_${o}_bm/$sample_ID 

							echo -e "\n working directory is ${PWD}\n"

							echo -e "\ngetting reads from ${reads}\n"
							
							# set strain_name to sample_ID
							strain_name=${sample_ID}
							
							# check if sample_ID contains period
							if echo ${sample_ID} | grep "\." ; then

								# replace . character in sample_ID with underscore and write to strain name
								strain_name=`echo "${sample_ID/\./_}"`

								# replace sample_ID with new strain name in .maf file
								# sed -i "s/${sample_ID}/${strain_name}/g" "${sample_ID}"_mira_assembly/"${sample_ID}"_mira_d_results/"${sample_ID}"_mira_out.maf

							fi

							if [ $NFS_STOP = "no" ]; then

								# run mitobim
								"$master_dir"/external_software/MITObim/MITObim.pl -start 0 -end 100 -sample "$strain_name" -ref "${MT_REFERENCE##*/}" -readpool "$reads" -quick "$MT_REFERENCE" --clean --NFS_warn_only

							fi

							if [ $NFS_STOP = "yes" ]; then

								# run mitobim
								"$master_dir"/external_software/MITObim/MITObim.pl -start 0 -end 100 -sample "$strain_name" -ref "${MT_REFERENCE##*/}" -readpool "$reads" -quick "$MT_REFERENCE" --clean

							fi
					
							# restore master working directory
							cd "${master_dir}"

							echo -e "\n working directory is ${PWD}\n"
			
						fi

						if [ "$o" = 2step ]; then

							## two-step ##

							# make manifest directory

							mkdir -p $PWD/mitobim_2step_bm/manifest

							# copy mt_reference to file with extension .fna #
					
							if [ ! -f "${MT_REFERENCE%.*}.fna" ]; then
								# store original file extension
								mt_extension="${MT_REFERENCE##*.}"
					
								# edit file extention
								cp ${MT_REFERENCE} "${MT_REFERENCE%.*}.fna"
							fi
	
							# overwrite MT_REFERENCE with new extension
							MT_REFERENCE="${MT_REFERENCE%.*}.fna"

							# get mt_reference name without extension
							mt_sans_path="${MT_REFERENCE##*/}"
							mt_name_only=${mt_sans_path%.*}
			
							if [ $NFS_STOP = "no" ]; then

								# make mira manifest
								echo -e "\n#manifest file for basic mapping assembly with illumina data using MIRA 4\n\nproject = "$sample_ID"_mira\n\njob=genome,mapping,accurate\n\nparameters = -NW:mrnl=0:cnfs=warn:cac=warn -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no\n\nreadgroup\nis_reference\ndata = "$MT_REFERENCE"\nstrain = "${mt_name_only}"\n\nreadgroup = reads\ndata = "$reads"\ntechnology = solexa\nstrain = "$sample_ID"\n" > $PWD/mitobim_2step_bm/manifest/"$sample_ID".conf

							fi

							if [ $NFS_STOP = "yes" ]; then

								# make mira manifest
								echo -e "\n#manifest file for basic mapping assembly with illumina data using MIRA 4\n\nproject = "$sample_ID"_mira\n\njob=genome,mapping,accurate\n\nparameters = -NW:mrnl=0:cac=warn -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no\n\nreadgroup\nis_reference\ndata = "$MT_REFERENCE"\nstrain = "${mt_name_only}"\n\nreadgroup = reads\ndata = "$reads"\ntechnology = solexa\nstrain = "$sample_ID"\n" > $PWD/mitobim_2step_bm/manifest/"$sample_ID".conf

							fi
					
							# create and change to output directory
							mkdir -p $PWD/mitobim_2step_bm/$sample_ID
							cd $PWD/mitobim_2step_bm/$sample_ID

							echo -e "\n working directory is ${PWD}\n"

							echo -e "\ngetting reads from ${reads}\n"		
								
							# run initial mira step #
							if [ "$RESUME_MITOBIM" == "no" ]; then
								mira "$master_dir"/mitobim_2step_bm/manifest/"$sample_ID".conf
							fi

							# check for success
							if ! ls "${master_dir}"/mitobim_2step_bm/${sample_ID}/${sample_ID}_mira_assembly/*results | grep -q .maf$ ; then
								printf "$sample_ID,2step_mitobim_mira_assembly\n" >> logs/${OUTPREFIX}failed_samples_log.csv
								printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
								cd "${master_dir}"
								echo -e "\nERROR: mira assembly of 2step mitobim failed\n"
								continue
							fi

							if ls "${master_dir}"/mitobim_2step_bm/${sample_ID}/${sample_ID}_mira_assembly/*results | grep -q .maf$ ; then

								# set strain_name to sample_ID
								strain_name=${sample_ID}

								# check if sample_ID contains period
								if echo ${sample_ID} | grep "\." ; then
									
									# replace . character in sample_ID with underscore and write to strain name
									strain_name=`echo "${sample_ID//\./_}"`

									# replace sample_ID with new strain name in .maf file
									sed -i "s/${sample_ID}/${strain_name}/g" "${sample_ID}"_mira_assembly/"${sample_ID}"_mira_d_results/"${sample_ID}"_mira_out.maf

								fi

								# run mitobim #

								if [ $NFS_STOP = "no" ]; then

									# run mitobim
									if [ "$RESUME_MITOBIM" == "no" ]; then
										"$master_dir"/external_software/MITObim/MITObim.pl -start 1 -end 100 -sample $strain_name -ref "${mt_name_only}" -readpool "$reads" -maf "${sample_ID}"_mira_assembly/"${sample_ID}"_mira_d_results/"${sample_ID}"_mira_out.maf --clean --NFS_warn_only
									else

										# locate last successful iteration folder - will be second folder with clean when sorted numerically from first to last
										last_successful_iteration=`ls | grep iteration | sort -V | tail -n +2 | head -n1`
										# remove unsuccessful iteration
										rm -r `ls | grep iteration | sort -V | tail -n 1`
										# get .maf file from last iteration
										last_successful_maf=`ls ${last_successful_iteration}/*_assembly/*_results/*.maf`
										# get new start number
										start_from="$(( ${last_successful_iteration#iteration} + 1 ))"

										# run mitobim
										"$master_dir"/external_software/MITObim/MITObim.pl -start ${start_from} -end 100 -sample $strain_name -ref "${mt_name_only}" -readpool "$reads" -maf ${last_successful_maf} --clean --NFS_warn_only
									fi

								fi

								if [ $NFS_STOP = "yes" ]; then

									if [ "$RESUME_MITOBIM" == "no" ]; then

										# run mitobim
										"$master_dir"/external_software/MITObim/MITObim.pl -start 1 -end 100 -sample $strain_name -ref "${mt_name_only}" -readpool "$reads" -maf "${sample_ID}"_mira_assembly/"${sample_ID}"_mira_d_results/"${sample_ID}"_mira_out.maf --clean
									else

										# locate last successful iteration folder - will be second folder with clean when sorted numerically from first to last                                                                                                                                                                                 
										last_successful_iteration=`ls | grep iteration | sort -V | tail -n +2 | head -n1`
										# remove unsuccessful iteration
                                                                                rm -r `ls | grep iteration | sort -V | tail -n 1`
                                                                                # get .maf file from last iteration
                                                                                last_successful_maf=`ls ${last_successful_iteration}/*_assembly/*_results/*.maf`
                                                                                # get new start number
                                                                                start_from="$(( ${last_successful_iteration#iteration} + 1 ))"

                                                                                # run mitobim
                                                                                "$master_dir"/external_software/MITObim/MITObim.pl -start ${start_from} -end 100 -sample $strain_name -ref "${mt_name_only}" -readpool "$reads" -maf ${last_successful_maf} --clean
									fi

								fi

								# restore master working directory
								cd "${master_dir}"

								echo -e "\n working directory is ${PWD}\n"

								# restore MT_REFERENCE to file with original extension #

								# overwrite MT_REFERENCE with old extension
								MT_REFERENCE="${MT_REFERENCE%.*}.${mt_extension}"

							fi
						# close if: option is 2step
						fi

						# confirm master working directory
						if [ ! "$master_dir" = $"PWD" ]; then
							cd "${master_dir}"
							echo -e "\n working directory is ${PWD}\n"
						fi

						# get last iteration folder
		
						iteration_folder=`ls $PWD/mitobim_${o}_bm/${sample_ID} | grep -i iteration* | sort -V | tail -n1`

						# check success -> look for out fasta file #

						if ! ls $PWD/mitobim_${o}_bm/${sample_ID}/${iteration_folder} | grep -q noIUPAC.fasta$ ; then
							printf "$sample_ID,${o}_mitobim_assembly\n" >> logs/${OUTPREFIX}failed_samples_log.csv
							printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
							echo -e "\nERROR: assembly via ${o} mitobim failed\n"
						fi

						if ls $PWD/mitobim_${o}_bm/${sample_ID}/${iteration_folder} | grep -q noIUPAC.fasta$ ; then
							
							echo -e "\n${o} mitobim assembly concluded\n"

###################################     MITObim assembly annotation via MitoZ     ##############

							## annotation via mitoz ##

							# make subdirectory

							mkdir -p $PWD/mitobim_${o}_bm/"$sample_ID"/mt_annotation

							# get mitobim out fasta

							mitobim_out=$(find $PWD/mitobim_${o}_bm/${sample_ID}/${iteration_folder} | grep -i noIUPAC.fasta$)

							# rename contigs in file -> too long for mitoz function

							# get header line numbers and store as string -> should only be one
							header_lines=$(printf -- '%s\n' | grep -n -i ">" $mitobim_out | cut -d : -f 1 | tr '\n' ' ')

							# convert string to array
							header_lines=($header_lines)

							# loop over contig header line array and replace contig names; sed may need sudo permissions
							for l in "${!header_lines[@]}" ;
							do
								sed -i "${header_lines[$l]} s/>.*/>mito_$l/" $mitobim_out

							done

							# run circle_check on mitobim out #

							# make subdirectory
							mkdir -p $PWD/mitobim_${o}_bm/${sample_ID}/circle_check

							# run circle check
							udocker run  --volume=$PWD:$PWD --workdir=$PWD --rm mitoz_v3.6u:latest  mitoz-tools circle_check -i "${mitobim_out}" -o "$PWD/mitobim_${o}_bm/${sample_ID}/circle_check/${sample_ID}"
					
							# if overlap information is empty then no overlap found, leave mitobim to be assumed linear

							# if overlap information is not empty then overlap found -> overwrite mitobim_out with circle_check output
							if [ -s $PWD/mitobim_${o}_bm/${sample_ID}/circle_check/${sample_ID}.overlap_information ]; then
								mitobim_out=$PWD/mitobim_${o}_bm/${sample_ID}/circle_check/${sample_ID}.mitogenome.fa
							fi

							# run mitoz
							if [ $UDOCKER = "no" ]; then
								docker run -v $PWD:$PWD -w $PWD --rm guanliangmeng/mitoz:3.6  mitoz annotate --workdir $PWD/mitobim_${o}_bm/"$sample_ID"/mt_annotation --fastafiles "$mitobim_out" --fq1 "$reads" --clade Chordata --outprefix "$sample_ID" --species_name "$sample_ID"
							fi
							if [ $UDOCKER = "yes" ]; then
								udocker run  --volume=$PWD:$PWD --workdir=$PWD --rm mitoz_v3.6u:latest  mitoz annotate --workdir $PWD/mitobim_${o}_bm/"$sample_ID"/mt_annotation --fastafiles "$mitobim_out" --fq1 "$reads" --clade Chordata --outprefix "$sample_ID" --species_name "$sample_ID"
							fi

							# check success
				
							if ! ls $PWD/mitobim_${o}_bm/"$sample_ID"/mt_annotation/"$sample_ID"."${mitobim_out##*/}".result | grep -q summary.txt; then

								printf "$sample_ID,mitobim_${o}_mitoz_annotation\n" >> logs/${OUTPREFIX}failed_samples_log.csv
								printf "$sample_ID\n" >> logs/${OUTPREFIX}failed_samples.txt
								echo -e "\nERROR: mitoz annotation of ${o} mitobim failed\n"

							fi

							if ls $PWD/mitobim_${o}_bm/"$sample_ID"/mt_annotation/"$sample_ID"."${mitobim_out##*/}".result | grep -q summary.txt; then				

								echo -e "\n${o} mitobim assembly mitoz annotation concluded\n"
		
								# store end date and time #
								end_date=$(date +"%Y-%m-%d")
								end_time=$(date +"%H:%M:%S")

								# log start and end times #
								printf "$sample_ID,mitobim_${o}_bm,$start_date,$start_time,$end_date,$end_time\n" >> logs/${OUTPREFIX}assembly_times_log.csv
					
							fi
						# close if: fasta exists
						fi

						# end tracking of elapsed time #
						} 2>&1
						) 2>${OUTPREFIX}time_elapsed.log

						# store elapsed time
						readarray -t time_elapsed < ${OUTPREFIX}time_elapsed.log
						printf "$sample_ID,mitobim_${o}_bm,${time_elapsed[0]},${time_elapsed[1]},${time_elapsed[2]}\n" >> logs/${OUTPREFIX}assembly_time_elapsed_log.csv
						rm ${OUTPREFIX}time_elapsed.log
			
					# end mitobim options loop #
					done

					# REMOVE DEDUPTRIMCAT
					rm $reads

				# close if: deduptrimCat reads exist
				fi
			# close if: mitobim is true
			fi

			# confirm master working directory
			if [ ! "$master_dir" = $"PWD" ]; then
				cd "${master_dir}"
				echo -e "\n working directory is ${PWD}\n"
			fi

# RUN TEST: mitoz assembly
fi # close if: method is not prep reads
done # close loop: sample ids
fi # close id: method is NOT none
#exit 0
# END TEST: mitoz assembly

###################################################################################################
#                                                                                                 #
# COMPILATION OF BENCHMARK STATISTICS                                                             #
#                                                                                                 #
# uses Rscripts to compile whole genome and individual gene stats for all _bm folders             #
#                                                                                                 #
###################################################################################################

if [ $ASSEMBLY_STATS == "yes" ]; then

        ## benchmark stats ##

        # make tmp directory
        mkdir -p ${OUTPREFIX}tmp

        # retrieve sample_IDs

        if [ $OVERRIDE_RAW = "no" ]; then

                if [ ! -d "${RAW_READS}" ] ; then

                        echo -r "${RAW_READS} does not exist. Supply path to raw reads to retrieve sample_IDs or supply sample_ID .txt file with --sample_id_list."
                        exit 1
                fi

                echo -e "\nretrieving sample_IDs from raw file names in ${RAW_READS}\n"

                # store raw_file names
                file_names=( $(ls ${RAW_READS}) )

                # get R1 file names and store as string; add .[IDENTIFIER] to _R1 to access certain type of reads e.g. trimmed vs untrimmed; paired vs unpaired
                R1_files=$(printf -- '%s\n' "${file_names[@]}" | grep _R1 | tr '\n' ' ')

                # convert string to array
                sample_array=($R1_files)

                # remove any existing sampleID file in tmp folder
                if [ -f ${OUTPREFIX}tmp/sampleIDs.txt ] ; then
                        rm ${OUTPREFIX}tmp/sampleIDs.txt
                fi

                # loop over sample array and write sample_IDs to txt file
                for i in "${sample_array[@]}";
                do
                        sample_ID="${i%_R1.*}"
                        echo -e "${sample_ID}" >> ${OUTPREFIX}tmp/sampleIDs.txt
                done

        fi

        if [ $OVERRIDE_RAW = "yes" ]; then

                echo -e "\nretrieving sample_IDs from supplied text file\n"

                # cp supplied text file to tmp directory
                cp "$SAMPLE_IDS" ${OUTPREFIX}tmp/sampleIDs.txt

        fi

        # sample_ID list is ${OUTPREFIX}tmp/sampleIDs.txt

        # remove any existing methods file in tmp folder
        if [ -f ${OUTPREFIX}tmp/benchmark_methods.txt ] ; then
                rm ${OUTPREFIX}tmp/benchmark_methods.txt
        fi

        # create list of _bm folders -> from folders in directory
        if [ ! -n $METHODS_LIST ]; then
                bm=`ls ${OUTDIR} | grep _bm`
                bm=( ${bm} )

                for b in "${bm[@]}" ;
                do
                        echo -e "${b}" >> ${OUTPREFIX}tmp/benchmark_methods.txt
                done
        fi

        # create list of _bm folders -> from given list
        if [ -n $METHODS_LIST ]; then
                cp ${METHODS_LIST} ${OUTPREFIX}tmp/benchmark_methods.txt
        fi


	# if expected_genes file exists
	if [ -s "${EXPECTED_GENES}" ]; then

	        # copy expected gene list and number of copies
        	cp ${EXPECTED_GENES} ${OUTPREFIX}tmp/expected_genes.csv
	else
			# create expected gene list and number of copies
	echo -e "gene,num_copies
ATP6,1
ATP8,1
COX1,1
COX2,1
COX3,1
CYTB,1
ND1,1
ND2,1
ND3,1
ND4,1
ND4L,1
ND5,1
ND6,1
l-rRNA,1
s-rRNA,1
trnA,1
trnR,1
trnN,1
trnD,1
trnC,1
trnQ,1
trnE,1
trnG,1
trnH,1
trnI,1
trnL,2
trnK,1
trnM,1
trnF,1
trnP,1
trnS,2
trnT,1
trnW,1
trnY,1
trnV,1" > ${OUTPREFIX}tmp/expected_genes.csv

	fi

	echo ${EXPECTED_GENES}

        # create result .csv files

        mkdir -p $PWD/benchmark_stats

        # create stats .csv files with header rows

        if [ ! -n "${EXCLUDE_GENES_READ_DEPTH}" ]; then

                printf "sample_ID,method,num_Ns,num_bp,num_bp_minus_Ns,Ns_per_bp,bp_per_refbp,bp_minus_Ns_per_refbp,num_genes,num_CDS_genes,num_rRNA_genes,num_tRNA_genes,circularity,num_contigs,potential_contamination,missing_stop_codons,num_internal_stops,avg_read_depth,med_read_depth,min_read_depth,max_read_depth\n" > benchmark_stats/${OUTPREFIX}whole_genome_criteria.csv

        fi

        if [ -n "${EXCLUDE_GENES_READ_DEPTH}" ]; then

                printf "sample_ID,method,num_Ns,num_bp,num_bp_minus_Ns,Ns_per_bp,bp_per_refbp,bp_minus_Ns_per_refbp,num_genes,num_CDS_genes,num_rRNA_genes,num_tRNA_genes,circularity,num_contigs,potential_contamination,missing_stop_codons,num_internal_stops,avg_read_depth,med_read_depth,min_read_depth,max_read_depth,avg_read_depth_excluded_genes,med_read_depth_excluded_genes,min_read_depth_excluded_genes,max_read_depth_excluded_genes\n" > benchmark_stats/${OUTPREFIX}whole_genome_criteria.csv

        fi

        printf "sample_ID,method,gene,num_Ns,num_bp,num_bp_minus_Ns,Ns_per_bp,bp_per_refbp,bp_minus_Ns_per_refbp,stop_codon_detected,num_copies_detected,num_internal_stops,avg_read_depth,med_read_depth,min_read_depth,max_read_depth,contig_name,contig_number,contig_length\n" > benchmark_stats/${OUTPREFIX}gene_criteria.csv

        # store present working directory
        master_dir="${PWD}"

        # make tmp logs directory -> stores parallel log files
        mkdir -p ${PWD}/${OUTPREFIX}tmp_logs/

        if [ $SKIP_FILE_COMPILATION = "no" ]; then

                if [ ! "${OUTDIR}" = "$PWD" ]; then

                        echo -e "changing working directory to ${OUTDIR} to retrieve file paths"
                        cd ${OUTDIR}

                fi

                echo -e "\nretrieving necessary file paths for benchmark statistics\n"

                # get file paths
		fasta_paths_tmp=`mktemp ${master_dir}/${OUTPREFIX}tmp_logs/find_fasta_XXXXX.txt`
		fasta_tmp_paths_tmp=`mktemp ${master_dir}/${OUTPREFIX}tmp_logs/find_tmp_fasta_XXXXX.txt`
		annotated_fasta_paths_tmp=`mktemp ${master_dir}/${OUTPREFIX}tmp_logs/find_annotated_fasta_XXXXX.txt`
		summary_paths_tmp=`mktemp ${master_dir}/${OUTPREFIX}tmp_logs/find_summary_XXXXX.txt`
		read_depth_paths_tmp=`mktemp ${master_dir}/${OUTPREFIX}tmp_logs/find_depth_XXXXX.txt`
		translated_fasta_paths_tmp=`mktemp ${master_dir}/${OUTPREFIX}tmp_logs/find_translated_fasta_XXXXX.txt`

		( find -maxdepth 5 -name "*_mitoscaf.fa.gbf.fasta" | grep -v "/tmp" > $fasta_paths_tmp ) &
		( find -maxdepth 5 -wholename "*_mitoscaf.fa/*_mitoscaf.fa" | grep "/tmp" > $fasta_tmp_paths_tmp ) &
		( find -maxdepth 5 -name "*_mitoscaf.fa.gbf.gene.fasta" | grep -v "/tmp" > $annotated_fasta_paths_tmp ) &
		( find -maxdepth 5 -name "summary.txt" | grep -v "/tmp" > $summary_paths_tmp ) &
		( find -maxdepth 5 -name "circos.depth.txt" > $read_depth_paths_tmp ) &
		( find -maxdepth 5 -name "*_mitoscaf.fa.gbf.cds_translation.fasta" | grep -v "/tmp" > $translated_fasta_paths_tmp ) &

		# wait for all paths to be found
		wait

		# read in paths
		readarray -t fasta_paths < $fasta_paths_tmp
		readarray -t fasta_tmp_paths < $fasta_tmp_paths_tmp
		readarray -t annotated_fasta_paths < $annotated_fasta_paths_tmp
		readarray -t summary_paths < $summary_paths_tmp
		readarray -t read_depth_paths < $read_depth_paths_tmp
		readarray -t translated_fasta_paths < $translated_fasta_paths_tmp

		# get fasta file paths
		#fasta_paths=( `find -maxdepth 5 -name "*_mitoscaf.fa.gbf.fasta" | grep -v "/tmp"` )

                # get fasta file tmp paths  -> if annotation is unsuccessful input fasta should still be available from mitoz tmp
		#fasta_tmp_paths=( `find -maxdepth 5 -wholename "*_mitoscaf.fa/*_mitoscaf.fa" | grep "/tmp"` )

                # get annotated fasta file paths
		#annotated_fasta_paths=( `find -maxdepth 5 -name "*_mitoscaf.fa.gbf.gene.fasta" | grep -v "/tmp"` )

                # get summary text file paths
		#summary_paths=( `find -maxdepth 5 -name "summary.txt" | grep -v "/tmp"` )

                # get read_depth text file paths
		#read_depth_paths=( `find -maxdepth 5 -name "circos.depth.txt"` )

                # get translated fasta file paths
		#translated_fasta_paths=( `find -maxdepth 5 -name "*_mitoscaf.fa.gbf.cds_translation.fasta" | grep -v "/tmp"` )

                if [ ! "${master_dir}" = "$PWD" ]; then

                        echo -e "restoring working directory to ${master_dir}"
                        cd ${master_dir}

                fi

                # create copies of necessary files with consistent and specific naming and paths

                readarray -t sample_IDs < ${PWD}/${OUTPREFIX}tmp/sampleIDs.txt
                readarray -t methods < ${PWD}/${OUTPREFIX}tmp/benchmark_methods.txt

                # create function to find and move files

                find_copy_files() {

                        # get inputs
                        #s=$1
                        #methods=$2
                        #fasta_paths=$3
                        #fasta_tmp_paths=$4
                        #annotated_fasta_paths=$5
                        #summary_paths=$6
                        #read_depth_paths=$7
                        #translated_fasta_paths=$8

                        echo -e "\n${s}\n"

                        for m in "${methods[@]}" ;
                        do

                                echo -e "\n${m}\n"

                                # create empty path and name variables -> prevents interference from previous loop iterations
                                fasta_path=""
                                annotated_fasta_path=""
                                summary_path=""
                                read_depth_path=""
                                fasta_name=""
                                annotated_fasta_name=""
                                summary_name=""
                                read_depth_name=""
                                translated_fasta_path=""
                                translated_fasta_name=""

                                # get sample and method specific fasta file path via null separated grep
				fasta_path=`printf '%s\0' "${fasta_paths[@]}" | grep -F -z "_mitoscaf.fa.gbf.fasta" | grep -z "/${s}/" | grep -z "/${m}/" | tr '\0' '\n'`
				
				# check for path and if found look for further results files
				if [ ! "${fasta_path}" = "" ]; then
					
					echo -e "found fasta path: ${fasta_path}"
					fasta_name=$(basename ${fasta_path})

					annotated_fasta_path=`printf '%s\0' "${annotated_fasta_paths[@]}" | grep -F -z "_mitoscaf.fa.gbf.gene.fasta" | grep -z "/${s}/" | grep -z "/${m}/" | tr '\0' '\n'`

					if [ ! "${annotated_fasta_path}" = "" ]; then
						echo -e "found annotated fasta path: ${annotated_fasta_path}"
                                                annotated_fasta_name=$(basename ${annotated_fasta_path})
					else
						echo -e "WARNING: no annotated fasta file found for ${s} in ${m} in mitoz directory, check manually for file\n"
					fi

					summary_path=`printf '%s\0' "${summary_paths[@]}" | grep -F -z "summary.txt" | grep -z "/${s}/" | grep -z "/${m}/" | tr '\0' '\n'`

					if [ ! "${summary_path}" = "" ]; then
						echo -e "found summary path: ${summary_path}"
						summary_name=$(basename ${summary_path})
					else
						echo -e "WARNING: no summary file found for ${s} in ${m} in mitoz directory, check manually for file\n"
					fi

					read_depth_path=`printf '%s\0' "${read_depth_paths[@]}" | grep -F -z "circos.depth.txt" | grep -z "/${s}/" | grep -z "/${m}/" | tr '\0' '\n'`

                                        if [ ! "${read_depth_path}" = "" ]; then
                                                echo -e "found read depth path: ${read_depth_path}"
                                                read_depth_name=$(basename ${read_depth_path})
                                        else
                                                echo -e "WARNING: no read depths file found for ${s} in ${m} in mitoz directory, check manually for file\n"
                                        fi

					translated_fasta_path=`printf '%s\0' "${translated_fasta_paths[@]}" | grep -F -z "_mitoscaf.fa.gbf.cds_translation.fasta" | grep -z "/${s}/" | grep -z "/${m}/" | tr '\0' '\n'`

                                        if [ ! "${translated_fasta_path}" = "" ]; then
                                                echo -e "found translated fasta path: ${translated_fasta_path}"
                                                translated_fasta_name=$(basename ${translated_fasta_path})
                                        else
                                                echo -e "WARNING: no translated fasta file found for ${s} in ${m} in mitoz directory, check manually for file\n"
                                        fi

				# else if fasta file not found in results directory, look in tmp
				else

                                        echo -e "\nWARNING: no input fasta file found for ${s} in ${m} now assuming mitoz annotation was unsuccessful\nChecking in mitoz tmp directory for input fasta\n"
					fasta_path=`printf '%s\0' "${fasta_tmp_paths[@]}" | grep -z "_mitoscaf.fa$" | grep -z "/${s}/" | grep -z "/${m}/" | tr '\0' '\n'`

	                                if [ "${fasta_path}" = "" ]; then
        	                                echo -e "WARNING: no input fasta file found for ${s} in ${m} in mitoz tmp directory now assuming annotation or assembly failed\n"
					else
						echo -e "found fasta path: ${fasta_path}"
						fasta_name=$(basename ${fasta_path})
					fi
				fi

                        if [ ! "${fasta_path}" = "" ]; then

                                # make unique directory to store copied files
                                mkdir -p ${PWD}/${OUTPREFIX}tmp/${s}/${m}

                                echo -e "\ncopying ${s} ${m} files into ${PWD}/${OUTPREFIX}tmp/${s}/${m}/"

                                # copy files
                                cp ${fasta_path} ${annotated_fasta_path} ${summary_path} ${read_depth_path} ${translated_fasta_path} ${PWD}/${OUTPREFIX}tmp/${s}/${m}/

                                # move and rename copied files

                                if [ ! "${fasta_name}" = "" ]; then
                                        echo -e "\nmoving ${PWD}/${OUTPREFIX}tmp/${s}/${m}/${fasta_name} to ${PWD}/${OUTPREFIX}tmp/${s}.${m}.input.fasta\n"
                                        mv "${PWD}/${OUTPREFIX}tmp/${s}/${m}/${fasta_name}" "${PWD}/${OUTPREFIX}tmp/${s}.${m}.input.fasta"
                                fi
                                if [ ! "${annotated_fasta_name}" = "" ]; then
                                        echo -e "\nmoving ${PWD}/${OUTPREFIX}tmp/${s}/${m}/${annotated_fasta_name} to ${PWD}/${OUTPREFIX}tmp/${s}.${m}.annotated.fasta\n"
                                        mv "${PWD}/${OUTPREFIX}tmp/${s}/${m}/${annotated_fasta_name}" "${PWD}/${OUTPREFIX}tmp/${s}.${m}.annotated.fasta"
                                fi
                                if [ ! "${summary_name}" = "" ]; then
                                        echo -e "\nmoving ${PWD}/${OUTPREFIX}tmp/${s}/${m}/${summary_name} to ${PWD}/${OUTPREFIX}tmp/${s}.${m}.summary.txt\n"
                                        mv "${PWD}/${OUTPREFIX}tmp/${s}/${m}/${summary_name}" "${PWD}/${OUTPREFIX}tmp/${s}.${m}.summary.txt"
                                fi
                                if [ ! "${read_depth_name}" = "" ]; then
                                        echo -e "\nmoving ${PWD}/${OUTPREFIX}tmp/${s}/${m}/${read_depth_name} to ${PWD}/${OUTPREFIX}tmp/${s}.${m}.read_depth.txt\n"
                                        mv "${PWD}/${OUTPREFIX}tmp/${s}/${m}/${read_depth_name}" "${PWD}/${OUTPREFIX}tmp/${s}.${m}.read_depth.txt"
                                fi
                                if [ ! "${translated_fasta_name}" = "" ]; then
                                        echo -e "\nmoving ${PWD}/${OUTPREFIX}tmp/${s}/${m}/${translated_fasta_name} to ${PWD}/${OUTPREFIX}tmp/${s}.${m}.translated.fasta\n"
                                        mv "${PWD}/${OUTPREFIX}tmp/${s}/${m}/${translated_fasta_name}" "${PWD}/${OUTPREFIX}tmp/${s}.${m}.translated.fasta"
                                fi

                        fi
                        # close methods loop
                        done

                        if [ -d ${PWD}/${OUTPREFIX}tmp/${s} ]; then
                                # remove empty ${s} directory
                                rm -r ${PWD}/${OUTPREFIX}tmp/${s}
                        fi
                # close function
                }


                # call function to find and move files

                for s in "${sample_IDs[@]}" ;
                do
			if [[ $(jobs -r -p | wc -l) -ge $num_processes ]]; then wait -n ; fi
                        #find_copy_files ${s} ${methods[@]} ${fasta_paths[@]} ${fasta_tmp_paths[@]} ${annotated_fasta_paths[@]} ${summary_paths[@]} ${read_depth_paths[@]} ${translated_fasta_paths[@]} > ${PWD}/${OUTPREFIX}tmp_logs/${s}.find_copy_files_log 2>&1 &
                        find_copy_files > ${PWD}/${OUTPREFIX}tmp_logs/${s}.find_copy_files_log 2>&1 &

                        # close sample loop
                done
                # wait for all child processes to finish
                wait

                # write logs to standard out
                cat ${PWD}/${OUTPREFIX}tmp_logs/*.find_copy_files_log

fi

        # create exclude genes from read depths file

        if [ -n "${EXCLUDE_GENES_READ_DEPTH}" ]; then

                # convert to array
                EXCLUDE_GENES_READ_DEPTH=( ${EXCLUDE_GENES_READ_DEPTH} )

                if [ -f ${PWD}/${OUTPREFIX}tmp/exclude_genes_read_depths.txt ]; then
                        rm "${PWD}/${OUTPREFIX}tmp/exclude_genes_read_depths.txt"
                fi

                for o in "${EXCLUDE_GENES_READ_DEPTH[@]}";
                do

                        printf "${o}\n" >> "${PWD}/${OUTPREFIX}tmp/exclude_genes_read_depths.txt"

                done

        fi


        #RSCRIPTS

        # export OUTPREFIX value for retrieval by Rscript
        export OUTPREFIX_VALUE_EXPORT=${OUTPREFIX}

        echo -e "\nsplitting summary.txt files into .csv files\n"

        # split summary.txt files and write to .csv files
		Rscript auto_scripts/split_summary_csv.R

        # export MT_REFERENCE value for retrieval by Rscript
        export MT_REFERENCE_VALUE_EXPORT=${MT_REFERENCE}

        echo -e "\ncompiling whole genome stats\n"

        # retrieve, calculate, and compile whole genome stats and write to .csv file
        Rscript auto_scripts/whole_genome_stats.R

        echo -e "\ncompiled whole genome stats found at $PWD/benchmark_stats/${OUTPREFIX}whole_genome_criteria.csv\n"

        # export MT_REFERENCE_ANNOTATED value for retrieval by Rscript
        export MT_REFERENCE_ANNOTATED_VALUE_EXPORT=${MT_REFERENCE_ANNOTATED}

        echo -e "\ncompiling individual gene stats\n"

        # retrieve, calculate, and compile individual gene stats and write to .csv file
        # Rscript auto_scripts/individual_gene_stats.R

        #conda activate mitogenome_stats_R
    Rscript auto_scripts/individual_gene_stats.R
  #conda deactivate

        echo -e "\ncompiled individual gene stats found at $PWD/benchmark_stats/${OUTPREFIX}gene_criteria.csv\n"

        echo -e "\ncompiling elapsed time logs\n"

        # get assembly_time_elapsed_log.csv files -> matching provided log pattern
        assembly_time_logs=`find -name "*assembly_time_elapsed_log.csv" | grep "${LOG_PATTERN}" | grep -v "${LOG_PATTERN_EXCLUDE}"`

        # convert to array
        assembly_time_logs=( $assembly_time_logs )

        # concatenate files
        cat ${assembly_time_logs[@]} > ${PWD}/${OUTPREFIX}tmp/elapsed_time_logs.csv

        # remove extra headers -> recreating file instead of inline editing voids sudo permissions
        sed '1!{/sample_ID,assembly_method,real,user,system/d}' ${PWD}/${OUTPREFIX}tmp/elapsed_time_logs.csv > ${PWD}/${OUTPREFIX}tmp/elapsed_time_logs1.csv

        # rename elapsed_time_logs1.csv to elapsed_time_logs.csv
        mv ${PWD}/${OUTPREFIX}tmp/elapsed_time_logs1.csv ${PWD}/${OUTPREFIX}tmp/elapsed_time_logs.csv

        # rename elapsed_time_logs.csv to CPU_time.csv
        mv ${PWD}/${OUTPREFIX}tmp/elapsed_time_logs.csv ${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}CPU_time.csv

        # move file to benchmark stats
        mv ${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}CPU_time.csv ${PWD}/benchmark_stats/${OUTPREFIX}CPU_time.csv


  # re-reading in sample and methods arrays because skipped when skipping file compilation
  readarray -t sample_IDs < ${PWD}/${OUTPREFIX}tmp/sampleIDs.txt
  readarray -t methods < ${PWD}/${OUTPREFIX}tmp/benchmark_methods.txt

  echo -e "\ncombining sequences into single fasta file\ncombining..."

        # create function to combine fasta files
        combine_fasta() {

        #s=$1
        #methods=$2

    # loop through methods
    for m in "${methods[@]}" ;
    do

      echo -e "\n${s}, ${m}\n"

      # rename header lines to reflect sample_ID, method, and contig number; append to combined fasta file


      # store file paths

      fasta_file_path="${PWD}/${OUTPREFIX}tmp/${s}.${m}.input.fasta"

      contig_summary_path="${PWD}/${OUTPREFIX}tmp/${s}.${m}.contig_summary.csv"


      # if fasta file does not exist, skip to next method
      if [ ! -s "$fasta_file_path" ]; then
        echo -e "\nWARNING: fasta file does not exist or is empty for ${s}, ${m}; skipping to next\n"

	printf "${s},${m}" >> "${PWD}/${OUTPREFIX}tmp_logs/no_fasta_file.txt"

        continue
      fi

      # if contig summary exists, order contig names by length

	if [ -s "$contig_summary_path" ]; then

      		# order contig names by length

      		# replace spaces with underscores in contig summary
      		sed 's/ /_/g' ${contig_summary_path} > ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.contig_summary.fasta

      		# remove quotes in contig summary
      		sed -i 's/\"//g' ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.contig_summary.fasta

      		contig_summary_path="${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.contig_summary.fasta"

      		# read contig summary as array without header row
      		readarray -t contig_names_lengths <<< `tail -n +2 ${contig_summary_path}`

      		# sort array by base pair length
      		contig_names_lengths_sorted=( `printf '%s\n' "${contig_names_lengths[@]}" | sort -t , -k 2 -g` )

      		# create or reset name object
      		contig_names=( )
      		contig_circularities=( )

      		# split names from rest of array
      		for i in "${!contig_names_lengths_sorted[@]}"; do

			contig_name=`echo ${contig_names_lengths_sorted[i]} | cut -d ',' -f 1`
        		contig_circularity=`echo ${contig_names_lengths_sorted[i]} | cut -d ',' -f 3`

		 	contig_names+=( "${contig_name}" )
        		contig_circularities+=( "${contig_circularity}" )

      		done

	# else if contig summary does not exist, take contigs in fasta file order and assume all linear
	else

		for i in `grep ">" ${fasta_file_path}` ; do

			# check for header marker as space before topology will lead to reading topology as separate array element
			if `echo $i | grep -q ">"` ; then
				contig_names+=( ${i} )
				contig_circularities+=( "no" )
			fi

		done

	fi

      # replace spaces with semicolons in fasta file
      sed 's/ /;/g' ${fasta_file_path} > ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.input.fasta

      fasta_file_path="${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.input.fasta"

      # get headers with line numbers and store as string
      header_lines=$(printf -- '%s\n' | grep -n -i ">" $fasta_file_path | tr '\n' ' ')

      # convert string to array
      header_lines=( $header_lines )


      # loop over ordered contig header lines

      for n in "${!contig_names[@]}" ;
      do

        # find header line indice that matches ordered contig name

        for l in "${!header_lines[@]}" ;
        do
          if echo ${header_lines[l]} | grep -q ">${contig_names[n]};" ; then
            header_line_index=${l}
            break
          fi
        done

        # get header line number
        header_line_number=`echo ${header_lines[${header_line_index}]} | cut -d : -f 1`

        # replace header line
        sed "${header_line_number} s/>.*/>${s}.${m}.$n/" $fasta_file_path > ${PWD}/${OUTPREFIX}tmp/tmp.2.${s}.${m}.input.fasta
        rm "${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.input.fasta"
        mv "${PWD}/${OUTPREFIX}tmp/tmp.2.${s}.${m}.input.fasta" "${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.input.fasta"

        # if contig summary exists, get read_depth header line -> mt1 is shortest contig, mt2 next shortest etc.
        
	if [ -s "$contig_summary_path" ]; then
		read_depth_label=`echo "mt$(($n+1))"`

       		# track header line changes -> new_header,original_line_number,original_contig_name,read_depth_name
	        printf "${s}.${m}.$n,${header_line_number},${contig_names[$n]},${read_depth_label}\n" >> "${PWD}/${OUTPREFIX}tmp/header_changes.csv"
	else
		# report NA  for read_depth_label
		printf "${s}.${m}.$n,${header_line_number},${contig_names[$n]},NA\n" >> "${PWD}/${OUTPREFIX}tmp/header_changes.csv"
	fi

      done

############## ADD REORDERING OF FILE HERE ##################

############## CHECK IF CONTIG IS CIRCULAR -> one entry in contig_names_lengths, plus yes value in circularity
# make subdirectory for reorder output
mkdir -p ${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}

if [ "${#contig_circularities[@]}" = 1 ] && [ "${contig_circularities[0]}" = "yes" ] ; then

	# reorder with circularity on -> output to subdirectory

	# if read_depth file exists, include header_changes.csv and read_depth.txt files
	if [ -s "${PWD}/${OUTPREFIX}tmp/${s}.${m}.read_depth.txt" ]; then
        	bash $PWD/auto_scripts/reorder_nucleotides.sh -c "yes" -f ${fasta_file_path} -r ${MT_REFERENCE} -i 6 -d ${PWD}/${OUTPREFIX}tmp/${s}.${m}.read_depth.txt -h ${PWD}/${OUTPREFIX}tmp/header_changes.csv -o ${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}
	else
		bash $PWD/auto_scripts/reorder_nucleotides.sh -c "yes" -f ${fasta_file_path} -r ${MT_REFERENCE} -i 6 -o ${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}
	fi
else
        # reorder with circularity off -> output to subdirectory

	# if read_depth file exists, include header_changes.csv and read_depth.txt files
	if [ -s "${PWD}/${OUTPREFIX}tmp/${s}.${m}.read_depth.txt" ]; then
	
        	bash $PWD/auto_scripts/reorder_nucleotides.sh -c "no" -f ${fasta_file_path} -r ${MT_REFERENCE} -i 6 -d ${PWD}/${OUTPREFIX}tmp/${s}.${m}.read_depth.txt -h ${PWD}/${OUTPREFIX}tmp/header_changes.csv -o ${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}
	else
		bash $PWD/auto_scripts/reorder_nucleotides.sh -c "no" -f ${fasta_file_path} -r ${MT_REFERENCE} -i 6 -o ${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}
	fi

fi

# check if realigned.fasta exists in subdirectory
if [ -f ${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}/realigned.fasta ] ; then

        #overwrite fasta file path with realigned.fasta
        #fasta_file_path="${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}/realigned.fasta"

        # rename realigned fasta with sample and method
        mv ${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}/realigned.fasta ${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}/${s}.${m}.realigned.fasta

        # copy realigned fasta to main tmp folder
        cp ${PWD}/${OUTPREFIX}tmp/reorder/${s}/${m}/${s}.${m}.realigned.fasta ${PWD}/${OUTPREFIX}tmp/

        # replace tmp.input fasta with realigned file
        mv ${PWD}/${OUTPREFIX}tmp/${s}.${m}.realigned.fasta ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.input.fasta

        #write sample and method to reordered_sequence.csv
        printf "${s},${m}\n" >> ${PWD}/${OUTPREFIX}tmp/reordered_sequence_list.csv

fi

############## END EDITING ##################################

      # add fasta file with edited headers to master fasta file

     # cat $fasta_file_path >> "${PWD}/${OUTPREFIX}tmp/combined_fasta.fasta"


    # close methods loop
    done
# close function
}

# call function to combine sequences into fasta file
for s in "${sample_IDs[@]}" ;
do

	if [[ $(jobs -r -p | wc -l) -ge $num_processes ]]; then wait -n ; fi
        #combine_fasta ${s} ${methods[@]} > ${PWD}/${OUTPREFIX}tmp_logs/${s}.combine_fasta_log 2>&1 &

	# skip samples with no input fasta
	#if ! [ -s ${PWD}/${OUTPREFIX}tmp/${s}.*.input.fasta ] ; then continue ; fi

        combine_fasta > ${PWD}/${OUTPREFIX}tmp_logs/${s}.combine_fasta_log 2>&1 &

# close sample loop
done
# wait for all child processes to finish
wait

# write edited fasta files to combined fasta file
cat ${PWD}/${OUTPREFIX}tmp/tmp.*.*.input.fasta >> "${PWD}/${OUTPREFIX}tmp/combined_fasta.fasta"

# write logs to standard out
cat ${PWD}/${OUTPREFIX}tmp_logs/*.combine_fasta_log

  echo -e "\ncombining concluded\nrunning mafft..."

  # run mafft alignment on combined fasta file

  combined_fasta_path="${PWD}/${OUTPREFIX}tmp/combined_fasta.fasta"

  #$PWD/external_software/mafft-linux64/mafft.bat --adjustdirection --6merpair --thread -1 --addfragments ${combined_fasta_path} ${MT_REFERENCE} > "${PWD}/${OUTPREFIX}tmp/combined_sequence_mafft.fasta"

# limiting threads to 8 for server run
$PWD/external_software/mafft-linux64/mafft.bat --adjustdirectionaccurately --6merpair --thread 8 --addfragments ${combined_fasta_path} ${MT_REFERENCE} > "${PWD}/${OUTPREFIX}tmp/combined_sequence_mafft.fasta"

  echo -e "\nmafft complete\nwriting read depths file and aligned nucleotides csv..."

  combined_mafft_path="${PWD}/${OUTPREFIX}tmp/combined_sequence_mafft.fasta"

  # create tmp mafft without reference header
  sed '1d' ${combined_mafft_path} > "${PWD}/${OUTPREFIX}tmp/tmp.combined_sequence_mafft.fasta"

  combined_mafft_path="${PWD}/${OUTPREFIX}tmp/tmp.combined_sequence_mafft.fasta"

  # loop through mafft headers

  # get headers with line numbers and store as string
  header_lines=$(printf -- '%s\n' | grep -n -i ">" $combined_mafft_path | tr '\n' ' ')

  # convert string to array
  header_lines=($header_lines)

# create function to map read depths to aligned nucleotides

        map_read_depths() {

        #l=$1
        #combined_mafft_path=$2

    # isolate contig sequence

    # store contig name -> remove line number and >
    header_text_with_indicator=`echo ${header_lines[l]} | cut -d : -f 2`
    header_text=`echo ${header_lines[l]} | cut -d : -f 2 | cut -d ">" -f 2`

        echo -e "\nisolating ${l} sequence..."

    # store contig header line number
    header_line_number=`echo ${header_lines[l]} | cut -d : -f 1`

    # get tail lines from header
    tail -n "+$header_line_number" ${combined_mafft_path} > "${PWD}/${OUTPREFIX}tmp/tmp.${l}.sequence_mafft_prepSeq.fasta"

    # find line number of next header
    scndheader=`grep -n -i -m2 ">" "${PWD}/${OUTPREFIX}tmp/tmp.${l}.sequence_mafft_prepSeq.fasta" | cut -d : -f 1 | tail -n1`

    # check if line number of second header equals first line
    # if no, then take head lines of file up to before second header
    if [ ! $scndheader = 1 ]; then
      head -n $((${scndheader}-1)) "${PWD}/${OUTPREFIX}tmp/tmp.${l}.sequence_mafft_prepSeq.fasta" > "${PWD}/${OUTPREFIX}tmp/tmp.${l}.isolated_sequence_mafft_prepSeq.fasta"
    fi

    # if yes, then send file 2 to 3 for naming consistency
    if [ $scndheader = 1 ]; then
      mv "${PWD}/${OUTPREFIX}tmp/tmp.${l}.sequence_mafft_prepSeq.fasta" "${PWD}/${OUTPREFIX}tmp/tmp.${l}.isolated_sequence_mafft_prepSeq.fasta"
    fi

    # remove header line
    sed '1d' "${PWD}/${OUTPREFIX}tmp/tmp.${l}.isolated_sequence_mafft_prepSeq.fasta" > "${PWD}/${OUTPREFIX}tmp/tmp.${l}.no_head_isolated_sequence_mafft_prepSeq.fasta"

    # split all characters onto new lines
    sed 's/*\?/\n/g' "${PWD}/${OUTPREFIX}tmp/tmp.${l}.no_head_isolated_sequence_mafft_prepSeq.fasta" > "${PWD}/${OUTPREFIX}tmp/tmp.${l}.split_nucleotides_isolated_sequence_prepSeq.fasta"

    # remove blank lines
    sed '/^[[:space:]]*$/d' "${PWD}/${OUTPREFIX}tmp/tmp.${l}.split_nucleotides_isolated_sequence_prepSeq.fasta" > "${PWD}/${OUTPREFIX}tmp/tmp.${l}.split_nucleotides_isolated_sequence.fasta"

    # remove intermediate files -> ending in _prepSeq
    rm ${PWD}/${OUTPREFIX}tmp/tmp.${l}.*_prepSeq.fasta

    # get matching read depth file -> from header

        # use reverse to cut method and contig number from end of header
        method=`echo ${header_text} | rev | cut -d . -f 2 | rev`
        contig_number=`echo ${header_text} | rev | cut -d . -f 1 | rev`

        method_contig_number="${method}.${contig_number}"

        # trim off method and contig number from end of header to get sample_ID -> prevents splitting of samples containing periods
        sample_ID="${header_text%.${method_contig_number}}"

    # if sample_ID contains _R_, remove
    if echo ">$sample_ID" | grep -q ">_R_" ; then sample_ID="${sample_ID#_R_*}" ; fi

        echo -e "retrieving ${l} depth information..."

    depths_file_path="${PWD}/${OUTPREFIX}tmp/${sample_ID}.${method}.read_depth.txt"

    # report if no depth file exists
    if [ ! -s "${depths_file_path}" ]; then
	    echo "${depths_file_path} does not exist or is empty"
	    return 1
    fi

    # GET REORDERED DEPTHS FILE IF EXISTING
    # check if sample and method recorded in reordered_sequences_list.csv
    if grep -iq "${sample_ID},${method}" ${OUTPREFIX}tmp/reordered_sequence_list.csv ; then

        # if yes, overwrite depths path with reordered read depth file
        depths_file_path=${PWD}/${OUTPREFIX}tmp/reorder/${sample_ID}/${method}/reordered_read_depths.csv

    fi

    # get contig label in read depths file
    read_depth_label=`echo "mt$((${contig_number}+1))"`

    # create tmp read depths file
    grep "${read_depth_label}" < ${depths_file_path} > ${OUTPREFIX}tmp/tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt

    depths_file_path="${PWD}/${OUTPREFIX}tmp/tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt"

    # replace spaces in read depths file with commas
    sed 's/ /,/g' ${depths_file_path} > ${PWD}/${OUTPREFIX}tmp/tmp.2.${sample_ID}.${method}.${contig_number}.read_depth.txt
    rm ${PWD}/${OUTPREFIX}tmp/tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt
    mv ${PWD}/${OUTPREFIX}tmp/tmp.2.${sample_ID}.${method}.${contig_number}.read_depth.txt ${PWD}/${OUTPREFIX}tmp/tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt

    # check if header in mafft file contains >_R_
    if echo ${header_text_with_indicator} | grep -q ">_R_" ; then

      # invert read depth file as mafft sequence adjusted for direction
      tac ${depths_file_path} > ${OUTPREFIX}tmp/tmp.2.${sample_ID}.${method}.${contig_number}.read_depth.txt
      rm ${PWD}/${OUTPREFIX}tmp/tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt
      mv ${PWD}/${OUTPREFIX}tmp/tmp.2.${sample_ID}.${method}.${contig_number}.read_depth.txt ${PWD}/${OUTPREFIX}tmp/tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt

    fi

    # read depths file for contig as array
    readarray -t read_depths < ${depths_file_path}

    # loop through nucleotides in split nucleotides file

    readarray -t nucleotides < ${PWD}/${OUTPREFIX}tmp/tmp.${l}.split_nucleotides_isolated_sequence.fasta

    # create bp value
    bp_value=0

    echo -e "updating ${l} depth information to reflect mafft bp positions..."

    for n in "${!nucleotides[@]}" ;
    do

      # check if nucleotide is a hyphen
      if [ "${nucleotides[n]}" = "-" ]; then
        # assign read_depth_value as 0
        read_depth_value=0
      fi

      # check if nucleotide is not a hyphen
      if [ "${nucleotides[n]}" != "-" ]; then

        # access read depth array index equal to bp_value -> value one less than actual bp to access line index
        read_depth_bp_entry=${read_depths[${bp_value}]}

        # add one to base pair value
        bp_value=$((${bp_value}+1))

        # get read_depth
        read_depth_value=`echo ${read_depth_bp_entry} | cut -d , -f 4`

        # get original bp -> if sequence is inverted then not equal to bp_value
        original_bp=`echo ${read_depth_bp_entry} | cut -d , -f 2`

        # add mafft bp vs original bp info to .csv file -> new header,mafft bp,original bp
        printf "${header_text},$((${n}+1)),${original_bp}\n" >> "${PWD}/${OUTPREFIX}tmp/position_changes.csv"

      fi

      # append header text, bp info, and read depth value combined read depth csv for header
      printf "${header_text},$((${n}+1)),$((${n}+1)),${read_depth_value}\n" >> "${PWD}/${OUTPREFIX}tmp/keep_tmp.${l}.combined_read_depths.csv"

    done

        # append combined read depth csv for header to full combined_read_depths.csv
        # cat "${PWD}/${OUTPREFIX}tmp/tmp.${l}.combined_read_depths.csv" >> "${PWD}/${OUTPREFIX}tmp/combined_read_depths.csv"

    #convert split nucleotides file to single string delineated by commas
    nucleotide_sequence=$(printf -- '%s\n' | cat "${PWD}/${OUTPREFIX}tmp/tmp.${l}.split_nucleotides_isolated_sequence.fasta" | tr '\n' ',')
    #remove trailing comma
    nucleotide_sequence=${nucleotide_sequence%?}

    #concatenate nucleotide string to header name and write to tmp file
    printf "${sample_ID},${method},${contig_number},${nucleotide_sequence}\n" >> "${PWD}/${OUTPREFIX}tmp/keep_tmp.${l}.aligned_nucleotides.csv"

        # remove header specific tmp files -> added period to prevent tmp.1. matching to tmp.11. etc
        rm ${PWD}/${OUTPREFIX}tmp/tmp.${l}.*

        echo -e "complete!\n"

# close function
}

# call read depth mapping function
  for l in "${!header_lines[@]}" ;
  do

	if [[ $(jobs -r -p | wc -l) -ge $num_processes ]]; then wait -n ; fi
        #map_read_depths        ${l} ${combined_mafft_path} > ${PWD}/${OUTPREFIX}tmp_logs/${l}.map_read_depths_log 2>&1 &
        map_read_depths > ${PWD}/${OUTPREFIX}tmp_logs/${l}.map_read_depths_log 2>&1 &

done
# wait for all child processes to finish
wait

# write logs to standard out
${PWD}/${OUTPREFIX}tmp_logs/*.map_read_depths_log

  echo -e "\nread depths file and aligned nucleotides csv complete"

  # append files to header lines where needed and copy to stats folder

  printf "contig,aligned_bp_position,original_bp_position\n" > "${PWD}/${OUTPREFIX}tmp/tmp.position_changes_header.csv"

  printf "aligned_header,original_header_line_position,original_contig_name,original_read_depth_label\n" > "${PWD}/${OUTPREFIX}tmp/tmp.header_changes_header.csv"

  printf "contig,bp_1,bp_2,read_depth\n" > "${PWD}/${OUTPREFIX}tmp/tmp.combined_read_depths_header.csv"

  cat "${PWD}/${OUTPREFIX}tmp/tmp.position_changes_header.csv" "${PWD}/${OUTPREFIX}tmp/position_changes.csv" > "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}position_changes.csv"
  rm "${PWD}/${OUTPREFIX}tmp/position_changes.csv" "${PWD}/${OUTPREFIX}tmp/tmp.position_changes_header.csv"

  cat "${PWD}/${OUTPREFIX}tmp/tmp.header_changes_header.csv" "${PWD}/${OUTPREFIX}tmp/header_changes.csv" > "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}header_changes.csv"
  rm "${PWD}/${OUTPREFIX}tmp/header_changes.csv" "${PWD}/${OUTPREFIX}tmp/tmp.header_changes_header.csv"

  cat "${PWD}/${OUTPREFIX}tmp/tmp.combined_read_depths_header.csv" ${PWD}/${OUTPREFIX}tmp/keep_tmp.*.combined_read_depths.csv > "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}combined_read_depths.csv"
  rm "${PWD}/${OUTPREFIX}tmp/tmp.combined_read_depths_header.csv" ${PWD}/${OUTPREFIX}tmp/keep_tmp.*.combined_read_depths.csv

  # create header for aligned_nucleotides.csv

  max_bp=`tail -n1 "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}combined_read_depths.csv" | cut -d , -f 2`

  printf "sample,method,contig_number" > ${PWD}/${OUTPREFIX}tmp/aligned_nucleotides_header.csv

  current_bp=0

  while [[ "${current_bp}" -lt "${max_bp}" ]] ; do current_bp=$(( current_bp+1 )) ; printf ",bp_${current_bp}" >> ${PWD}/${OUTPREFIX}tmp/aligned_nucleotides_header.csv ; done

  printf "\n" >> ${PWD}/${OUTPREFIX}tmp/aligned_nucleotides_header.csv

  cat ${PWD}/${OUTPREFIX}tmp/aligned_nucleotides_header.csv ${PWD}/${OUTPREFIX}tmp/keep_tmp.*.aligned_nucleotides.csv > ${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}aligned_nucleotides.csv
  rm ${PWD}/${OUTPREFIX}tmp/aligned_nucleotides_header.csv ${PWD}/${OUTPREFIX}tmp/keep_tmp.*.aligned_nucleotides.csv

mv "${PWD}/${OUTPREFIX}tmp/combined_sequence_mafft.fasta" "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}combined_sequence_mafft.fasta"
mv "${PWD}/${OUTPREFIX}tmp/combined_fasta.fasta" "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}combined_fasta.fasta"

  echo -e "\ncopying output to stats folder..."

  cp "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}position_changes.csv" "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}header_changes.csv" "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}combined_read_depths.csv" "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}aligned_nucleotides.csv" "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}combined_sequence_mafft.fasta" "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}combined_fasta.fasta" ${PWD}/benchmark_stats/

  echo -e "\ncopying complete"

# remove remaining .tmp files
  rm ${PWD}/${OUTPREFIX}tmp/tmp.*

        # collate individual genes, run alignments, map read depths

        # read expected genes into an array, sans header row
        expected_genes=$(tail -n +2 ${PWD}/${OUTPREFIX}tmp/expected_genes.csv | tr '\n' ' ')
        expected_genes=( ${expected_genes} )

        # write function to separate gene sequences from annotated fasta files
        separate_fasta_gene() {
                #s=$1
                #method=$2
                #expected_genes=$3

                # loop through methods
                for m in "${methods[@]}" ;
                do

                        echo -e "\n${s}, ${m}\n"

                        # get annotated fasta file path
                        annotated_fasta_path="${PWD}/${OUTPREFIX}tmp/${s}.${m}.annotated.fasta"

                        # check path to annotated fasta exists, else continue
                        if [ ! -s ${PWD}/${OUTPREFIX}tmp/${s}.${m}.annotated.fasta ]; then
                                echo -e "\nWARNING: annotated fasta file not found for ${s}, ${m}\n"
                                continue
                        fi

                        echo -e "getting annotated fasta from: ${annotated_fasta_path}"

                        # create tmp file
                        cp ${annotated_fasta_path} "${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.annotated.fasta"
                        annotated_fasta_path="${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.annotated.fasta"

                        # rename contig name in header lines to new name -> reflects sample, method, contig number

                        # store header lines and line numbers
                        header_lines=$(printf -- '%s\n' | grep -n -i ">" ${annotated_fasta_path} | tr '\n' ' ')
                        # store as array
                        header_lines=( ${header_lines} )

                        echo -e "\nupdating contig names...\n"

                        # loop through header lines
                        for l in "${header_lines[@]}";
                        do

                                # get header line number
                                header_line_number=`echo ${l} | cut -d : -f 1`

                                # get original contig name
                                original_contig_name=`echo ${l} | cut -d ";" -f 1 | cut -d ">" -f 2`

                                # find line with original header name, sample and method in header change.csv
                                new_header_match=`grep -i ",${original_contig_name}," ${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}header_changes.csv | grep -i ".${m}." | grep -i "^${s}."`

                                # get new header
                                new_header=`echo ${new_header_match} | cut -d , -f 1`

                                # cleave original header name from header line
                                header_remainder="${l##${header_line_number}:>${original_contig_name};}"

                                # paste new header and header remainder together
                                new_header_full=`printf '%s' "${new_header};${header_remainder}"`

                                #echo -e "\nupdating ${original_contig_name} to ${new_header_full}"

                                # replace header line
                                sed "${header_line_number} s/>.*/>${new_header_full}/" ${annotated_fasta_path} > ${PWD}/${OUTPREFIX}tmp/tmp.2.${s}.${m}.annotated.fasta
                                rm "${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.annotated.fasta"
                                mv "${PWD}/${OUTPREFIX}tmp/tmp.2.${s}.${m}.annotated.fasta" "${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.annotated.fasta"
                        done

                        # loop through genes
                        for g in "${expected_genes[@]}" ;
                        do

                                # separate gene from expected frequency
                                gene=`echo "${g}" | cut -d , -f 1`

                                echo -e "\nlocating header lines for: ${gene}...\n"

                                # find header lines containing gene
                                gene_header_lines=$(printf -- '%s\n' | grep -n -i ">" ${annotated_fasta_path} | grep "${gene}" | tr '\n' ' ')

                                # convert string to array
                                gene_header_lines=( ${gene_header_lines} )

                                echo -e "\nconfirming target gene, ${gene}, does not match to any other expected gene...\n"

                                # check if gene matches to any other expected gene
                                gene_matches=$(printf -- '%s\n' | grep -i "${gene}" ${OUTPREFIX}tmp/expected_genes.csv | tr '\n' ' ')
                                # convert string to array
                                gene_matches=( ${gene_matches} )

                                # if there is more than one gene match
                                if [[ "${#gene_matches[@]}" -gt 1 ]]; then

                                        echo -e "\n${gene} matches to ${#gene_matches[@]} expected genes\nfiltering incorrect gene matches..."

                                        # create empty, incorrect match storage array
                                        incorrect_matches=( )

                                        # get incorrect gene matches
                                        for i in "${!gene_matches[@]}";
                                        do

                                                # separate gene from expected frequency
                                                query_gene=`echo ${gene_matches[${i}]} | cut -d , -f 1`

                                                # check if gene matches exactly to target gene, if yes continue
                                                if echo "${query_gene}" | grep -q "^${gene}$" ; then
                                                        continue
                                                else
                                                        # write off target gene to incorrect match array
                                                        incorrect_matches+=( ${query_gene} )
                                                fi

                                        done

                                        # separate header lines matching incorrect genes from those matching target gene

                                        # create tmp storage matching header
                                        tmp_gene_headers=( ${gene_header_lines[@]} )

                                        # loop through incorrect gene matches
                                        for i in "${incorrect_matches[@]}";
                                        do

                                                # create second tmp storage
                                                filtered_matches=( )

                                                # loop through header lines in tmp storage
                                                for h in "${tmp_gene_headers[@]}";
                                                do

                                                        # check if header line matches to incorrect gene
                                                        if echo ${h} | grep -q ${i} ; then
                                                                # if yes continue
                                                                continue
                                                        else
                                                                # if not add to filtered matches
                                                                filtered_matches+=( ${h} )
                                                        fi

                                                done

                                                # overwrite tmp match storage with filtered matches to enable filtering of next incorrect gene
                                                tmp_gene_headers=( ${filtered_matches[@]} )

                                        done

                                        # overwite header lines with filtered headers
                                        gene_header_lines=( ${tmp_gene_headers[@]} )

                                        echo -e "\nfiltering complete!\n"

                                # close if multiple gene matches
                                fi

                                echo -e "\nseparating ${gene} sequences from annotated fasta file...\n"

				# if NO header lines found for gene, skip to next gene -> likely unnecessary
				#if [ "${gene_header_lines[0]}" == "" ]; then
				#	echo "WARNING: No $g annotations found for $s in $m"
				#	continue
				#fi

                                # loop through header lines
                                for l in "${gene_header_lines[@]}";
                                do

                                        # get header line number
                                        header_line_number=$(echo ${l} | cut -d : -f 1)

                                        # get tail of annotated fasta file from header line
                                        tail -n +${header_line_number} ${annotated_fasta_path} > ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.${gene}.fasta

                                        # find second header line number
                                        scnd_header_line=$(grep -ni -m2 ">" ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.${gene}.fasta | tail -n1 | cut -d : -f 1)

                                        # if second header line number not equal to 1, remove lines from second header
                                        if [ ! ${scnd_header_line} = 1 ] ; then

                                                # get head of file up to second header
                                                head -n $((${scnd_header_line}-1)) ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.${gene}.fasta > ${PWD}/${OUTPREFIX}tmp/tmp.2.${s}.${m}.${gene}.fasta
                                                rm ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.${gene}.fasta
                                                mv ${PWD}/${OUTPREFIX}tmp/tmp.2.${s}.${m}.${gene}.fasta ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.${gene}.fasta

                                        fi

                                        # append isolated gene to combined gene specific file for sample and method
                                        cat ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.${gene}.fasta >> ${PWD}/${OUTPREFIX}tmp/keep_tmp.${s}.${m}.${gene}.fasta
                                        rm ${PWD}/${OUTPREFIX}tmp/tmp.${s}.${m}.${gene}.fasta

                                done

                                echo -e "\nisolation complete!\n"

                        # close expected gene loop
                        done

                # close method loop
                done

        # close function
        }

# call separate gene fasta function
for s in "${sample_IDs[@]}" ;
do

	if [[ $(jobs -r -p | wc -l) -ge $num_processes ]]; then wait -n ; fi
        #separate_fasta_gene ${s} ${methods[@]} ${expected_genes[@]} > ${PWD}/${OUTPREFIX}tmp_logs/${s}.separate_fasta_gene_log 2>&1 &
        separate_fasta_gene > ${PWD}/${OUTPREFIX}tmp_logs/${s}.separate_fasta_gene_log 2>&1 &

done
# wait for all child processes to finish
wait

# write logs to standard out
cat ${PWD}/${OUTPREFIX}tmp_logs/*.separate_fasta_gene_log

echo -e "\ncompiling isolated gene sequences across assemblies...\n"

# loop through expected genes, concatenate tmp gene files together into gene specific file
for g in "${expected_genes[@]}" ;
do
        # separate gene from expected frequency
        gene=$(echo "${g}" | cut -d , -f 1)

        #cat sample gene fasta files together into one gene fasta file
        cat ${PWD}/${OUTPREFIX}tmp/keep_tmp.*.*.${gene}.fasta >> ${PWD}/${OUTPREFIX}tmp/combined_${gene}_fasta.fasta
        rm ${PWD}/${OUTPREFIX}tmp/keep_tmp.*.*.${gene}.fasta

done

echo -e "\ncompilation complete!\n"

# if trna fastas exist, combine then
if ls ${PWD}/${OUTPREFIX}tmp/combined_trn*_fasta.fasta >/dev/null 2>&1; then
echo -e "\nadditionally compiling tRNAs into single fasta...\n"

cat ${PWD}/${OUTPREFIX}tmp/combined_trn*_fasta.fasta >> ${PWD}/${OUTPREFIX}tmp/combined_tRNA_fasta.fasta
fi

echo -e "\ncompilation complete!\n"

# run alignments on gene fasta files

for g in "${expected_genes[@]}";
do

        # separate gene name from expected frequency
        gene=`echo ${g} | cut -d , -f 1`

        # isolate gene sequence from annotated reference genome

        # get prefix from annotated reference file -> remove file path and type
        annotated_prefix=$(basename ${MT_REFERENCE_ANNOTATED%.*})

        # make reference directory if needed
        mkdir -p ${PWD}/reference

        # check that gene reference does not exist or is empty in reference directory

        if [ ! -s ${PWD}/reference/${annotated_prefix}_${gene}.fasta ]; then

                echo -e "\nseparating ${gene} annotated reference from ${MT_REFERENCE_ANNOTATED}\n"

                # find header lines containing gene
                gene_header_lines=$(printf -- '%s\n' | grep -n -i ">" ${MT_REFERENCE_ANNOTATED} | grep -i "${gene}" | tr '\n' ' ')

                # convert string to array
                gene_header_lines=( ${gene_header_lines} )

                echo -e "\nconfirming target gene, ${gene}, does not match to any other expected gene...\n"

                # check if gene matches to any other expected gene
                gene_matches=$(printf -- '%s\n' | grep -i "${gene}" ${OUTPREFIX}tmp/expected_genes.csv | tr '\n' ' ')

                # convert string to array
                gene_matches=( ${gene_matches} )

                # if there is more than one gene match
                if [[ "${#gene_matches[@]}" -gt 1 ]]; then

                        echo -e "\n${gene} matches to ${#gene_matches[@]} expected genes\nfiltering incorrect gene matches..."

                        # create empty, incorrect match storage array
                        incorrect_matches=( )

                        # get incorrect gene matches
                        for i in "${!gene_matches[@]}";
                        do
                                # separate gene from expected frequency
                                query_gene=`echo ${gene_matches[${i}]} | cut -d , -f 1`

                                # check if gene matches exactly to target gene, if yes continue
                                if echo "${query_gene}" | grep -q "^${gene}$" ; then
                                        continue
                                else
                                        # write off target gene to incorrect match array
                                        incorrect_matches+=( ${query_gene} )
                                fi
                        done

                        # separate header lines matching incorrect genes from those matching target gene

                        # create tmp storage matching header
                        tmp_gene_headers=( ${gene_header_lines[@]} )

                        # loop through incorrect gene matches
                        for i in "${incorrect_matches[@]}";
                        do

                                # create second tmp storage
                                filtered_matches=( )

                                # loop through header lines in tmp storage
                                for h in "${tmp_gene_headers[@]}";
                                do

                                        # check if header line matches to incorrect gene
                                        if echo ${h} | grep -q ${i} ; then
                                                # if yes continue
                                                continue
                                        else
                                                # if not add to filtered matches
                                                filtered_matches+=( ${h} )
                                        fi
                                done

                                # overwrite tmp match storage with filtered matches to enable filtering of next incorrect gene
                                tmp_gene_headers=( ${filtered_matches[@]} )

                        done

                        # overwite header lines with filtered headers
                        gene_header_lines=( ${tmp_gene_headers[@]} )

                        echo -e "\nfiltering complete!\n"

                # close if multiple gene matches
                fi

                # get only first header for annotated reference -> duplicate genes represented by one gene header
                l=${gene_header_lines[0]}

                echo -e "\nseparating ${gene} reference sequence from annotated reference fasta...\n"

                # get header line number
                header_line_number=$(echo ${l} | cut -d : -f 1)

                # get tail of annotated fasta file from header line
                tail -n +${header_line_number} ${MT_REFERENCE_ANNOTATED} > ${PWD}/${OUTPREFIX}tmp/tmp.reference.${gene}.fasta

                # find second header line number
                scnd_header_line=$(grep -ni -m2 ">" ${PWD}/${OUTPREFIX}tmp/tmp.reference.${gene}.fasta | tail -n1 | cut -d : -f 1)

                # if second header line number not equal to 1, remove lines from second header
                if [ ! ${scnd_header_line} = 1 ] ; then

                        # get head of file up to second header
                        head -n $((${scnd_header_line}-1)) ${PWD}/${OUTPREFIX}tmp/tmp.reference.${gene}.fasta > ${PWD}/${OUTPREFIX}tmp/tmp.reference.2.${gene}.fasta
                        mv ${PWD}/${OUTPREFIX}tmp/tmp.reference.2.${gene}.fasta ${PWD}/${OUTPREFIX}tmp/tmp.reference.${gene}.fasta

                fi

                # move isolated sequence to reference folder
                mv ${PWD}/${OUTPREFIX}tmp/tmp.reference.${gene}.fasta ${PWD}/reference/${annotated_prefix}_${gene}.fasta

                echo -e "\nisolation complete!\n"

                # remove tmp files
                rm ${PWD}/${OUTPREFIX}tmp/tmp.reference.*

        # if reference gene sequence already exists in reference folder
        else

                echo -e "\nreference ${gene} sequence already extracted from ${MT_REFERENCE_ANNOTATED}\n"

        fi

        echo -e "\nrunning mafft for ${gene}...\n"

        # run mafft for gene against reference with fixed direction
        $PWD/external_software/mafft-linux64/mafft.bat --6merpair --thread -1 --addfragments ${PWD}/${OUTPREFIX}tmp/combined_${gene}_fasta.fasta ${PWD}/reference/${annotated_prefix}_${gene}.fasta > ${PWD}/${OUTPREFIX}tmp/combined_${gene}_mafft.fasta

        echo -e "\nmafft complete!\n"

done

# run mafft on combined tRNA sequences -> all reference sequences

# check that trna inputs exist
if ls ${PWD}/${OUTPREFIX}tmp/combined_tRNA_fasta.fasta >/dev/null 2>&1; then

echo -e "\nadditionally running mafft on combined tRNA sequences, including reference sequences\n"

cat ${PWD}/reference/${annotated_prefix}_trn*.fasta ${PWD}/${OUTPREFIX}tmp/combined_tRNA_fasta.fasta >> ${PWD}/${OUTPREFIX}tmp/combined_tRNA_incRef_fasta.fasta

$PWD/external_software/mafft-linux64/mafft.bat --maxiterate 5 --thread 8 ${PWD}/${OUTPREFIX}tmp/combined_tRNA_incRef_fasta.fasta > ${PWD}/${OUTPREFIX}tmp/combined_tRNA_mafft.fasta

echo -e "\nmafft complete!\n"

fi

# create function to map read depths for individual genes to aligned nucleotides

map_gene_read_depths() {
        # header line index
        #l=$1
        #combined_mafft_path=$2
        #gene=$3

        # create tmp mafft -> replace : with - which enables removal of header line number without removal of bp information BUT NOTE THAT HEADER LINES ARE STORED BEFORE CHANGE IS MADE
        #sed 's/:/-/' ${combined_mafft_path} > ${PWD}/${OUTPREFIX}tmp/tmp.combined_${gene}_mafft.fasta
        #combined_mafft_path=${PWD}/${OUTPREFIX}tmp/tmp.combined_${gene}_mafft.fasta

        # isolate contig sequence

        echo -e "\nisolating ${l} sequence..."

        # store contig name -> remove line number and >
        header_text_with_indicator=`echo ${header_lines[l]} | cut -d : -f 2`
        header_text=`echo ${header_lines[l]} | cut -d : -f 2 | cut -d ">" -f 2`

        # store contig header line number
        header_line_number=`echo ${header_lines[l]} | cut -d : -f 1`

        # get tail lines from header
        #tail -n "+$header_line_number" ${combined_mafft_path} > "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.sequence_mafft_prepSeq.fasta"

        # find line number of next header
        #scndheader=`grep -n -i -m2 ">" "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.sequence_mafft_prepSeq.fasta" | cut -d : -f 1 | tail -n1`

        # check if line number of second header equals first line
        # if no, then take head lines of file up to before second header
        #if ! [ "$scndheader" == 1 ]; then
        #        head -n $((${scndheader}-1)) "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.sequence_mafft_prepSeq.fasta" > "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.isolated_sequence_mafft_prepSeq.fasta"
        #fi

        # if yes, then send file 2 to 3 for naming consistency
        #if [ "$scndheader" == 1 ]; then
        #        mv "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.sequence_mafft_prepSeq.fasta" "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.isolated_sequence_mafft_prepSeq.fasta"
        #fi

        # remove header line
        #sed '1d' "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.isolated_sequence_mafft_prepSeq.fasta" > "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.no_head_isolated_sequence_mafft_prepSeq.fasta"

	# isolate sequence (need number lines in file)
	# isolate sequence
	isolated_sequence=`tail -n "+${header_line_number}" ${combined_mafft_path} | grep ">" -m2 -B 2000000000 | grep -v ">"`

	# if isolated sequence is blank, just get tail of file as last header
	if [ "${isolated_sequence}" == "" ] ; then
		isolated_sequence=`tail -n "+${header_line_number}" ${ALIGNMENT_FILE} | grep -v ">" | tr -cd 'a-zA-Z\-'`
	fi

	# write to file
	echo ${isolated_sequence} > "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.no_head_isolated_sequence_mafft_prepSeq.fasta"

        # split all characters onto new lines
        sed 's/*\?/\n/g' "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.no_head_isolated_sequence_mafft_prepSeq.fasta" > "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.split_nucleotides_isolated_sequence_prepSeq.fasta"

        # remove blank lines
        sed '/^[[:space:]]*$/d' "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.split_nucleotides_isolated_sequence_prepSeq.fasta" > "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.split_nucleotides_isolated_sequence.fasta"

        # remove intermediate files -> ending in _prepSeq
        rm ${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.*_prepSeq.fasta

        # get matching read depth file -> from header

        # extract contig origin header information from gene information
        contig_origin=`echo ${header_text} | cut -d ";" -f 1`

        # use reverse to cut method and contig number from end of contig_origin header
        method=`echo ${contig_origin} | rev | cut -d . -f 2 | rev`
        contig_number=`echo ${contig_origin} | rev | cut -d . -f 1 | rev`
        method_contig_number="${method}.${contig_number}"

        # trim off method and contig number from end of contig_origin header to get sample_ID -> prevents splitting of samples containing periods
        sample_ID="${contig_origin%.${method_contig_number}}"

        # if sample_ID contains _R_, remove
        if echo ">$sample_ID" | grep -q ">_R_" ; then sample_ID="${sample_ID#_R_*}" ; fi

                echo -e "retrieving ${l} depth information..."

                # store original depths file path -> reordering of complete sequence not relevant
                depths_file_path="${PWD}/${OUTPREFIX}tmp/${sample_ID}.${method}.read_depth.txt"

                # get contig label in read depths file
                read_depth_label=`echo "mt$((${contig_number}+1))"`

                # if contig specific depths file not existing or empty, create contig specific depths file
                if [ ! -s ${OUTPREFIX}tmp/keep_tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt ]; then

                # create tmp read depths file
                grep "${read_depth_label}" < ${depths_file_path} > ${OUTPREFIX}tmp/keep_tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt
                depths_file_path="${PWD}/${OUTPREFIX}tmp/keep_tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt"

                # replace spaces in read depths file with commas
                sed 's/ /,/g' ${depths_file_path} > ${PWD}/${OUTPREFIX}tmp/tmp.2.${sample_ID}.${method}.${contig_number}.read_depth.txt
                rm ${PWD}/${OUTPREFIX}tmp/keep_tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt
                mv ${PWD}/${OUTPREFIX}tmp/tmp.2.${sample_ID}.${method}.${contig_number}.read_depth.txt ${PWD}/${OUTPREFIX}tmp/keep_tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt

        else
                depths_file_path="${PWD}/${OUTPREFIX}tmp/keep_tmp.${sample_ID}.${method}.${contig_number}.read_depth.txt"
        fi

        # get first bp and last bp info from header -> originally stored header has : and not - separating bp information
        first_bp=`echo ${header_lines[l]} | cut -d ";" -f 4 | cut -d "[" -f 2 | cut -d "]" -f 1 | cut -d ":" -f 1`
        last_bp=`echo ${header_lines[l]} | cut -d ";" -f 4 | cut -d "[" -f 2 | cut -d "]" -f 1 | cut -d ":" -f 2`

        # if bp has < or > information, remove it
        if echo ${first_bp} | grep -q "<" ; then
                first_bp=`echo ${first_bp} | cut -d "<" -f 2`
        fi
        if echo ${first_bp} | grep -q ">" ; then
                first_bp=`echo ${first_bp} | cut -d ">" -f 2`
        fi
        if echo ${last_bp} | grep -q ">" ; then
                last_bp=`echo ${last_bp} | cut -d ">" -f 2`
        fi
        if echo ${last_bp} | grep -q "<" ; then
                last_bp=`echo ${last_bp} | cut -d "<" -f 2`
        fi

        gene_length=$(( ${last_bp} - ${first_bp} ))

        # create tmp depths file limited to bp lines
        tail -n +${first_bp} ${depths_file_path} | head -n ${gene_length} > ${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.read_depth.txt
        depths_file_path=${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.read_depth.txt

        # get gene direction information
        gene_direction=`echo $header_text_with_indicator | cut -d "(" -f 2 | cut -d ")" -f 1`

        # check if gene direction is '-'
        if [ "${gene_direction}" = "-" ] ; then
                # invert read depth file as gene is on reverse strand
                tac ${depths_file_path} > ${PWD}/${OUTPREFIX}tmp/tmp.2.${gene}.${l}.read_depth.txt
                rm ${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.read_depth.txt
                mv ${PWD}/${OUTPREFIX}tmp/tmp.2.${gene}.${l}.read_depth.txt ${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.read_depth.txt
        fi

        # read depths file for contig as array
        readarray -t read_depths < ${depths_file_path}

        # loop through nucleotides in split nucleotides file
        readarray -t nucleotides < ${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.split_nucleotides_isolated_sequence.fasta

        # create bp value
        bp_value=0

        echo -e "updating ${l} depth information to reflect mafft bp positions..."

        for n in "${!nucleotides[@]}" ;
        do
                # check if nucleotide is a hyphen
                if [ "${nucleotides[n]}" = "-" ]; then
                        # assign read_depth_value as 0
                        read_depth_value=0
                fi

                # check if nucleotide is not a hyphen
                if [ "${nucleotides[n]}" != "-" ]; then
                        # access read depth array index equal to bp_value -> value one less than actual bp to access line index
                        read_depth_bp_entry=${read_depths[${bp_value}]}

                        # add one to base pair value
                        bp_value=$((${bp_value}+1))

                        # get read_depth
                        read_depth_value=`echo ${read_depth_bp_entry} | cut -d , -f 4`

                        # get original bp -> if sequence is inverted then not equal to bp_value
                        original_bp=`echo ${read_depth_bp_entry} | cut -d , -f 2`

                        # add mafft bp vs original bp info to .csv file -> new header,mafft bp,original bp
                        printf "${header_text},$((${n}+1)),${original_bp}\n" >> "${PWD}/${OUTPREFIX}tmp/${gene}_position_changes.csv"
                fi

                # append header text, bp info, and read depth value combined read depth csv for header
                printf "${header_text},$((${n}+1)),$((${n}+1)),${read_depth_value}\n" >> "${PWD}/${OUTPREFIX}tmp/keep_tmp.${gene}.${l}.combined_read_depths.csv"

        done

        # append combined read depth csv for header to full combined_read_depths.csv
        # cat "${PWD}/${OUTPREFIX}tmp/tmp.${l}.combined_read_depths.csv" >> "${PWD}/${OUTPREFIX}tmp/combined_read_depths.csv"

        #convert split nucleotides file to single string delineated by commas
        nucleotide_sequence=$(printf -- '%s\n' | cat "${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.split_nucleotides_isolated_sequence.fasta" | tr '\n' ',')

        #remove trailing comma
        nucleotide_sequence=${nucleotide_sequence%?}

        # concatenate nucleotide string to header name and write to tmp file
        printf "${sample_ID},${method},${contig_number},${nucleotide_sequence}\n" >> "${PWD}/${OUTPREFIX}tmp/keep_tmp.${gene}.${l}.aligned_nucleotides.csv"

        # remove header specific tmp files -> added period to prevent tmp.1. matching to tmp.11. etc
        rm ${PWD}/${OUTPREFIX}tmp/tmp.${gene}.${l}.*

        echo -e "complete!\n"

# close function
}

for g in "${expected_genes[@]}";
do

        #separate gene from expected frequency
        gene=`echo ${g} | cut -d "," -f 1`

        echo -e "\nwriting read depths file and aligned nucleotides csv for ${gene}...\n"

        combined_mafft_path="${PWD}/${OUTPREFIX}tmp/combined_${gene}_mafft.fasta"

        # create tmp mafft without reference header
        sed '1d' ${combined_mafft_path} > "${PWD}/${OUTPREFIX}tmp/tmp.noRefHead.combined_${gene}_mafft.fasta"

        combined_mafft_path="${PWD}/${OUTPREFIX}tmp/tmp.noRefHead.combined_${gene}_mafft.fasta"

        # loop through mafft headers

        # get headers with line numbers and store as string then convert to array
        header_lines=$(printf -- '%s\n' | grep -n -i ">" $combined_mafft_path | tr '\n' ' ')
        header_lines=($header_lines)

        for l in "${!header_lines[@]}" ;
        do
                if [[ $(jobs -r -p | wc -l) -ge $num_processes ]]; then wait -n ; fi

                # call gene read depth mapping function
                #map_gene_read_depths ${l} ${combined_mafft_path} ${gene} > ${PWD}/${OUTPREFIX}tmp_logs/${gene}.${l}.map_gene_read_depths_log 2>&1 &
                map_gene_read_depths > ${PWD}/${OUTPREFIX}tmp_logs/${gene}.${l}.map_gene_read_depths_log 2>&1 &
        done
        # wait for all child processes to finish
        wait

        # write logs to standard out
        cat ${PWD}/${OUTPREFIX}tmp_logs/${gene}.*.map_gene_read_depths_log

        # append files to header lines where needed and copy to stats folder
        printf "contig,aligned_bp_position,original_bp_position\n" > "${PWD}/${OUTPREFIX}tmp/tmp.position_changes_header.csv"
        printf "contig,bp_1,bp_2,read_depth\n" > "${PWD}/${OUTPREFIX}tmp/tmp.combined_read_depths_header.csv"

        cat "${PWD}/${OUTPREFIX}tmp/tmp.position_changes_header.csv" "${PWD}/${OUTPREFIX}tmp/${gene}_position_changes.csv" > "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}${gene}_position_changes.csv"
        rm "${PWD}/${OUTPREFIX}tmp/tmp.position_changes_header.csv" "${PWD}/${OUTPREFIX}tmp/${gene}_position_changes.csv"

        cat "${PWD}/${OUTPREFIX}tmp/tmp.combined_read_depths_header.csv" ${PWD}/${OUTPREFIX}tmp/keep_tmp.${gene}.*.combined_read_depths.csv > "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}${gene}_combined_read_depths.csv"
        rm "${PWD}/${OUTPREFIX}tmp/tmp.combined_read_depths_header.csv" ${PWD}/${OUTPREFIX}tmp/keep_tmp.${gene}.*.combined_read_depths.csv

        # create header for aligned_nucleotides.csv

        max_bp=`tail -n1 "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}${gene}_combined_read_depths.csv" | cut -d , -f 2`

        printf "sample,method,contig_number" > ${PWD}/${OUTPREFIX}tmp/aligned_nucleotides_header.csv

        current_bp=0

        while [[ "${current_bp}" -lt "${max_bp}" ]] ; do current_bp=$(( current_bp+1 )) ; printf ",bp_${current_bp}" >> ${PWD}/${OUTPREFIX}tmp/aligned_nucleotides_header.csv ; done

        cat ${PWD}/${OUTPREFIX}tmp/aligned_nucleotides_header.csv ${PWD}/${OUTPREFIX}tmp/keep_tmp.${gene}.*.aligned_nucleotides.csv > ${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}${gene}_aligned_nucleotides.csv
        rm ${PWD}/${OUTPREFIX}tmp/aligned_nucleotides_header.csv ${PWD}/${OUTPREFIX}tmp/keep_tmp.${gene}.*.aligned_nucleotides.csv

        mv "${PWD}/${OUTPREFIX}tmp/combined_${gene}_mafft.fasta" "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}combined_${gene}_mafft.fasta"
        mv "${PWD}/${OUTPREFIX}tmp/combined_${gene}_fasta.fasta" "${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}combined_${gene}_fasta.fasta"

        echo -e "\ncopying output to stats folder..."

        cp ${PWD}/${OUTPREFIX}tmp/${OUTPREFIX}*${gene}* ${PWD}/benchmark_stats/

        echo -e "\ncopying complete"

        if ls ${PWD}/${OUTPREFIX}tmp/ | grep *tmp.${gene}.* ; then
                # remove gene tmp files
                rm ${PWD}/${OUTPREFIX}tmp/*tmp.${gene}.*
        fi

done

# close assembly stats if
fi

echo "script concluded"

exit 0
