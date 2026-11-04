#!/usr/bin/bash

# script to reorder nucleotides according to reference

TEMP=`getopt -o f:r:o:c:i:d:h: --long fasta:,reference:,outdir:,circular:,max_iterate:,read_depths:,header_match: -- "$@"`
eval set -- "$TEMP"
# extract options and their arguments into variables.
while true ; do
    case "$1" in
	-f|--fasta)
            case "$2" in
                "") shift 2 ;;
                *) FASTA=$2 ; shift 2 ;;
            esac ;;
	-r|--reference)
            case "$2" in
                "") shift 2 ;;
                *) MT_REFERENCE=$2 ; shift 2 ;;
            esac ;;
	-o|--outdir)
            case "$2" in
                "") shift 2 ;;
                *) OUTDIR=$2 ; shift 2 ;;
            esac ;;
	-c|--circular)
            case "$2" in
                "") CIRCULAR="no"; shift 2 ;;
                *) CIRCULAR=$2 ; shift 2 ;;
            esac ;;
	-i|--max_iterate)
            case "$2" in
                "") MAX_I=5; shift 2 ;;
                *) MAX_I=$2 ; shift 2 ;;
            esac ;;
	-d|--read_depths)
            case "$2" in
                "") READ_DEPTHS=""; shift 2 ;;
                *) READ_DEPTHS=$2 ; shift 2 ;;
            esac ;;
	-h|--header_match)
            case "$2" in
                "") HEADER_MATCH=""; shift 2 ;;
                *) HEADER_MATCH=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

if [ ! -n "$OUTDIR" ]; then
	OUTDIR="$PWD"
fi
if [ ! -n "$MAX_I" ]; then
	MAX_I=5
fi

# check fasta file has been given
if [ ! -n "$FASTA" ]; then
	echo -e "\nERROR: cannot access fasta file\n"
	exit 1
fi

echo -e "\ngetting fasta file from $FASTA\n"

# check that reference has been given

if [ ! -n "$MT_REFERENCE" ]; then
	echo -e "\nERROR: cannot access mitochondrial reference fasta file"
	exit 1
fi

echo -e "getting reference fasta from $MT_REFERENCE\n"

# count expected number of nucleotides
number_characters=`tail -n +2 ${MT_REFERENCE} | wc -m`
number_newlines=`tail -n +2 ${MT_REFERENCE} | wc -l`
expected_nucleotides=$((${number_characters} - ${number_newlines}))

echo -e "\nsupplied mitochondrial reference has ${expected_nucleotides} nucleotides\n"

# set circular default
if [ ! -n "$CIRCULAR" ]; then
	CIRCULAR="no"
fi

# check if read depths file given
if [ -n "$READ_DEPTHS" ]; then
	# check if file exists
	if [ ! -f "$READ_DEPTHS" ]; then
		echo -e "\nERROR: cannot access read depth file"
		exit 1
	fi
	# check if header match file exists
	if [ ! -n "$HEADER_MATCH" ] || [ ! -f "$HEADER_MATCH" ]; then
		echo -e "\nERROR: cannot access header match file"
		exit 1
	fi

	echo -e "getting read depth file from $READ_DEPTHS\n"
	echo -e "getting header match file from $HEADER_MATCH\n"

fi

# run initial mafft

echo -e "\nrunning initial mafft alignment...\n"

# run mafft against single reference
$PWD/external_software/mafft-linux64/mafft.bat --adjustdirectionaccurately --6merpair --thread -1 --addfragments $FASTA $MT_REFERENCE > ${OUTDIR}/initial_mafft.fasta

# assign initial mafft to $MAFFT
MAFFT=${OUTDIR}/initial_mafft.fasta

# check for mafft failure
if [ ! -s ${OUTDIR}/initial_mafft.fasta ]; then
	echo -e "\nERROR: initial mafft failed"
	exit 1
fi

# if read depths included
if [ -n "$READ_DEPTHS" ]; then
	
	echo -e "converting read depth file to .csv\n"

	# replace spaces with commas in read depths file
	sed 's/ /,/g' ${READ_DEPTHS} > ${OUTDIR}/tmp.reorder.read_depth.csv

	if [ ! -f ${OUTDIR}/tmp.reorder.read_depth.csv ] ; then
		echo -e "failed to convert read depths file to .csv"
		exit 1
	fi

	# overwrite read depths variable
	READ_DEPTHS="${OUTDIR}/tmp.reorder.read_depth.csv"
			
fi

# remove any existing reordered_contigs list
if [ -f ${OUTDIR}/list_reordered_contigs.txt ] ; then
	rm ${OUTDIR}/list_reordered_contigs.txt
fi

# Get number of hyphens at beginning and end of reference

# remove first header line
sed '1d' $MAFFT > ${OUTDIR}/hyph1;

# remove all lines from second header line
header=`grep -m1 -n -i ">" ${OUTDIR}/hyph1 | cut -d : -f 1`;
head -n $((${header}-1)) ${OUTDIR}/hyph1 > ${OUTDIR}/hyph2;

# count expected number of nucleotides
#number_characters=`wc -m < ${OUTDIR}/hyph2`
#number_newlines=`wc -l < ${OUTDIR}/hyph2`
#expected_nucleotides=$((${number_characters} - ${number_newlines}))

#echo -e "\nsupplied mitochondrial reference has ${expected_nucleotides} nucleotides\n"

# split all characters onto new lines
sed 's/*\?/\n/g' ${OUTDIR}/hyph2 > ${OUTDIR}/hyph3

# remove blank lines
sed '/^[[:space:]]*$/d' ${OUTDIR}/hyph3 > ${OUTDIR}/hyph4

# find line number of first non-hyphen
non_hyph=`grep -m1 -n -v "-" ${OUTDIR}/hyph4 | cut -d : -f 1`

# Get number of hyphens at start of reference
leading_hyph_ref=$((non_hyph - 1))

echo -e "\nThe number of leading hyphens on the reference is: $leading_hyph_ref\n"

# set iterate to true
iterate="true"

# set iterate counter to 0
iterate_counter=0

# set master_reordered_contig_exists to false
master_reordered_contig_exists="false"

while [ "${iterate}" = "true" ]; do

	# invert nucleotide order
	tac ${OUTDIR}/hyph4 > ${OUTDIR}/hyph5

	# find line number of first non-hyphen
	non_hyph=`grep -m1 -n -v "-" ${OUTDIR}/hyph5 | cut -d : -f 1`

	# Get number of hyphens at end of reference
	trailing_hyph_ref=$((non_hyph - 1))

	echo -e "\nThe number of trailing hyphens on the reference is: $trailing_hyph_ref\n"

	
	# if circular is no -> in circular, causes to abort early instead of treating as if aligns right
	
	if [ "${CIRCULAR}" = "no" ]; then

	# check if leading and trailing on reference are 0 -> all contigs encompassed by reference
	if [ "${leading_hyph_ref}" = 0 ] && [ "${trailing_hyph_ref}" = 0 ]; then

		echo -e "\nall contigs are encompassed by reference, no reordering required\n"	

		# check if master reordered contig exists
		if [ ${master_reordered_contig_exists} = "false" ]; then

			# remove intermediate files
			rm ${OUTDIR}/hyph*

			# remove reordered depths file
			if [ -n "$READ_DEPTHS" ]; then

				# remove tmp.reorder files
				rm ${OUTDIR}/tmp.reorder.*

			fi

			echo -e "Renaming ${OUTDIR}/initial_mafft.fasta, to ${OUTDIR}/final_mafft.fasta\n"
			mv ${OUTDIR}/initial_mafft.fasta ${OUTDIR}/final_mafft.fasta

			# Exit script 
			exit 0

		fi

		if [ ${master_reordered_contig_exists} = "true" ]; then
			echo -e "Renaming ${OUTDIR}/realigned_mafft.fasta, to ${OUTDIR}/final_mafft.fasta\n"
			mv ${OUTDIR}/realigned_mafft.fasta ${OUTDIR}/final_mafft.fasta

			# remove intermediate files
			rm ${OUTDIR}/hyph*

			# if read depths given
			if [ -n "$READ_DEPTHS" ]; then
				echo -e "Renaming ${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv, to ${OUTDIR}/reordered_read_depths.csv\n"
				mv "${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv" ${OUTDIR}/reordered_read_depths.csv
				
				# remove intemediate files
				rm ${OUTDIR}/tmp.reorder.*
			fi

			# Exit script 
			exit 0
		fi


	fi
	# close if not circular
	fi

	# remove intermediate files
	rm -r ${OUTDIR}/hyph*

	# Remove any existing realigned.fasta
	if ls ${OUTDIR} | grep -q realigned.fasta ; then
		rm ${OUTDIR}/realigned.fasta
	fi

	# remove first header line -> tmp files in loop stored in CHecK CoNtig, chkcn
	sed '1d' $MAFFT > ${OUTDIR}/no_ref_head;

	# store remaining headers in string
	header_lines=$(printf -- '%s\n' | grep -n -i ">" ${OUTDIR}/no_ref_head | tr '\n' ' ')

	# convert string to array
	header_lines=($header_lines)

	# set reordered_contig_exists to false
	reordered_contig_exists="false"

	# loop over contig header array
	for c in "${header_lines[@]}" ;
	do
	
		# isolate contig sequence

		# store contig name
		header_text=`echo $c | cut -d : -f 2`

		echo -e "\nchecking contig: $header_text\n"

		if [ -n "$READ_DEPTHS" ]; then

			# store contig name without > marker
			header_text_unmarked=`echo "$header_text" | cut -d ">" -f 2`

			# remove _R_ prefix from header_text_unmarked if existing
			if echo ${header_text} | grep -q ">_R_" ; then
				header_text_unmarked="${header_text##>_R_}"
			fi
		
			# locate line containing header in header change file -> including comma stops contig 1 matching to contig 15
			header_change_line=`grep -i "${header_text_unmarked}," "${HEADER_MATCH}"`

			# get contig name in read depth file
			read_depth_contig_name=`echo $header_change_line | cut -d , -f 4`

			echo -e "\ngetting header specific read depth lines\n"

			# store header specific read depths in tmp file
			grep "${read_depth_contig_name}," ${READ_DEPTHS} > ${OUTDIR}/tmp.reorder.${header_text_unmarked}.read_depth.csv

			if [ ! -f ${OUTDIR}/tmp.reorder.${header_text_unmarked}.read_depth.csv ] ; then
				echo -e "failed to create header specific read depths .csv"
				exit 1
			fi

		fi

		# store contig header line number
		header_line_number=`echo $c | cut -d : -f 1`

		# get tail lines from header
		tail -n "+$header_line_number" ${OUTDIR}/no_ref_head > ${OUTDIR}/chkcn2

		# find line number of next header
		scndheader=`grep -n -i -m2 ">" ${OUTDIR}/chkcn2 | cut -d : -f 1 | tail -n1`

		# check if line number of second header equals first line
		# if no, then take head lines of file up to before second header
		if [ ! $scndheader = 1 ]; then
			head -n $((${scndheader}-1)) ${OUTDIR}/chkcn2 > ${OUTDIR}/chkcn3;
		fi

		# if yes, then send file 2 to 3 for naming consistency
		if [ $scndheader = 1 ]; then
			mv ${OUTDIR}/chkcn2 ${OUTDIR}/chkcn3
		fi

		# remove header line
		sed '1d' ${OUTDIR}/chkcn3 > ${OUTDIR}/chkcn3a;

		# count number of nucleotides on query
		number_hyphens=`tr -dc '-' < "${OUTDIR}/chkcn3a" | awk '{ print length; }'`
		number_characters=`wc -m < ${OUTDIR}/chkcn3a`
		number_newlines=`wc -l < ${OUTDIR}/chkcn3a`
		observed_nucleotides=$((${number_characters} - ${number_newlines} - ${number_hyphens}))

		echo -e "\nthe contig, $header_text, has ${observed_nucleotides} nucleotides\n"

		# split all characters onto new lines
		sed 's/*\?/\n/g' ${OUTDIR}/chkcn3a > ${OUTDIR}/chkcn4
	
		# remove blank lines
		sed '/^[[:space:]]*$/d' ${OUTDIR}/chkcn4 > ${OUTDIR}/chkcn5

		# get number of hyphens at start and end of contig

		# find line number of first non-hyphen
		non_hyph=`grep -m1 -n -v "-" ${OUTDIR}/chkcn5 | cut -d : -f 1`

		# Get number of hyphens at start of contig
		leading_hyph=$((non_hyph - 1))

		echo -e "\nThe number of leading hyphens on the contig is: $leading_hyph\n"

		# invert nucleotide order
		tac ${OUTDIR}/chkcn5 > ${OUTDIR}/chkcn6
	
		# find line number of first non-hyphen
		non_hyph=`grep -m1 -n -v "-" ${OUTDIR}/chkcn6 | cut -d : -f 1`

		# Get number of hyphens at end of contig
		trailing_hyph=$((non_hyph - 1))

		echo -e "\nThe number of trailing hyphens on the contig is: $trailing_hyph\n"

		# Isolate contig from leading and trailing hyphens

		# remove leading hyphens
		tail -n "+$((${leading_hyph}+1))" ${OUTDIR}/chkcn5 > ${OUTDIR}/chkcn7a

		# invert
		tac ${OUTDIR}/chkcn7a > ${OUTDIR}/chkcn7b

		# remove trailing hyphens
		tail -n "+$((${trailing_hyph}+1))" ${OUTDIR}/chkcn7b > ${OUTDIR}/chkcn7_inverted

		# restore order
		tac ${OUTDIR}/chkcn7_inverted > ${OUTDIR}/chkcn7

		# Check if hyphens at start of query equal to or greater than number of hyphens at start of reference

		if [ ${leading_hyph} -ge ${leading_hyph_ref} ]; then
			leading_ge="true"
		else
			leading_ge="false"
		fi

		# Check if hyphens at end of query equal to or greater than number of hyphens at end of reference

		if [ ${trailing_hyph} -ge ${trailing_hyph_ref} ]; then
			trailing_ge="true"
		else
			trailing_ge="false"
		fi
	
		# If start equal to or greater AND end equal to or greater	
		if [ ${leading_ge} = "true" ] && [ ${trailing_ge} = "true" ]; then

			# if circular and start is greater than, not equal to
			if [ ${CIRCULAR} = "yes" ] && [ ${leading_hyph} -gt ${leading_hyph_ref} ]; then
				# treat as if aligns to right of reference; move number of nucleotides as there are hyphens at start of query
				echo -e "\ncontig is circular but does not start with reference, treating as if contig aligns to right of reference\n"
					
				# subtract number of leading hyphens on reference from number of hyphens on query
				leading_diff=$((${leading_hyph} - ${leading_hyph_ref}))

				# get same number of nucleotides from end of query as is the difference
				head -n "${leading_diff}" ${OUTDIR}/chkcn7_inverted > ${OUTDIR}/file_trailing_inverted

				# invert inverted file
				tac ${OUTDIR}/file_trailing_inverted > ${OUTDIR}/file_trailing

				# get remaining nucleotides from start of query
				tail -n "+$((${leading_diff}+1))" ${OUTDIR}/chkcn7_inverted > ${OUTDIR}/file_leading_inverted
	
				# invert inverted file
				tac ${OUTDIR}/file_leading_inverted > ${OUTDIR}/file_leading

			fi

			# if not circular or leading equal, assume contig encompassed by reference - does not solve misalignment when end of query aligns to end of reference but is actually start
			if [ ${CIRCULAR} = "no" ] || [ ${leading_hyph} = ${leading_hyph_ref} ]; then
				# contig is encompassed by reference, append to realigned fasta
				echo -e "\ncontig is circular and starts with reference or is encompassed by reference\n"	

				# remove hyphens from isolated contig file
				sed 's/-//g' ${OUTDIR}/chkcn3 > ${OUTDIR}/contig_file
		
				# concatenate to realigned.fasta
				cat ${OUTDIR}/contig_file >> ${OUTDIR}/realigned.fasta

				# add newline to end of realigned fasta if doesn't exist to ensure concatenated contig headers start on newline  -> https://unix.stackexchange.com/questions/31947/how-to-add-a-newline-to-the-end-of-a-file
				sed -i -e '$a\' ${OUTDIR}/realigned.fasta

				# remove intermediate files

				rm -r ${OUTDIR}/contig_file ${OUTDIR}/chkcn*

				# if read depths given, rename unchanged file to ensure capture when compiling depths together
				if [ -n ${READ_DEPTHS} ] ; then
					mv ${OUTDIR}/tmp.reorder.${header_text_unmarked}.read_depth.csv ${OUTDIR}/tmp.reorder.${header_text_unmarked}.reordered.read_depth.csv
				fi

				# continue to next contig
				continue
			fi

		fi
	
		# If start NOT equal to or greater AND end equal to or greater	
		if [ ${leading_ge} = "false" ] && [ ${trailing_ge} = "true" ]; then
			# contig aligns to the left of the reference
			echo -e "\ncontig aligns to left of reference\n"

			# subtract number of hyphens on query from number of hyphens on reference
			leading_diff=$((${leading_hyph_ref} - ${leading_hyph}))

			# get same number of nucleotides from start of query as is the difference
			head -n "${leading_diff}" ${OUTDIR}/chkcn7 > ${OUTDIR}/file_leading

			# get remaining nucleotides from end of query
			tail -n "+$((${leading_diff}+1))" ${OUTDIR}/chkcn7 > ${OUTDIR}/file_trailing

			# if not circular, calculate number of Ns to join reordered contigs
			if [ ${CIRCULAR} = "no" ]; then
				# subtract number of hyphens at end of reference from number of hyphens at end of query
				# num_Ns=$((${trailing_hyph} - ${trailing_hyph_ref}))

				#subtract number of nucleotides observed on query from number of nucleotides included in reference
				num_Ns=$((${expected_nucleotides} - ${observed_nucleotides}))
			fi
	
		fi

		# If start equal to or greater AND end NOT equal to or greater	
		if [ ${leading_ge} = "true" ] && [ ${trailing_ge} = "false" ]; then
			# check if start is equal to -> contig is longer than reference
			if [ ${leading_hyph} = ${leading_hyph_ref} ]; then
				echo -e "\ncontig starts at reference, but extends past reference\n"

				# remove hyphens from isolated contig file
				sed 's/-//g' ${OUTDIR}/chkcn3 > ${OUTDIR}/contig_file

				# concatenate to realigned.fasta
	                        cat ${OUTDIR}/contig_file >> ${OUTDIR}/realigned.fasta

				# add newline to end of realigned fasta if doesn't exist to ensure concatenated contig headers start on newline  -> https://unix.stackexchange.com/questions/31947/how-to-add-a-newline-to-the-end-of-a-file
				sed -i -e '$a\' ${OUTDIR}/realigned.fasta

				# remove intermediate files
				rm -r ${OUTDIR}/contig_file ${OUTDIR}/chkcn*

				# if read depths given, rename unchanged file to ensure capture when compiling depths together
				if [ -n ${READ_DEPTHS} ] ; then
					mv ${OUTDIR}/tmp.reorder.${header_text_unmarked}.read_depth.csv ${OUTDIR}/tmp.reorder.${header_text_unmarked}.reordered.read_depth.csv
				fi

				# continue to next contig
				continue
			fi

			if [ ! ${leading_hyph} = ${leading_hyph_ref} ]; then
				# contig aligns to the right of the reference	
				echo -e "\ncontig aligns to right of reference\n"

				# subtract number of hyphens on query from number of hyphens on reference
				trailing_diff=$((${trailing_hyph_ref} - ${trailing_hyph}))

				# get same number of nucleotides from end of query as is the difference
				head -n "${trailing_diff}" ${OUTDIR}/chkcn7_inverted > ${OUTDIR}/file_trailing_inverted

				# invert inverted file
				tac ${OUTDIR}/file_trailing_inverted > ${OUTDIR}/file_trailing

				# get remaining nucleotides from start of query
				tail -n "+$((${trailing_diff}+1))" ${OUTDIR}/chkcn7_inverted > ${OUTDIR}/file_leading_inverted
	
				# invert inverted file
				tac ${OUTDIR}/file_leading_inverted > ${OUTDIR}/file_leading

				# if not circular, calculate number of Ns to join reordered contigs
				if [ ${CIRCULAR} = "no" ]; then
					# subtract number of hyphens at start of reference from number of hyphens at start of query
					# num_Ns=$((${leading_hyph} - ${leading_hyph_ref}))

					#subtract number of nucleotides observed on query from number of nucleotides included in reference
                                	num_Ns=$((${expected_nucleotides} - ${observed_nucleotides}))
				fi
			fi

		fi

		# If start NOT equal to or greater AND end NOT equal to or greater
		if [ ${leading_ge} = "false" ] && [ ${trailing_ge} = "false" ]; then
			# contig (or alignment) is longer than reference, move excess start to end (treat as if aligns to left)
			echo -e "\ncontig or alignment is longer than reference, treating as if aligns to left of reference\n"

			# subtract number of hyphens on query from number of hyphens on reference
			leading_diff=$((${leading_hyph_ref} - ${leading_hyph}))

			# get same number of nucleotides from start of query as is the difference
			head -n "${leading_diff}" ${OUTDIR}/chkcn7 > ${OUTDIR}/file_leading

			# get remaining nucleotides from end of query
			tail -n "+$((${leading_diff}+1))" ${OUTDIR}/chkcn7 > ${OUTDIR}/file_trailing

			# if not circular, join with 3 Ns
			if [ ${CIRCULAR} = "no" ]; then
				num_Ns=3
			fi

		fi

		# set reordered_contig_exists to true
		reordered_contig_exists="true"

		# set master_reordered_contig_exists to true
		master_reordered_contig_exists="true"

		# if not circular, create Ns file

		if [ ${CIRCULAR} = "no" ]; then
			
			# if num_Ns is negative, correct to 3
			if [ ${num_Ns} -lt 0 ] ; then
				num_Ns=3
			fi

			# for num_Ns iterations, add N to Ns file
			echo -e "\ncontig is not circular, will rejoin contig ends with ${num_Ns} Ns\n"
			for ((i=1;i<=${num_Ns};i++)); do
    				printf '%s' "N" >> ${OUTDIR}/file_Ns
			done

		fi

		# reorder nucleotides

		# store nucleotides on single line

		# convert newline characters to space
		file_trailingA=`tr '\n' ' ' < ${OUTDIR}/file_trailing`
		file_leadingA=`tr '\n' ' ' < ${OUTDIR}/file_leading`

		# convert to file
		echo $file_trailingA > ${OUTDIR}/file_trailingA
		echo $file_leadingA > ${OUTDIR}/file_leadingA

		# remove spaces
		sed 's/ //g' ${OUTDIR}/file_trailingA > ${OUTDIR}/file_trailingB
		sed 's/ //g' ${OUTDIR}/file_leadingA > ${OUTDIR}/file_leadingB

		# remove hyphens
		sed 's/-//g' ${OUTDIR}/file_trailingB > ${OUTDIR}/file_trailingC
		sed 's/-//g' ${OUTDIR}/file_leadingB > ${OUTDIR}/file_leadingC

                # check if reordered_contigs_list exists
		if [ -f ${OUTDIR}/list_reordered_contigs.txt ]; then
			# look for contig header in reordered list, $ should anchor pattern to end of line and ^ to the beginning, preventing mismatches
			if grep -qi '^${header_text}$' ${OUTDIR}/list_reordered_contigs.txt ; then
				echo -e "\ncontig has been reordered previously, will not rejoin using Ns if otherwise applicable\n"
				contig_reordered="true"
			else
				contig_reordered="false"
			fi
		else
			# if reordered_contigs_list.txt does not exist, set contig_reordered to false
			contig_reordered="false"
		fi

# if read depths included, count number of nucleotides in file_leadingC -> get head of read depth lines matching number of nucleotides; get tail of read depth lines
# if file_Ns exist, create tmp read depth file with contig,bp_1,bp_2,read_depth equal to contig,NA,NA,0
# cat tail, Ns, head together into new tmp header file -> need to append to combined file at some point
		if [ -n "$READ_DEPTHS" ]; then

			# count number of nucleotides moved from start of sequence
			num_nucleotides=`wc -m < ${OUTDIR}/file_leadingC`

			# get equal number of head lines from read depths file
			head -n ${num_nucleotides} ${OUTDIR}/tmp.reorder.${header_text_unmarked}.read_depth.csv > ${OUTDIR}/tmp.reorder.${header_text_unmarked}.leading.read_depth.csv

			# get remaining from read depths file
			tail -n "+$((${num_nucleotides}+1))" ${OUTDIR}/tmp.reorder.${header_text_unmarked}.read_depth.csv > ${OUTDIR}/tmp.reorder.${header_text_unmarked}.trailing.read_depth.csv

			# if Ns file does not exist, concatenate read depth files in new order
			if [ ! -f ${OUTDIR}/file_Ns ] || [ "${contig_reordered}" = "true" ]; then
				cat ${OUTDIR}/tmp.reorder.${header_text_unmarked}.trailing.read_depth.csv ${OUTDIR}/tmp.reorder.${header_text_unmarked}.leading.read_depth.csv > ${OUTDIR}/tmp.reorder.${header_text_unmarked}.reordered.read_depth.csv
			fi

			# if Ns file exists AND contig not reordered previously, create read depths Ns file and concatenate read depth files in new order
			if [ -f ${OUTDIR}/file_Ns ]; then
				if [ "${contig_reordered}" = "false" ] ; then
					# for num_Ns iterations, add contig,NA,NA,0 to read depths Ns file
					for ((i=1;i<=${num_Ns};i++)); do
    						printf '%s\n' "${read_depth_contig_name},NA,NA,0" >> ${OUTDIR}/tmp.reorder.${header_text_unmarked}.Ns.read_depth.csv
					done

					# create reordered file
					cat ${OUTDIR}/tmp.reorder.${header_text_unmarked}.trailing.read_depth.csv ${OUTDIR}/tmp.reorder.${header_text_unmarked}.Ns.read_depth.csv ${OUTDIR}/tmp.reorder.${header_text_unmarked}.leading.read_depth.csv > ${OUTDIR}/tmp.reorder.${header_text_unmarked}.reordered.read_depth.csv
			
				fi
			fi
		fi
# CURRENT EDITING ENDS HERE #

		# store header in file
		echo -e "$header_text" > ${OUTDIR}/file_header

		# reorder contig and remove newline

		# if Ns file exists AND contig not reordered previously, add Ns file
		if [ -f ${OUTDIR}/file_Ns ] && [ "${contig_reordered}" = "false" ] ; then
			cat ${OUTDIR}/file_trailingC ${OUTDIR}/file_Ns ${OUTDIR}/file_leadingC > ${OUTDIR}/reordered_contig | tr '\n' ' '

			#add header to reordered_contigs_list.txt
			printf '%s\n' "${header_text}" >> ${OUTDIR}/list_reordered_contigs.txt

		fi

		# if Ns file does not exist OR contig reordered previously, exclude Ns file
		if [ ! -f ${OUTDIR}/file_Ns ] || [ "${contig_reordered}" = "true" ] ; then
			cat ${OUTDIR}/file_trailingC ${OUTDIR}/file_leadingC > ${OUTDIR}/reordered_contig | tr '\n' ' '

                        #add header to reordered_contigs_list.txt
			printf '%s\n' "${header_text}" >> ${OUTDIR}/list_reordered_contigs.txt

		fi

		# remove space introduced by newline removal
		sed 's/ //g' ${OUTDIR}/reordered_contig > ${OUTDIR}/reordered_contigS

		# add header
		cat ${OUTDIR}/file_header ${OUTDIR}/reordered_contigS > ${OUTDIR}/reordered_contigH

		echo -e "\nreordered contig $header_text\n"
		
		# remove intermediate files
		rm -r ${OUTDIR}/chkcn*
		rm -r ${OUTDIR}/file_*

		# append to fasta file of other contigs (starting with reference file)

		cat ${OUTDIR}/reordered_contigH >> ${OUTDIR}/realigned.fasta 

		# add newline to end of file to ensure contig headers start on newline
		sed -i -e '$a\' ${OUTDIR}/realigned.fasta
	
		# remove individual reordered contig file and intermediate files
		rm -r ${OUTDIR}/reordered_contig*
		
	# end contig loop
	done

	# remove intermediate files
	rm -r ${OUTDIR}/no_ref_head

	# remove blank lines from reordered fasta
	sed -i '/^[[:space:]]*$/d' ${OUTDIR}/realigned.fasta

	# check if no reordered contigs exist
	if [ ${reordered_contig_exists} = "false" ]; then

		echo -e "\nno reordering required\n"

		# check if master reordered contig exists
		if [ ${master_reordered_contig_exists} = "false" ]; then
			# remove realigned.fasta
			rm -r ${OUTDIR}/realigned.fasta
			
			# remove reordered depths file
			if [ -n "$READ_DEPTHS" ]; then
				# remove tmp.reorder files
				rm ${OUTDIR}/tmp.reorder.*
			fi

			echo -e "Renaming ${OUTDIR}/initial_mafft.fasta, to ${OUTDIR}/final_mafft.fasta\n"
			mv ${OUTDIR}/initial_mafft.fasta ${OUTDIR}/final_mafft.fasta

			# Exit script 
			exit 0
		fi

		if [ ${master_reordered_contig_exists} = "true" ]; then
			echo -e "Renaming ${OUTDIR}/realigned_mafft.fasta, to ${OUTDIR}/final_mafft.fasta\n"
			mv ${OUTDIR}/realigned_mafft.fasta ${OUTDIR}/final_mafft.fasta

			# if read depths given
			if [ -n "$READ_DEPTHS" ]; then
				echo -e "Renaming ${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv, to ${OUTDIR}/reordered_read_depths.csv\n"
				mv "${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv" ${OUTDIR}/reordered_read_depths.csv
			
				if [ -f ${OUTDIR}/tmp.reorder.* ]; then	
					# remove tmp.reorder files
					rm ${OUTDIR}/tmp.reorder.*
				fi
			fi

			# Exit script 
			exit 0
		fi

	fi

	if [ -n "$READ_DEPTHS" ]; then

		# if combined read depths file exists, remove it
		if [ -f "${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv" ]; then
			rm "${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv"
		fi	

		# if read depth file given, concat individual reordered contig read depth files into single file
		echo -e "Combining reordered depths files into single file: ${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv\n"
		cat ${OUTDIR}/tmp.reorder.*.reordered.read_depth.csv > ${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv

		# overwrite READ_DEPTHS with new base file
		READ_DEPTHS="${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv"

		# remove tmp.reorder files
		if [ -f ${OUTDIR}/tmp.reorder.* ]; then
                        # remove tmp.reorder files
                	rm ${OUTDIR}/tmp.reorder.*
                fi

	fi

	# Run mafft on realigned.fasta
	echo -e "\nrunning mafft on realigned.fasta\n"

	# run mafft against single reference and do not allow adjust direction
	$PWD/external_software/mafft-linux64/mafft.bat --6merpair --thread -1 --addfragments ${OUTDIR}/realigned.fasta $MT_REFERENCE > ${OUTDIR}/realigned_mafft.fasta

	echo -e "\nmafft complete: realigned_mafft.fasta\n"

	# Check if there are leading hyphens on the reference

	MAFFT=${OUTDIR}/realigned_mafft.fasta

	# remove first header line
	sed '1d' $MAFFT > ${OUTDIR}/hyph1;

	# remove all lines from second header line
	header=`grep -m1 -n -i ">" ${OUTDIR}/hyph1 | cut -d : -f 1`;
	head -n $((${header}-1)) ${OUTDIR}/hyph1 > ${OUTDIR}/hyph2;

	# split all characters onto new lines
	sed 's/*\?/\n/g' ${OUTDIR}/hyph2 > ${OUTDIR}/hyph3

	# remove blank lines
	sed '/^[[:space:]]*$/d' ${OUTDIR}/hyph3 > ${OUTDIR}/hyph4

	# find line number of first non-hyphen
	non_hyph=`grep -m1 -n -v "-" ${OUTDIR}/hyph4 | cut -d : -f 1`

	# Get number of hyphens at start of reference
	leading_hyph_ref=$((non_hyph - 1))

	echo -e "\nThe number of leading hyphens on the reference is: $leading_hyph_ref\n"

	# if there are no leading hyphens, check if a contig starts with a non-hyphen
	if [ "${leading_hyph_ref}" = 0 ]; then

		# check that a contig starts with a non-hyphen
		
		# remove first header line -> tmp files in loop stored in CHecK CoNtig, chkcn
		sed '1d' $MAFFT > ${OUTDIR}/chkcn1;

		# store remaining headers in string
		header_lines=$(printf -- '%s\n' | grep -n -i ">" ${OUTDIR}/chkcn1 | tr '\n' ' ')

		# convert string to array
		header_lines=($header_lines)

		# loop over contig header array
		for c in "${header_lines[@]}" ;
		do
	
			# isolate contig sequence
	
			# store contig name
			header_text=`echo $c | cut -d : -f 2`

			# store contig header line number
			header_line_number=`echo $c | cut -d : -f 1`

			# get tail lines from header
			tail -n "+$header_line_number" ${OUTDIR}/chkcn1 > ${OUTDIR}/chkcn2

			# remove header line
			sed '1d' ${OUTDIR}/chkcn2 > ${OUTDIR}/chkcn3a;

			# get first line of nucleotides
			head -n1 ${OUTDIR}/chkcn3a > ${OUTDIR}/chkcn3b

			# split all characters onto new lines
			sed 's/*\?/\n/g' ${OUTDIR}/chkcn3b > ${OUTDIR}/chkcn4
	
			# remove blank lines
			sed '/^[[:space:]]*$/d' ${OUTDIR}/chkcn4 > ${OUTDIR}/chkcn5
			
			# find line number of first non-hyphen
			non_hyph=`grep -m1 -n -v "-" ${OUTDIR}/chkcn5 | cut -d : -f 1`

			# if line number equals 1, set iterate to false and break
			if [ "${non_hyph}" = 1 ]; then
				echo -e "\ncontig ${header_text} starts with reference\n"
				iterate="false"
				break
			fi
		
		done

		if [ "${iterate}" = "true" ]; then
			echo -e "\nno contig starts with reference\n"
		fi

		# remove intermediate files
		rm -r ${OUTDIR}/chkcn*

	fi

	# check if iterate counter has reached max
	if [ "${iterate_counter}" = "${MAX_I}" ]; then
		iterate="false"
		echo -e "\nWARNING: reached maximum number of iterations: ${MAX_I}\n"
		printf "WARNING: reached maximum number of iterations: ${MAX_I}" >> ${OUTDIR}/WARNING_reorder_nucleotides.txt
	fi

	if [ "${iterate}" = "true" ]; then
		iterate_counter=$((iterate_counter+1))
		echo -e "\nITERATION ${iterate_counter}\n"
	fi

# close while loop
done

echo -e "Renaming ${OUTDIR}/realigned_mafft.fasta, to ${OUTDIR}/final_mafft.fasta\n"
mv ${OUTDIR}/realigned_mafft.fasta ${OUTDIR}/final_mafft.fasta

if [ -n "$READ_DEPTHS" ]; then
	# rename tmp reordered read depths file
	mv "${OUTDIR}/keep_tmp.reorder.combined.reordered.read_depth.csv" ${OUTDIR}/reordered_read_depths.csv
fi

# remove intermediate files
rm -r ${OUTDIR}/hyph*

