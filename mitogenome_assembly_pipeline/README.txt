Pipeline for mitogenome assembly from NGS bycatch using MitoZ (MEGAHIT), MITObim (MIRA), and NOVOPlasty

This pipeline was developed specifically for our purposes.
While some effort was made to make it reusable, it was not developed to be a fully portable and supported pipeline for adaptation to other projects.

Software dependency installation instructions:

All software dependencies should be installed in a directory labelled 'external_software' under the top level project directory (the pipeline should be executed from this project directory).

BBMap (v39.06):
	Download from sourceforge:
	wget https://sourceforge.net/projects/bbmap/files/BBMap_39.06.tar.gz

	unpack with:
	tar xfvz BBMap_39.06.tar.gz

Trimmomatic (v0.39):
	install from github repository:
	wget https://github.com/usadellab/Trimmomatic/archive/refs/tags/v0.39.tar.gz

	unpack with:
	tar zxvf v0.39.tar.gz

MAFFT (v7.526):
	Download from mafft website:
	https://mafft.cbrc.jp/alignment/software/mafft-7.526-linux.tgz

	unpack with:
	tar xfvz mafft-7.526-linux.tgz


MITObim (v1.9.1):
	clone the .git repository using:
	git clone --recursive git://github.com/chrishah/MITObim.git

	you may need to make the MITObim.pl file executable with chmod a+x MITObim.pl

MitoZ (v3.6):
	We have followed the instructions on the MitoZ github issues page (https://github.com/linzhi2013/MitoZ/issues/212) to address a bug that affected our workflow.
	For convenience, the edited MitoZ docker image is supplied here as 'mitoz_v3.6u.tar.gz'.

	You will need to decompress (but not unpack) the .tar.gz archive (e.g. with gunzip)

	When the official MitoZ GitHub updates MitoZ to address this bug, the version here will be removed.

NOVOPlasty (v4.3.5):
	Download .pl file from github e.g. with:
	wget https://github.com/ndierckx/NOVOPlasty/blob/master/NOVOPlasty4.3.5.pl

	You may need to make the NOVOPlasty4.3.5.pl file executable with chmod a+x NOVOPlasty4.3.5.pl

Udocker (v1.3.16):

	Install from github repository:
	wget https://github.com/indigo-dc/udocker/releases/download/1.3.16/udocker-1.3.16.tar.gz

	unpack with:
	tar zxvf udocker-1.3.16.tar.gz

Software dependency versions and file paths have been hard-coded into the pipeline.
Use of other versions or software not installed as directed in the 'external_software' directory will require editing of the pipeline code.


Compilation of assembly statistics with the --assembly_stats option relies on R scripts.

Pipeline has been successfully run with R v4.3.3 and relies on the following packages:
	tidyverse(v2.0.0)
	sjmisc(v2.8.11)
	janitor(v2.2.0)
	

Example command:
bash auto_scripts/mitogenome_assembly_pipeline.sh \
	-r "0_raw" \
	-m "all" \
	-p \
	--max_memory 50 \
	--udocker \
	--mt_reference reference/Mus_musculus_reference.fasta \
	--mt_annotated_ref reference/Mus_musculus_annotated_reference.fasta \
	--sample_id_list sample_list.txt \
	--stats_methods_list methods_list.txt \
	--assembly_stats \
	--outprefix "test"  \
	--exclude_genes_read_depth "COX1 CYTB"
	
Details on arguments and accepted options can be shown with:
bash auto_scripts/mitogenome_assembly_pipeline.sh --help


Important notes on pipline function:
While the mitogenome_assembly_pipeline.sh file can be moved to another folder,
all .sh and .R files invoked by mitogenome_assembly_pipeline.sh must be located in an 'auto_scripts' directory under the top level project directory (the pipeline should be executed from this project directory).

Paths supplied to --mt_reference and --mt_annotated_ref must be the canonical path to the file, it cannot include symbolic links.

The file supplied to --mt_annotated_ref should be a GenBank gene_features .fasta file, ideally with all 37 loci included.
Header rows must contain MitoZ formatted loci names, you may need to manually edit some of these for successful function (e.g. replace [gene=mt-Rnr2] with [gene=s-rRNA])

For use with taxa or purposes where you do not expect annotation of the standard 37 animal mitochondrial loci, an expected genes text/csv file can be supplied.
It must follow the following format (header included):
gene,num_copies
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
trnV,1