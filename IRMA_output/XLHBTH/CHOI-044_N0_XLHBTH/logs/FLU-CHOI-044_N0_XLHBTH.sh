### BACKGROUND INFO ###
# Iterative Refinement Meta-Assembler (IRMA), v1.2.0, 23 Aug 2024
# Samuel S. Shepard (CDC/DDID/NCIRD/ID), vfn4@cdc.gov
# Last commit: e72b7a9d51ab9ad023c6fe69f53dcc95aac92510

### START CONFIG ###
PARAM_FILE_NAME="FLU"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2015-01-29"

### PERFORMANCE ###
GRID_ON=0	# grid computation on [1,0] for on or off
GRID_PATH=""	# optional grid path for working directory, leave "" for no option
LIMIT_BLAT=60000	# threshold before grid
LIMIT_SSW=80000		# threshold before grid
LIMIT_SAM=500		# threshold before grid
SINGLE_LOCAL_PROC=16	# local maximum processes
DOUBLE_LOCAL_PROC=8	# local maximum processes (double this number)
ALLOW_TMP=1		# if GRID_ON=0, try to use /tmp for working directory
TMP=/tmp		# the scratch/tmpfs for working on the assemblies
PACKAGED_FASTQ=1		# Create *.tar.gz for final fastq, otherwise use *.gz for each file
USE_IRMA_CORE=0	# Use IRMA CORE where available for improved performance

### REFERENCE ###
MIN_FA=1		# no alternative reference [0..1]
MIN_CA=20		# minimum count for alternative finished assembly
SKIP_E=1		# skip reference elongation
REF_SET=/Users/ashleysobelleonard/code/flu-amd/IRMA_RES/modules/FLU/reference/consensus.fasta	# Starting reference, usually default for $DEF_SET
ASSEM_REF=0

### READ GATHERING ###
MAX_ROUNDS=5		# round of read gathering
USE_MEDIAN=1		# use the median quality or the average [1,0]
QUAL_THRESHOLD=30	# minimum read statistic
MIN_LEN=125		# minimum read length
INCL_CHIM=0		# includes chimera or not [0,1]
			# transposase adapter, clips 5' of the adapter on the forward strand and 3' on the reverse strand
			# applicable to nexttera pair-end reads 
ADAPTER="AGATGTGTATAAGAGACAG"
FUZZY_ADAPTER=1			# Allow 1 adapter mismatch
ENFORCE_CLIPPED_LENGTH=0	# Reads are filtered for minimum length post adapter trimming.	
MERGE_SECONDARY=0	# merge secondary data after the first round, good if no co-infections
BLAT_IDENTITY=80	# blat identity for read gathering, usually 80
DO_SECONDARY=0		# do secondary assembly
RESIDUAL_ASSEMBLY_FACTOR=0	# the primary match to alternative match (secondary) ratio limit. If less than this number, do the residual.

## MATCH STEP
MATCH_PROC=20		# grid maximum processes for the MATCH
MATCH_PROG="BLAT"	# match (all or any match) program [BLAT]
MIN_RP=15		# minimum read pattern count to continue doing primary assembly
MIN_RC=15		# minimum read count to continue doing primary assembly
MIN_RP_RESIDUAL=150		# minimum read pattern count to continue doing residual assembly
MIN_RC_RESIDUAL=150		# minimum read count to continue doing residual assembly
MIN_BLAT_MATCH=0	# minimum blat match, default settings within the program practically limit to 30 bp, only useful if set higher.

## SORT STEP 
SORT_PROG="BLAT"	# [LABEL,BLAT]
SORT_PROC=80		# currently not used
NONSEGMENTED=0		# segmented! [0,1]
# LABEL
LFASTM=		# LABEL sorting fast-mode
GENE_GROUP="HA,NA:OG"

## ALIGN STEP ##
ALIGN_PROG="SAM"	# rough assembly / alignment to working reference [SAM,BLAT]
ALIGN_PROC=20		# grid maximum processes for the rough align
DEL_TYPE=""	# rough assembly deletion handling: NNN, REF, or DEL

### FINISHING ASSEMBLY ###
MAX_ITER_ASSEM=5	# max assembly iteration [5]
NO_MERGE=0		# do not merge read pairs [0]
ASSEM_PROG="SSW"	# assembly program [SSW]
ASSEM_PROC=20		# grid maximum processes for assembly
INS_T=0.25		# minimum frquenncy threshold for insertion refinement
DEL_T=0.60		# minimum frequency threshold for deletion refinement 
INS_T_DEPTH=1		# minimum coverage depth for insertion refinement
DEL_T_DEPTH=1		# minimum coverage depth for deletion refinement 
SILENCE_COMPLEX_INDELS=0	# silences reads with complex indels (having 4 indels with at least one greater than 9)
MIN_AMBIG=0.20		# minimum called SNV frequency for mixed base in amended consensus folder
SSW_M=2			# smith-waterman match score
SSW_X=5			# smith-waterman mismatch penalty
SSW_O=10			# smith-waterman gap open penalty
SSW_E=1			# smith-waterman gap extension penalty
MM2_A=2		# minimap2 match
MM2_B=5		# minimap2 mismatch
MM2_O=10		# minimap2 gap open
MM2_E=1		# minimap2 gap extend
MIN_CONS_SUPPORT=1	# minimum consensus support for amended consensus
MIN_CONS_QUALITY=0	# minimum consensus average quality for amended consensus
ALIGN_AMENDED=0		# do global alignment of the amended consensus to the HMM profile
PADDED_CONSENSUS=0	# attempt to pad amended_consensus with Ns for amplicaton dropout: requires ALIGN_AMENDED=1 and ASSEM_REF=1
MIN_DROPOUT_EDGE_DEPTH=0 # minimum threshold before dropouts (91+) with flanking regions (6+) are masked during final iterative refinement

### VARIANT CALLING ###
# HEURISTICS
AUTO_F=1		# auto-adjust frequency threshold [1,0]
MIN_FI=0.005		# minimum insertion variant frequency
MIN_FD=0.005		# minimum deletion variant frequency
MIN_F=0.008		# minimum frequency for single nucleotide variants
MIN_C=2			# minimum count for variants
MIN_AQ=24		# minimum average variant quality, does not apply to deletions
MIN_TCC=100		# minimum non-ambiguous column coverage
MIN_CONF=0.80		# minimum confidence not machine error

# CONFIDENCE INTERVALS
SIG_LEVEL=0.999		# significance test level for variant calling (.90,.95,.99,.999). 
