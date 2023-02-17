#!/usr/bin/bash
echo -n "Enter Deliverables Directory PATH: "
read deliverables_dir
#deliverables_dir=/storage/colddata/basesolve/pipeline/NCGM_Methylseq_DD_09_01_2023/OUTPUT

echo -n "Enter Deliverables Directory PATH at scratch: "
read scratch_deliverables_dir
#scratch_deliverables_dir=/storage/scratch2/NCGM-DrParth_DD_28_12_2022

echo -n "Enter SymbolicLink to fastq files Directory PATH: "
read symDir
#symDir=/storage/scratch2/NCGM-DrParth_DD_28_12_2022/00_rawdata

echo -n "Enter Slurm Partition: MidP OR HighP: "
read partition
#partition="MidP"

for i in `ls ${symDir} | grep R1`;do echo $i;done
echo -n "Please look into the read names and provide the extension after sample name Example: _R1_001.fastq.gz: "
read "R1_extn"
#R1_extn="_R1_001.fastq.gz"

for i in `ls ${symDir} | grep R2`;do echo $i;done
echo -n "Please look into the read names and provide the extension after sample name Example: _R2_001.fastq.gz: "
read "R2_extn"
#R2_extn="_R2_001.fastq.gz"

#echo -n "number of threads: "
#read cpu
#cpu=256

echo -n "Enter reference file: "
read "ref"
#ref=/storage/scratch2/NCGM-DrParth_DD_28_12_2022/DATABASE/grch38_core+bs_controls.fa

## Coldata storage

echo -n "Enter the directory containing all the containers: "
read "tools"


################
#Defining tools
################
#tools="/storage/colddata/basesolve/tools/emseq-tools"
fastp="${tools}/fastp_0.23.2.sif fastp"
bwa_meth="${tools}/bwameth_0.2.2.sif bwameth.py"
mark_NC="${tools}/mark-nonconverted-reads.sif mark-nonconverted-reads.py"
sambamba="${tools}/sambamba_0.8.2.sif sambamba"
samtools="${tools}/samtools_1_16_1.sif samtools"
fastqc="${tools}/fastqc_0.11.9.sif fastqc"
picard="${tools}/picard_2.27.4.sif picard"
goleft="${tools}/goleft_0.2.4.sif goleft"
samblaster="${tools}/samblaster_0.1.20.sif samblaster"


touch ${scratch_deliverables_dir}/LOG_FILE.txt
LOG=${scratch_deliverables_dir}/LOG_FILE.txt

SAMPLES_DONE=()

		######################################################################################################

#						STEP 2: RUN FASTP ON THE RAW FASTQ FILES

		######################################################################################################


declare -A fastp_jobs

for i in $(ls $symDir | grep ${R1_extn}); do
	sample=`basename $i ${R1_extn}`

	# Enter the input Directory
	
	mkdir -p ${sample}
	mkdir -p $scratch_deliverables_dir/${sample}/00_logs
	mkdir -p $scratch_deliverables_dir/${sample}/02_PROCESSED_DATA

	logpath="${scratch_deliverables_dir}/${sample}/00_logs"

	R1_in=${symDir}/${sample}${R1_extn}
	R2_in=${symDir}/${sample}${R2_extn}

	# Enter the Output Directory

	R1_out=${scratch_deliverables_dir}/${sample}/02_PROCESSED_DATA/${sample}_R1.fastp.fastq.gz
	R2_out=${scratch_deliverables_dir}/${sample}/02_PROCESSED_DATA/${sample}_R2.fastp.fastq.gz

	echo "Start the Analysis of ${sample}" >> ${LOG}
	echo ${R1_in} >> ${LOG}
	echo ${R2_in} >> ${LOG}

	html_report="${scratch_deliverables_dir}/${sample}/02_PROCESSED_DATA/${sample}_fastp.html"
	json_report="${scratch_deliverables_dir}/${sample}/02_PROCESSED_DATA/${sample}_fastp.json"

	inst_name=$(zcat ${R1_in} | head -n 1 | cut -f 1 -d ':' | sed 's/^@//')
	echo $inst_name >> ${LOG}    
	fastq_barcode=$(zcat ${R1_in} | head -n 1 | sed -r 's/.*://')
	echo $fastq_barcode >> ${LOG}
	    if [[ "${inst_name:0:2}" == 'A0' ]] || [[ "${inst_name:0:2}" == 'NS' ]] || \
       		[[ "${inst_name:0:2}" == 'NB' ]] || [[ "${inst_name:0:2}" == 'VH' ]] ; then
	       trim_polyg='--trim_poly_g'
	       echo '2-color instrument: poly-g trim mode on' >> ${LOG}
	    else
	       trim_polyg=''
	    fi
	
	echo "Start the FASTP Analysis of ${sample}" >> ${LOG}
	
	# Execute the Slurm Script

	run_fastp=$(sbatch --parsable \
	--mem=10G \
	--cpus-per-task=16 \
	--partition=${partition} \
	-J fastp_${sample} \
	-o ${logpath}/${sample}_fastp.out \
	-e ${logpath}/${sample}_fastp.err \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${fastp} --in1 ${R1_in} --in2 ${R2_in} --out1 ${R1_out} --out2 ${R2_out} --thread 16 -l 2 -Q ${trim_polyg} --split=16  --html ${html_report} --json ${json_report} --report_title ${sample}_fastp_report 2>&1 | tee ${logpath}/${sample}_fastp.log && touch ${logpath}/${sample}_fastp.ok")
	
	# Check the Status of the Job
	
	ST="PENDING"
	while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
		ST=$(sacct -j ${run_fastp##* } -o State | awk 'FNR == 3 {print $1}')
		EXIT=$(sacct -j ${run_fastp##* } -o ExitCode | awk 'FNR == 3 {print $1}')
		sleep 10s
	done
	
	echo "FINAL STATUS $ST" >> ${LOG}
	echo "EXIT CODE $EXIT" >> ${LOG}

	if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
		echo "FastP Failed" >> ${LOG}
		exit 1
	elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
		echo "Successfully Running!!" >> ${LOG}
		fastp_jobs[$sample]=${run_fastp}
	elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
		echo "Successfully Completed!!" >> ${LOG}
		fastp_jobs[$sample]=${run_fastp}
	fi
	
	echo "Completed FastP for ${sample}" >> ${LOG}
done

	
		#############################################################################################################

#						STEP 3: RUN BWAMETH ON THE FASTP PROCESSED FILES							

		#############################################################################################################




for keys in "${!fastp_jobs[@]}"; do
	echo "Associative for sample ${keys} is ${fastp_jobs[$keys]}" >> ${LOG}
done


declare -A bwameth_jobs

for i in $(ls $symDir | grep ${R1_extn}); do
        sample=`basename $i ${R1_extn}`
	echo "Starting to run BwaMeth for ${sample}" >> ${LOG}

	# Get the information about read groups

	header=$(zcat ${R1_in} | head -n 1)
	id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
	sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
#	echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"

	# Input Directory:
	
	FASTP_DIR=${scratch_deliverables_dir}/${sample}/02_PROCESSED_DATA/	

	# Create Output Directories
	
	mkdir -p ${scratch_deliverables_dir}/${sample}/03_BWAMETH
	mkdir -p ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES
	BWAOUT=${scratch_deliverables_dir}/${sample}/03_BWAMETH

	logpath="${scratch_deliverables_dir}/${sample}/00_logs"

	# Get the FastP file extensions #

	R1_fastp_extn=_R1.fastp.fastq.gz
        R2_fastp_extn=_R2.fastp.fastq.gz
	
	logpath="${scratch_deliverables_dir}/${sample}/00_logs"

	 for j in `ls ${FASTP_DIR} | grep ${R1_fastp_extn}`; do
		SUBSAMPLE=`basename $j ${R1_fastp_extn}`
		
		FASTP_R1_in=${scratch_deliverables_dir}/${sample}/02_PROCESSED_DATA/${SUBSAMPLE}${R1_fastp_extn}
                FASTP_R2_in=${scratch_deliverables_dir}/${sample}/02_PROCESSED_DATA/${SUBSAMPLE}${R2_fastp_extn}

		# Enter the Commands to be executed

		echo "${bwa_meth} -p --threads 16 --read-group $(echo \"@RG\\tID:$id\\tSM:$id"_"$sm\\tLB:$id"_"$sm\\tPL:ILLUMINA\") --reference ${ref} ${FASTP_R1_in} ${FASTP_R2_in} > ${scratch_deliverables_dir}/${sample}/03_BWAMETH/${SUBSAMPLE}.bwameth.sam" > ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/${SUBSAMPLE}.sh
        done
	
#	Execute the Slurm Script	

	run_bwameth=$(sbatch --parsable -d afterok:${fastp_jobs[$sample]} \
	--mem=50G \
	--cpus-per-task=256 \
	--partition=${partition} \
	-J ${sample}_bwameth \
	-o ${logpath}/${sample}_bwameth.out \
	-e ${logpath}/${sample}_bwameth.log \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" parallel bash ::: ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0001.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0002.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0003.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0004.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0005.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0006.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0007.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0008.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0009.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0010.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0011.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0012.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0013.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0014.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0015.${sample}.sh ${scratch_deliverables_dir}/${sample}/BWAMETH_FILES/0016.${sample}.sh 2>&1 | tee ${logpath}/bwameth.log && touch ${logpath}/${sample}_BWAMETH.ok")


#	Check the Status of the Job

	ST="PENDING"
	while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "CANCELLED+" && "$ST" != "RUNNING"  ]] ; do
		ST=$(sacct -j ${run_bwameth##* } -o State | awk 'FNR == 3 {print $1}')
		EXIT=$(sacct -j ${run_bwameth##* } -o ExitCode | awk 'FNR == 3 {print $1}')
		sleep 10s
	done
	
	echo "FINAL STATUS $ST" >> ${LOG}
	echo "EXIT CODE $EXIT" >> ${LOG}
	if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
		echo "Bwa Meth Failed" >> ${LOG}
		exit 1
	elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
		echo "Successfully Completed!!" >> ${LOG}
		bwameth_jobs[$sample]=${run_bwameth}
	elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
		echo "Successfully Running!!!" >> ${LOG}
		bwameth_jobs[$sample]=${run_bwameth}
	fi
	
	echo "Executed jobs for Bwa-Meth for Samples: ${sample}" >> ${LOG}

done

		######################################################################################################################

	#				STEP 4: RUN MARK NONCONVERTED READS ON BWAMETH ALIGNED SAM FILES

		######################################################################################################################

declare -A marknc_jobs

for keys in "${!bwameth_jobs[@]}"; do
        echo "Associative for sample ${keys} is ${bwameth_jobs[$keys]}" >> ${LOG}
done

for i in $(ls $symDir | grep ${R1_extn}); do
        sample=`basename $i ${R1_extn}`
	echo "Now running Mark Nonconverted on the Bwameth output files in ${sample}" >> ${LOG}

	# Enter Input Directory:

	BWAMETH_DIR="${scratch_deliverables_dir}/${sample}/03_BWAMETH"

	sam_extn=".bwameth.sam"

	## Create Output Directory:

	mkdir -p ${scratch_deliverables_dir}/${sample}/04_MARKNC
        mkdir -p ${scratch_deliverables_dir}/${sample}/MARKNC_FILES

	MARKNC=${scratch_deliverables_dir}/${sample}/04_MARKNC
	logpath="${scratch_deliverables_dir}/${sample}/00_logs"

	for j in `ls ${BWAMETH_DIR} | grep ${sam_extn}`; do
                SUBSAMPLE=`basename $j ${sam_extn}`

		SAM_IN=${BWAMETH_DIR}/${SUBSAMPLE}${sam_extn}

		# Enter the Commands to be executed		

		echo "${mark_NC} --reference ${ref} --bam ${SAM_IN} --out ${MARKNC}/${SUBSAMPLE}_nonconverted.sam  2> ${MARKNC}/${SUBSAMPLE}_nonconverted.tsv" > ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/${SUBSAMPLE}.sh
	done
	
	# Execute the Slurm Script

	run_marknc=$(sbatch --parsable -d afterok:${bwameth_jobs[$sample]} \
	--mem=100G \
        --cpus-per-task=256 \
        --partition=${partition} \
        -J ${sample}_markNC \
        -o ${logpath}/${sample}_markNC.out \
        -e ${logpath}/${sample}_markNC.log \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" parallel bash ::: ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0001.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0002.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0003.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0004.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0005.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0006.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0007.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0008.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0009.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0010.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0011.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0012.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0013.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0014.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0015.${sample}.sh ${scratch_deliverables_dir}/${sample}/MARKNC_FILES/0016.${sample}.sh 2>&1 | tee ${logpath}/marknc.log && touch ${logpath}/${sample}_MARKNC.ok && rm -fr ${BWAMETH_DIR}")

## 	Check the Status of the Job ##


	ST="PENDING"
	while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
		ST=$(sacct -j ${run_marknc##* } -o State | awk 'FNR == 3 {print $1}')
		EXIT=$(sacct -j ${run_marknc##* } -o ExitCode | awk 'FNR == 3 {print $1}')
		sleep 10s
	done

        echo "FINAL STATUS $ST" >> ${LOG}
	echo "EXIT CODE $EXIT" >> ${LOG}
	if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
		echo "Mark Nonconverted Failed" >> ${LOG}
		exit 1
	elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
		echo "Successfully Running Mark Nonconverted Reads" >> ${LOG}
		marknc_jobs[$sample]=${run_marknc}
	elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
		echo "Successfully completed Mark Nonconverted Reads" >> ${LOG}
		marknc_jobs[$sample]=${run_marknc}
	fi

	echo "Running Mark NonConverted reads for ${sample}" >> ${LOG}
done

		#########################################################################################################################
#		STEP 5 CONVERT SAM TO BAM AND SORT THEM FOR EACH OF THE SPLITTED FILES
		#######################################################################################################################

declare -A sam2bam_jobs

for keys in "${!marknc_jobs[@]}"; do
        echo "Associative for sample ${keys} is ${marknc_jobs[$keys]}" >> ${LOG}
done

for i in $(ls $symDir | grep ${R1_extn}); do
        sample=`basename $i ${R1_extn}`
        echo "Now Running Sambama to Convert Sam to Bam on the MarkNC output files in ${sample}" >> ${LOG}

        # Enter Input Directory:

	MARKNC=${scratch_deliverables_dir}/${sample}/04_MARKNC
	
	marknc_extn="_nonconverted.sam"

        ## Create Output Directory:

	mkdir -p ${scratch_deliverables_dir}/${sample}/05_SAMBAMBA_OUT
        mkdir -p ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES	
	
	SAMBAMBA=${scratch_deliverables_dir}/${sample}/05_SAMBAMBA_OUT

	logpath="${scratch_deliverables_dir}/${sample}/00_logs"
	
	for j in `ls ${MARKNC} | grep ${marknc_extn}`; do
                SUBSAMPLE=`basename $j ${marknc_extn}`
                MARKNC_IN=${MARKNC}/${SUBSAMPLE}${marknc_extn}
                
		# Enter the Commands to be executed

		echo "${scratch_deliverables_dir}/tools/sambamba_0.8.2.sif sambamba view -t 16 -S ${MARKNC_IN} -f bam -o ${SAMBAMBA}/${SUBSAMPLE}.bam" > ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/${SUBSAMPLE}.sh		

                echo "${sambamba} sort -m 20GB -t 16 --tmpdir=${scratch_deliverables_dir}/${sample} --sort-by-name -o ${SAMBAMBA}/${SUBSAMPLE}_sorted.bam ${SAMBAMBA}/${SUBSAMPLE}.bam" >> ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/${SUBSAMPLE}.sh

                echo "rm ${SAMBAMBA}/${SUBSAMPLE}.bam" >> ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/${SUBSAMPLE}.sh
	done

	# Execute the Slurm Script		

	run_sambamba=$(sbatch --parsable -d afterok:${marknc_jobs[$sample]} \
	--mem=320G \
	--cpus-per-task=256 \
	--partition=${partition} \
	-J ${sample}_Sam2Bam_splitted \
	-o ${logpath}/${sample}_Sam2Bam_splitted.out \
	-e ${logpath}/${sample}_Sam2Bam_splitted.log \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" parallel bash ::: ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0001.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0002.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0003.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0004.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0005.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0006.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0007.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0008.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0009.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0010.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0011.${sample}.sh  ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0012.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0013.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0014.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0015.${sample}.sh ${scratch_deliverables_dir}/${sample}/SAMBAMBA_FILES/0016.${sample}.sh 2>&1 | tee ${logpath}/01_splitted_sam2bam.log && rm ${MARKNC}/*.sam")
	
#	Check the Status of the Job

	ST="PENDING"
	while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
		ST=$(sacct -j ${run_sambamba##* } -o State | awk 'FNR == 3 {print $1}')
		EXIT=$(sacct -j ${run_sambamba##* } -o ExitCode | awk 'FNR == 3 {print $1}')
		sleep 10s
	done

	echo "FINAL STATUS $ST" >> ${LOG}
	echo "EXIT CODE $EXIT" >> ${LOG}
	if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
		echo "SAM 2 BAM Failed" >> ${LOG}
		exit 1
	elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
		echo "Successfully Running SAM 2 BAM" >> ${LOG}
		sam2bam_jobs[$sample]=${run_sambamba}
	elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
		echo "Successfully completed SAM 2 BAM" >> ${LOG}
		sam2bam_jobs[$sample]=${run_sambamba}
	fi
	
	echo "Completed SAM 2 BAM for Splitted BAM files for ${sample}" >> ${LOG}
done

		##############################################################################################
#						STEP 6: Merge the BAM Files
		##############################################################################################

declare -A mergebam_jobs

for keys in "${!sam2bam_jobs[@]}"; do
	echo "Associative for sample ${keys} is ${sam2bam_jobs[$keys]}" >> ${LOG}
done

for i in $(ls $symDir | grep ${R1_extn}); do
	sample=`basename $i ${R1_extn}`
	echo "Now Merging the Splitted Bam Files for ${sample}" >> ${LOG}
	
	# Enter Input Directory:
	
	SAMBAMBA_IN=${scratch_deliverables_dir}/${sample}/05_SAMBAMBA_OUT
	
	# Enter Output Directory:

	mkdir -p ${scratch_deliverables_dir}/${sample}/05_MERGED_BAM
	MERGED_BAM=${scratch_deliverables_dir}/${sample}/05_MERGED_BAM
	
	logpath="${scratch_deliverables_dir}/${sample}/00_logs"

	# Enter the Commands to be executed

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/sambamba_0.8.2.sif sambamba merge ${MERGED_BAM}/${sample}.bam ${SAMBAMBA_IN}/0001.${sample}_sorted.bam ${SAMBAMBA_IN}/0002.${sample}_sorted.bam ${SAMBAMBA_IN}/0003.${sample}_sorted.bam ${SAMBAMBA_IN}/0004.${sample}_sorted.bam ${SAMBAMBA_IN}/0005.${sample}_sorted.bam ${SAMBAMBA_IN}/0006.${sample}_sorted.bam ${SAMBAMBA_IN}/0007.${sample}_sorted.bam ${SAMBAMBA_IN}/0008.${sample}_sorted.bam ${SAMBAMBA_IN}/0009.${sample}_sorted.bam ${SAMBAMBA_IN}/0010.${sample}_sorted.bam ${SAMBAMBA_IN}/0011.${sample}_sorted.bam ${SAMBAMBA_IN}/0012.${sample}_sorted.bam ${SAMBAMBA_IN}/0013.${sample}_sorted.bam ${SAMBAMBA_IN}/0014.${sample}_sorted.bam ${SAMBAMBA_IN}/0015.${sample}_sorted.bam ${SAMBAMBA_IN}/0016.${sample}_sorted.bam -t 256 2>&1 | tee ${logpath}/02_splitted_merge.log" >> ${MERGED_BAM}/${sample}_Merge.sh

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/sambamba_0.8.2.sif sambamba sort -m 200GB -t 64 --tmpdir=${scratch_deliverables_dir}/${sample} --sort-by-name -o ${MERGED_BAM}/${sample}_sorted.bam ${MERGED_BAM}/${sample}.bam 2>&1 | tee ${logpath}/03_sort_merged.log" >> ${MERGED_BAM}/${sample}_Merge.sh
	
	## Execute the Slurm Script:

	run_mergebam=$(sbatch --parsable -d afterok:${sam2bam_jobs[$sample]} \
	--mem=200G \
	--cpus-per-task=64 \
	--partition=${partition} \
	-J ${sample}_MergeBam_splitted \
	-o ${logpath}/${sample}_MergeBam.out \
	-e ${logpath}/${sample}_MergeBam.log \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" bash ${MERGED_BAM}/${sample}_Merge.sh && rm -fr ${SAMBAMBA_IN}")

	# Check the Status of the job

	ST="PENDING"
	while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
		ST=$(sacct -j ${run_mergebam##* } -o State | awk 'FNR == 3 {print $1}')
		EXIT=$(sacct -j ${run_mergebam##* } -o ExitCode | awk 'FNR == 3 {print $1}')
		sleep 10s
	done

	echo "FINAL STATUS $ST" >> ${LOG}
	echo "EXIT CODE $EXIT" >> ${LOG}
	if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
		echo "MergeBam Failed" >> ${LOG}
		exit 1
	elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
		echo "Successfully Running Merge BAM" >> ${LOG}
		mergebam_jobs[$sample]=${run_mergebam}
	elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
		echo "Successfully completed Merge BAM" >> ${LOG}
		mergebam_jobs[$sample]=${run_mergebam}
	fi

        echo "Completed Merging BAM files for ${sample}" >> ${LOG}
	
#	rm ${MERGED_BAM}/${sample}_Merge.sh
done
	###################################################################################################

#	STEP 7 : MARK DUPLICATES ON THE MERGED BAM FILES

	##################################################################################################

declare -A markdup_jobs

for keys in "${!mergebam_jobs[@]}"; do
	echo "Associative for sample ${keys} is ${mergebam_jobs[$keys]}" >> ${LOG}
done

for i in $(ls $symDir | grep ${R1_extn}); do
	sample=`basename $i ${R1_extn}`
	echo "Now Running Mark Duplicates on  ${sample}" >> ${LOG}

	# Enter Input Directory:
	
	MERGED_BAM_IN=${scratch_deliverables_dir}/${sample}/05_MERGED_BAM

	bam_extn="_sorted.bam"
	Bam_in=${MERGED_BAM_IN}/${sample}${bam_extn}

	# Enter Output Directory:

	mkdir -p ${scratch_deliverables_dir}/${sample}/06_MERGE_MARK_DUPLICATES
	MARKDUP_DIR=${scratch_deliverables_dir}/${sample}/06_MERGE_MARK_DUPLICATES
	DUPOUT=${scratch_deliverables_dir}/${sample}/06_MERGE_MARK_DUPLICATES/${sample}.sam
	BAMDUPOUT=${scratch_deliverables_dir}/${sample}/06_MERGE_MARK_DUPLICATES/${sample}.bam
	SORTOUT=${scratch_deliverables_dir}/${sample}/06_MERGE_MARK_DUPLICATES/${sample}_SORTED.bam

	logpath="${scratch_deliverables_dir}/${sample}/00_logs"

	# Enter the Commands to be executed

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/samtools_1_16_1.sif samtools cat -b <(find ${Bam_in}) | samtools view -h | ${scratch_deliverables_dir}/tools/samblaster_0.1.20.sif samblaster -o ${DUPOUT} 2>&1 | tee ${logpath}/01_markdup_samblaster.log" >> ${MARKDUP_DIR}/${sample}_MarkDuplicates.sh

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/sambamba_0.8.2.sif sambamba view -t 64 -S -f bam ${DUPOUT} -o ${BAMDUPOUT} 2>&1 | tee ${logpath}/02_markdup_sam2bam.log" >> ${MARKDUP_DIR}/${sample}_MarkDuplicates.sh

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/sambamba_0.8.2.sif sambamba sort -m 200GB --tmpdir=${scratch_deliverables_dir}/${sample} -t 64 -o ${SORTOUT} ${BAMDUPOUT} 2>&1 | tee ${logpath}/03_markdup_sambamba_sort.log" >> ${MARKDUP_DIR}/${sample}_MarkDuplicates.sh
	
	echo "rm ${DUPOUT}" >> ${MARKDUP_DIR}/${sample}_MarkDuplicates.sh

	echo "rm ${BAMDUPOUT}" >>  ${MARKDUP_DIR}/${sample}_MarkDuplicates.sh 

	## Execute the Slurm Script
	
	run_markdup=$(sbatch --parsable -d afterok:${mergebam_jobs[$sample]} \
	--mem=200G \
	--cpus-per-task=64 \
	--partition=${partition} \
	-J ${sample}_MarkDuplicates \
	-o ${logpath}/${sample}_MarkDuplicates.out \
	-e ${logpath}/${sample}_MarkDuplicates.log \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" bash ${MARKDUP_DIR}/${sample}_MarkDuplicates.sh && rm -fr ${MERGED_BAM_IN}")
	
	# Check the Status of the job

	ST="PENDING"
	while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
		ST=$(sacct -j ${run_markdup##* } -o State | awk 'FNR == 3 {print $1}')
		EXIT=$(sacct -j ${run_markdup##* } -o ExitCode | awk 'FNR == 3 {print $1}')
		sleep 10s
	done

        echo "FINAL STATUS $ST" >> ${LOG}
        echo "EXIT CODE $EXIT" >> ${LOG}
        if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
                echo "MarkDuplicates Failed" >> ${LOG}
                exit 1
        elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
                echo "Successfully Running Mark Duplicates" >> ${LOG}
                markdup_jobs[$sample]=${run_markdup}
        elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
                echo "Successfully completed Mark Duplicates" >> ${LOG}
                markdup_jobs[$sample]=${run_markdup}
        fi

	echo "Completed Mark Duplicates for ${sample}" >> ${LOG}

done


		##################################################################################################################

#			STEP 8 RUN METHYL DACKEL ON THE FINAL SORTED BAM FILES AFTER MARK DUPLICATION

		#################################################################################################################


declare -A methyldackel_jobs
declare -A process_human_jobs
declare -A bamfastqc_jobs
declare -A samstats_jobs
declare -A picard_jobs
declare -A goleft_jobs


for keys in "${!markdup_jobs[@]}"; do
	echo "Associative for sample ${keys} is ${markdup_jobs[$keys]}" >> ${LOG}
done

for i in $(ls $symDir | grep ${R1_extn}); do
	sample=`basename $i ${R1_extn}`
	echo "Now Running Methyl Dackel on  ${sample}" >> ${LOG}

	# Enter Input Directory:

	MARKDUP_DIR_IN=${scratch_deliverables_dir}/${sample}/06_MERGE_MARK_DUPLICATES
	
	BAM_extn="_SORTED.bam"
	
	BAM_in=${MARKDUP_DIR_IN}/${sample}${BAM_extn}

	# Enter Output Directory:

	mkdir -p ${scratch_deliverables_dir}/${sample}/07_METHYL_DACKEL

	COMMAND=${scratch_deliverables_dir}/${sample}/07_METHYL_DACKEL/${sample}_METHACKEL.sh
	METHACK_OUT=${scratch_deliverables_dir}/${sample}/07_METHYL_DACKEL/${sample}_mbias.tsv

	logpath="${scratch_deliverables_dir}/${sample}/00_logs"

	# Enter the Commands
	echo "echo -e \"chr\\tcontext\\tstrand\\tRead\\tPosition\\tnMethylated\\tnUnmethylated\\tnMethylated(+dups)\\tnUnmethylated(+dups)\" > ${METHACK_OUT}" >> ${COMMAND}
	echo "chrs=(\`samtools view -H ${BAM_in} | grep @SQ | cut -f2 | sed 's/SN://' | grep -v _random | grep -v chrUn | sed 's/|/\\|/'\`)" >> ${COMMAND}
	echo "for chr in \${chrs[*]}; do" >> ${COMMAND}
	echo -e "\tfor context in CHH CHG CpG; do" >> ${COMMAND}
	echo -e "\t\tjoin -t \$'\\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \\" >> ${COMMAND}
	echo -e "\t\t<( \\" >> ${COMMAND} 
	echo -e "\t\t/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/methyldackel_0_6_1.sif MethylDackel mbias --noSVG \$arg -@ 64 -r \$chr ${ref} ${BAM_in} ${METHACK_OUT} | \\" >> ${COMMAND}
	echo -e "\t\ttail -n +2 | awk '{print \$1\"-\"\$2\"-\"\$3\"\\t\"\$0}' | sort -k 1b,1" >> ${COMMAND}
	echo -e "\t\t) \\" >> ${COMMAND}
	echo -e "\t\t<( \\" >> ${COMMAND}
	echo -e "\t\t/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/methyldackel_0_6_1.sif MethylDackel mbias --noSVG --keepDupes \$arg -F 2816 -@ 64 -r \$chr ${ref} ${BAM_in} ${METHACK_OUT} | \\" >> ${COMMAND}
	echo -e "\t\ttail -n +2 | awk '{print \$1\"-\"\$2\"-\"\$3\"\\t\"\$0}' | sort -k 1b,1" >> ${COMMAND}
	echo -e "\t\t) \\" >> ${COMMAND}
	echo -e "\t\t| sed \"s/^/\${chr}\\t\${context}\\t/\" \\" >> ${COMMAND}
	echo -e "\t\t>> ${METHACK_OUT}" >> ${COMMAND}
	echo -e "\tdone" >> ${COMMAND}
	echo "done" >> ${COMMAND}
	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/methyldackel_0_6_1.sif MethylDackel mbias -@ 64 --noCpG --CHH --CHG -r \${chrs[0]} ${ref} ${BAM_in} ${METHACK_OUT}_chn 2>&1 | tee ${logpath}/01_methyldackel.log" >> ${COMMAND}

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/methyldackel_0_6_1.sif MethylDackel mbias -@ 64 -r \${chrs[0]} ${ref} ${BAM_in} ${METHACK_OUT}_cpg 2>&1 | tee ${logpath}/02_methyldackel.log" >> ${COMMAND}

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/methyldackel_0_6_1.sif MethylDackel extract --methylKit --OT 0,0,0,95 --OB 0,0,5,0 -@ 64 --CHH --CHG  -o ${METHACK_OUT}_EXTRACTED ${ref} ${BAM_in} 2>&1 | tee ${logpath}/03_methyldackel.log" >> ${COMMAND}

	
#	Execute the Slurm Script

	run_methyldackel=$(sbatch --parsable -d afterok:${markdup_jobs[$sample]} \
        --mem=96G \
        --cpus-per-task=64 \
	--partition=${partition} \
	-J ${sample}_MethylDackel \
	-o ${logpath}/${sample}_MethylDackel.out \
	-e ${logpath}/${sample}_MethylDackel.log \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" bash ${COMMAND}")

# 	Check the Status of the job

	ST="PENDING"
	while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
		ST=$(sacct -j ${run_methyldackel##* } -o State | awk 'FNR == 3 {print $1}')
		EXIT=$(sacct -j ${run_methyldackel##* } -o ExitCode | awk 'FNR == 3 {print $1}')
		sleep 10s
	done

	echo "FINAL STATUS $ST" >> ${LOG}
	echo "EXIT CODE $EXIT" >> ${LOG}
	if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
		echo "Methyl Dackel Failed" >> ${LOG}
		exit 1
	elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
		echo "Successfully Running Methyl Dackel" >> ${LOG}
		methyldackel_jobs[$sample]=${run_methyldackel}
	elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
		echo "Successfully completed Methyl Dackel" >> ${LOG}
		methyldackel_jobs[$sample]=${run_methyldackel}
	fi

	echo "Completed Running Methyl Dackel for Sample: ${sample}" >> ${LOG}

		##############################################################################################################

#						STEP 9: RUN PROCESS HUMAN READS

		##############################################################################################################


	for keys in "${!methyldackel_jobs[@]}"; do
		echo "Associative for sample ${keys} is ${methyldackel_jobs[$keys]}" >> ${LOG}
	done
	
	echo "Now running Process Human reads for ${sample}" >> ${LOG}
	
	
#	Enter the output Directory

	mkdir -p ${scratch_deliverables_dir}/${sample}/08_PROCESS_HUMAN_READS
	PROCESS_HUMAN_OUT=${scratch_deliverables_dir}/${sample}/08_PROCESS_HUMAN_READS

	logpath="${scratch_deliverables_dir}/${sample}/00_logs"

#	Execute the Slurm Script
	
	run_process_human=$(sbatch --parsable -d afterok:${markdup_jobs[$sample]} \
	--mem=32G \
	--cpus-per-task=64 \
	--partition=${partition} \
	-J ${sample}_Process_Human_reads \
	-o ${logpath}/${sample}_Process_Human_reads.out \
	-e ${logpath}/${sample}_Process_Human_reads.log \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/sambamba_0.8.2.sif sambamba view -t 64 -l 0 -f bam ${BAM_in} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${PROCESS_HUMAN_OUT}/${sample}_human.bam 2>&1 | tee ${logpath}/process_human_reads.log && touch ${logpath}/${sample}_process_humanreads.ok")

#	Check the status of the Job

	ST="PENDING"
        while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
                ST=$(sacct -j ${run_process_human##* } -o State | awk 'FNR == 3 {print $1}')
                EXIT=$(sacct -j ${run_process_human##* } -o ExitCode | awk 'FNR == 3 {print $1}')
                sleep 10s
        done

        echo "FINAL STATUS $ST" >> ${LOG}
        echo "EXIT CODE $EXIT" >> ${LOG}
        if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
                echo "Process Human Reads Failed" >> ${LOG}
                exit 1
        elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
                echo "Successfully Running Process Human Reads" >> ${LOG}
                process_human_jobs[$sample]=${run_process_human}
        elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
                echo "Successfully completed Process Human Reads" >> ${LOG}
                process_human_jobs[$sample]=${run_process_human}
        fi

	echo "Completed Running Process Human Reads" >> ${LOG}

		#############################################################################################################################

#						STEP 10		RUN FASTQC OF BAM FILES

		############################################################################################################################# 

	for keys in "${!process_human_jobs[@]}"; do
		echo "Associative for sample ${keys} is ${process_human_jobs[$keys]}" >> ${LOG}
	done

	echo "Now running fastQC on BAM for ${sample}" >> ${LOG}
	
#       Enter the output Directory

	mkdir -p ${scratch_deliverables_dir}/${sample}/09_FASTQC_BAM
	FASTQC_BAM=${scratch_deliverables_dir}/${sample}/09_FASTQC_BAM

	run_bamfastqc=$(sbatch --parsable -d afterok:${markdup_jobs[$sample]} \
	--mem=32G \
	--cpus-per-task=16 \
	--partition=${partition} \
	-J ${sample}_Bam_FastQC \
	-o ${logpath}/${sample}_Bam_FastQC.out \
	-e ${logpath}/${sample}_Bam_FastQC.log \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/fastqc_0_11_9.sif fastqc -f bam ${BAM_in} -o ${FASTQC_BAM} 2>&1 | tee ${logpath}/fastqc_bam.log && touch ${logpath}/${sample}_bamfastqc.ok")

#	Check the status of the Job

	ST="PENDING"
        while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
                ST=$(sacct -j ${run_bamfastqc##* } -o State | awk 'FNR == 3 {print $1}')
                EXIT=$(sacct -j ${run_bamfastqc##* } -o ExitCode | awk 'FNR == 3 {print $1}')
                sleep 10s
        done

        echo "FINAL STATUS $ST" >> ${LOG}
        echo "EXIT CODE $EXIT" >> ${LOG}
        if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
                echo "BAM FastQC Failed" >> ${LOG}
                exit 1
        elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
                echo "Successfully Running Fastqc on BAM" >> ${LOG}
                bamfastqc_jobs[$sample]=${run_bamfastqc}
        elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
                echo "Successfully completed Fastqc on BAM" >> ${LOG}
                bamfastqc_jobs[$sample]=${run_bamfastqc}
        fi

        echo "Completed Running BAM FastQC for ${sample}" >> ${LOG}

		#######################################################################################################################################

#							STEP 11	CALCULATE SAM STATS FROM BAM FILES

		#######################################################################################################################################

	echo "Now Running Calculate SamStats for ${sample}" >> ${LOG}

#	Enter the output Directory

	mkdir -p ${scratch_deliverables_dir}/${sample}/10_SAM_STATS
	SAM_STATS=${scratch_deliverables_dir}/${sample}/10_SAM_STATS

#	Enter the Commands to be executed

	COMMAND_11=${SAM_STATS}/${sample}_SAMSTATS.sh

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/samtools_1_16_1.sif samtools flagstats -@ 16 ${BAM_in} > ${SAM_STATS}/${sample}.flagstats 2>&1 | tee ${logpath}/01_flagstats.log" >> ${COMMAND_11}
	
	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/samtools_1_16_1.sif samtools idxstats -@ 16 ${BAM_in} > ${SAM_STATS}/${sample}.idxstats 2>&1 | tee ${logpath}/02_idxstats.log" >> ${COMMAND_11}

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/samtools_1_16_1.sif samtools stats -@ 16 ${BAM_in} > ${SAM_STATS}/${sample}.stats 2>&1 | tee ${logpath}/03_stats.log" >> ${COMMAND_11}

# 	Execute the Slurm Script

	run_samstats=$(sbatch --parsable -d afterok:${markdup_jobs[$sample]} \
        --mem=32G \
        --cpus-per-task=16 \
        --partition=${partition} \
        -J ${sample}_Samstats \
        -o ${logpath}/${sample}_Samstats.out \
        -e ${logpath}/${sample}_Samstats.log \
	--wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" bash ${COMMAND_11} && touch ${logpath}/${sample}_samstats.ok")

#	Check the status of the job

	ST="PENDING"
        while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
                ST=$(sacct -j ${run_samstats##* } -o State | awk 'FNR == 3 {print $1}')
                EXIT=$(sacct -j ${run_samstats##* } -o ExitCode | awk 'FNR == 3 {print $1}')
                sleep 10s
        done

        echo "FINAL STATUS $ST" >> ${LOG}
        echo "EXIT CODE $EXIT" >> ${LOG}
        if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
                echo "Samstats Failed" >> ${LOG}
                exit 1
        elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
                echo "Successfully Running Sam stats" >> ${LOG}
                samstats_jobs[$sample]=${run_samstats}
        elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
                echo "Successfully completed Sam stats" >> ${LOG}
                samstats_jobs[$sample]=${run_samstats}
        fi

	echo "Completed running Calculate Sam stats for ${sample}" >> ${LOG}
	
		#####################################################################################################################

#						STEP 12: RUN PICARD FOR GC BIAS AND INSERT SIZE

		#####################################################################################################################

	echo "Now Running Picard to calculate GC Bias and Insert size for Sample:$sample" >> ${LOG}

#	Enter the output Directory

	mkdir -p ${scratch_deliverables_dir}/${sample}/11_PICARD_GCBIAS_IS
	PICARD=${scratch_deliverables_dir}/${sample}/11_PICARD_GCBIAS_IS

	COMMAND_12=${PICARD}/${sample}_PICARD.sh

#	Enter the Commands

	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/Picard_2.27.4.sif picard -Xmx32g CollectGcBiasMetrics -BS true --VALIDATION_STRINGENCY LENIENT -I ${BAM_in} -O ${PICARD}/${sample}.gcmetrics -S ${PICARD}/${sample}.gc_summary_metrics -CHART ${PICARD}/${sample}.gc.pdf -R ${ref} 2>&1 | tee ${logpath}/01_picard_gcbias.log" >> ${COMMAND_12}
	echo "/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/Picard_2.27.4.sif picard -Xmx32g CollectInsertSizeMetrics --VALIDATION_STRINGENCY LENIENT -I ${BAM_in} -O ${PICARD}/${sample}.insertsize_metrics --MINIMUM_PCT 0 -H /dev/null 2>&1 | tee ${logpath}/02_picard_IS.log" >> ${COMMAND_12}

#	Execute the Slurm Script

	run_picard=$(sbatch --parsable -d afterok:${markdup_jobs[$sample]} \
	--mem=64G \
	--cpus-per-task=16 \
	--partition=${partition} \
	-J ${sample}_Picard \
	-o ${logpath}/${sample}_Picard.out \
	-e ${logpath}/${sample}_Picard.log \
        --wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" bash ${COMMAND_12} && touch ${logpath}/${sample}_picard.ok")

#       Check the status of the job

	ST="PENDING"
        while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "CANCELLED+"  ]] ; do
                ST=$(sacct -j ${run_picard##* } -o State | awk 'FNR == 3 {print $1}')
                EXIT=$(sacct -j ${run_picard##* } -o ExitCode | awk 'FNR == 3 {print $1}')
                sleep 10s
        done

        echo "FINAL STATUS $ST" >> ${LOG}
        echo "EXIT CODE $EXIT" >> ${LOG}
        if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
                echo "Picard Failed" >> ${LOG}
                exit 1
        elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
                echo "Successfully Running Picard" >> ${LOG}
                picard_jobs[$sample]=${run_picard}
        elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
                echo "Successfully completed Picard" >> ${LOG}
                picard_jobs[$sample]=${run_picard}
        fi

	echo "Completed Running Picard on ${sample}" >> ${LOG}

		######################################################################################################

#					STEP 13: RUN GOLEFT TO CALCULATE THE COVERAGE

		#####################################################################################################

	echo "Now Running GOLeft to get the Coverage for: $sample" >> ${LOG}

#	Enter the output Directory

	mkdir -p ${scratch_deliverables_dir}/${sample}/12_GOLEFT
	GOLEFT=${scratch_deliverables_dir}/${sample}/12_GOLEFT

#	Execute the Slurm Script

	run_goleft=$(sbatch --parsable -d afterok:${markdup_jobs[$sample]} \
        --mem=12G \
        --cpus-per-task=16 \
        --partition=${partition} \
        -J ${sample}_GOLEFT \
        -o ${logpath}/${sample}_GOLEFT.out \
        -e ${logpath}/${sample}_GOLEFT.log \
        --wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" /usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" ${scratch_deliverables_dir}/tools/goleft_0.2.4.sif goleft indexcov --directory ${GOLEFT}/${sample} ${BAM_in} 2>&1 | tee ${logpath}/Goleft.log  && touch ${logpath}/${sample}_goleft.ok")


#	Check the status of the Job

	ST="PENDING"
        while  [[ "$ST" != "COMPLETED" && "$ST" != "FAILED" && "$ST" != "RUNNING" && "$ST" != "CANCELLED+"  ]] ; do
                ST=$(sacct -j ${run_goleft##* } -o State | awk 'FNR == 3 {print $1}')
                EXIT=$(sacct -j ${run_goleft##* } -o ExitCode | awk 'FNR == 3 {print $1}')
                sleep 10s
        done

        echo "FINAL STATUS $ST" >> ${LOG}
        echo "EXIT CODE $EXIT" >> ${LOG}
        if [[ $ST == "COMPLETED" && $EXIT != "0:0" ]]; then
                echo "GOLeft Failed" >> ${LOG}
                exit 1
        elif [[ $ST == "RUNNING" && $EXIT == "0:0" ]]; then
                echo "Successfully Running GOLeft" >> ${LOG}
                goleft_jobs[$sample]=${run_picard}
        elif [[ $ST == "COMPLETED" && $EXIT == "0:0" ]]; then
                echo "Successfully completed GOLeft" >> ${LOG}
                goleft_jobs[$sample]=${run_picard}
        fi

        echo "Completed Running GoLEFT on ${sample}" >> ${LOG}

	rm -fr ${Scratch_deliverables_dir}/${sample}/sambamba-pid*

	## Move the Sample Folder to the Deliverables Directory
	
	run_move=$(sbatch --parsable -d afterok:${picard_jobs[$sample]} \
        --mem=10G \
        --cpus-per-task=10 \
        --partition=${partition} \
        -J Move_${sample} \
        -o ${sample}_move.out \
        -e ${sample}_move.log \
        --wrap="/usr/bin/time -f \"mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U\" mv ${scratch_deliverables_dir}/${sample} ${deliverables_dir}")

done


echo "################# COMPLETED EM-SEQ PIPELINE EXECUTION ############################"
#############################################################################################################################################

#						EMSEQ PIPELINE EXECUTED

############################################################################################################################################
