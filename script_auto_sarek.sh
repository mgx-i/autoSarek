#!/bin/bash

# set -x
# set -e

# load variables
source auto_sarek.conf  #TODO: modify path


############################ FUNCTIONS

# create sarek_runs.txt file containing run names and status
make_runs_file () {
    # argument 1 runs ouput directory
    # argument 2 run list
    for RUN in ${$2}
        do
        #TODO: improve header, verify format
        HEADER="[HEADER]
        A list of all the runs in this directory and their status regarding the sarek workflow
            0: to ignore
            1: to be treated
            2: done
        [FILES]"
        echo "${HEADER}" >>${RUNS_FILE}
        echo "${RUN} 0" >>${$1}${RUNS_FILE}
        done
}

# get bed file from samplesheet. If does not exist, modify status to 1
get_bed_file () {
    local DESCRIPTION=$(grep "Description.*\.bed" "${SEQ_PATH}${RUN}/${SAMPLESHEET_SEQ}")
    if [ -z ${DESCRIPTION} ]
    then 
        sed -i "s/${RUN_NAME} .*/${RUN_NAME} 1/" "${SEQ_PATH}${RUNS_FILE}"
        LINE=$(grep "${RUN}" "${SEQ_PATH}${RUNS_FILE}")
    else
        local tmp=${DESCRIPTION#*,}
        BEDFILE=${tmp%#*}
        #TODO: check if bed exists in valid beds list
    fi
}

create_samplesheet () {
    # header
    echo "patient,sample,lane,fastq_1,fastq_2" >"${SEQ_PATH}${RUN}/${SAMPLESHEET_SAREK}"
    # loop on Samplesheet.csv on Data after header
    sed -n "/\[Data\]/ {n; :a; n; p; ba}" "${SEQ_PATH}${RUN}/${SAMPLESHEET_SEQ}" | \
    while read SAMPLESHEET_LINE #TODO: check if works
        do
        local SAMPLE_NAME=${SAMPLESHEET_LINE%%,*}
        local FASTQ1=$(ls ${SEQ_PATH}${RUN}/Alignement_*/${RUN}*/Fastq/ | grep "^${SAMPLE_NAME}_.*_R1_.*")
        local FASTQ2=$(ls ${SEQ_PATH}${RUN}/Alignement_*/${RUN}*/Fastq/ | grep "^${SAMPLE_NAME}_.*_R2_.*")
        echo "${SAMPLE_NAME},${SAMPLE_NAME},1,${FASTQ1},${FASTQ2}" >>"${SEQ_PATH}${RUN}/${SAMPLESHEET_SAREK}"
    done
}

launch_sarek () {
    #TODO: find which reference genomes from the bed name or exhaustive list 
    # should all be in the name except CF_panel_v2_3395481.bed (to check)
    
    # get genome version
    VGENOME=$(echo ${BEDFILE} | grep -oE 'hg[0-9]{2}')
    if [ ${VGENOME} = hg38 || ${BEDFILE} in ${BEDLIST38} ]
    then 
        # string with genome-specific parameters
        PARAMS="--genome ${GENOME_38} \
                --igenomes_base ${IGENOMES_BASE_38} \
                --fasta ${FASTA_38} \
                --fasta_fai ${FASTA_FAI_38}"
    fi
    if [ ${VGENOME} = hg19 || ${BEDFILE} in ${BEDLIST19} ]
    then 
        # string with genome-specific parameters
        PARAMS=""
    fi
    #TODO: if want to add an else to notify unknown bed, do cases ? or second in first's else
    
    #TODO: necessary to get every parameter through a variable ?
    COMMAND="nextflow run ${SAREK_PATH} -profile singularity \
        -c ${CONFIG_PATH} \
        -bg \
        -w ${WORKDIR_PATH} \
        --input "${SEQ_PATH}${RUN}/${SAMPLESHEET_SAREK}" \
        --outdir "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}" \
        --step ${STEP} \
        --tools ${TOOLS} \
        --intervals ${BEDFILE} \
        --wes ${WES} \
        --dict ${DICT} \
        --save_output_as_bam ${SAVE_OUTPUT_AS_BAM} \
        --concatenate_vcfs ${CONCATENATE_VCFS} \
        --trim_fastq ${TRIM_FASTQ} \
        --split_fastq ${SPLIT_FASTQ} \
        ${PARAMS}"
    ${COMMAND} 2>"${SEQ_PATH}${RUNS_FILE}/${NF_ERROR_FILE}" > #TODO: output file
}



############################ MAIN

for SEQ_PATH in ${MINISEQ_PATH} ${MISEQ_PATH} ${NEXTSEQ_PATH}
    do
    RUN_LIST=$(ls ${SEQ_PATH} | egrep '^[0-9]{6}_')

    # checks runs file
    if [ ! -f "${SEQ_PATH}${RUNS_FILE}" ]
    then
        make_runs_file ${SEQ_PATH} ${RUN_LIST}
    fi

    for RUN in ${RUN_LIST}
        do
        LINE=$(grep "${RUN}" "${SEQ_PATH}${RUNS_FILE}")
        # adds the new run to runs file
        if [ -z "${LINE}" ]
        then
              echo "${RUN} 1" >>"${SEQ_PATH}${RUNS_FILE}"
        fi
        get_bed_file 
        STATUS=${LINE##* }
        # checks if the run is ready (status 1 and trigger file exists)
        if [ STATUS = 1 && -f ${SEQ_PATH}${RUNS_FILE}/${TRIGGER_FILE} ]
        then
            create_samplesheet
            conda activate NF-core
            launch_sarek
            # check if run executed properly and change status to 2
            # check for Pipeline completed successfully when 

            conda deactivate
        done

done

