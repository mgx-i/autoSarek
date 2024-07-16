#!/bin/bash

set -e

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

change_status () {
    # argument : new status
    sed -i "s/${RUN_NAME} .*/${RUN_NAME} $1/" "${SEQ_PATH}${RUNS_FILE}"
}

# get bed file from samplesheet. If does not exist, modify status to 0 (ignore)
get_bed_file () {
    local DESCRIPTION=$(grep "Description.*\.bed" "${SEQ_PATH}${RUN}/${SAMPLESHEET_SEQ}")
    if [ -z ${DESCRIPTION} ]
    then 
        change_status 0
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
        local FASTQ1=$(ls ${SEQ_PATH}${RUN}/Alignement_*/${RUN}*/Fastq/ | grep "^${SAMPLE_NAME}_.*_R1_.*") #TODO: check if path right for every sequencer
        local FASTQ2=$(ls ${SEQ_PATH}${RUN}/Alignement_*/${RUN}*/Fastq/ | grep "^${SAMPLE_NAME}_.*_R2_.*")
        echo "${SAMPLE_NAME},${SAMPLE_NAME},1,${FASTQ1},${FASTQ2}" >>"${SEQ_PATH}${RUN}/${SAMPLESHEET_SAREK}"
    done
}

launch_sarek () { 
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
        -w ${SEQ_PATH}${RUN}/${WORKDIR} \
        --input "${SEQ_PATH}${RUN}/${SAMPLESHEET_SAREK}" \
        --outdir "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}" \
        --step ${STEP} \
        --tools ${TOOLS} \
        --intervals ${BEDFILE} \
        --wes ${WES} \
        --save_output_as_bam ${SAVE_OUTPUT_AS_BAM} \
        --concatenate_vcfs ${CONCATENATE_VCFS} \
        --trim_fastq ${TRIM_FASTQ} \
        --split_fastq ${SPLIT_FASTQ} \
        ${PARAMS}"
    ${COMMAND} >"${SEQ_PATH}${RUNS_FILE}/${NF_STDOUT_FILE}"
}


############################ MAIN

for SEQ_PATH in ${MINISEQ_PATH} ${MISEQ_PATH} ${NEXTSEQ_PATH}
    do
    RUN_LIST=$(ls ${SEQ_PATH} | egrep '^[0-9]{6}_')

    # check runs file
    if [ ! -f "${SEQ_PATH}${RUNS_FILE}" ]
    then
        make_runs_file ${SEQ_PATH} ${RUN_LIST}
    fi

    for RUN in ${RUN_LIST}
        do
        LINE=$(grep "${RUN}" "${SEQ_PATH}${RUNS_FILE}")
        # add the new run to runs file
        if [ -z "${LINE}" ]
        then
              echo "${RUN} 1" >>"${SEQ_PATH}${RUNS_FILE}"
        fi
        get_bed_file 
        STATUS=${LINE##* }
        # check if the run is ready (status 1 and trigger file exists)
        if [ ${STATUS} = 1 && -f ${SEQ_PATH}${RUNS_FILE}/${TRIGGER_FILE} ]
        then
            create_samplesheet
            conda activate NF-core && launch_sarek
            # get error
            FAIL_LINES=$(grep FAILED "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}pipeline_info/execution_trace_*.txt" ) 
            # change trace name
            mv "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}pipeline_info/execution_trace_*.txt" "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}pipeline_info/_execution_trace_$(date +%F).txt"
            # if error is from custom dump software versions, continue, else, throw error
            if [ $(echo "${FAIL_LINES}" | wc -l) -eq 1 ] && [ -z "$(echo "${FAIL_LINES}" | grep DSV)" ]
            then
                # get temporary version file
                VERSIONS_FILE_PATH="${SEQ_PATH}${RUN}/${WORKDIR}$(echo $FAIL_LINES | awk '{print $2}')*/collated_versions.yml"
                TMP_VERSIONS_FILE_PATH="${SEQ_PATH}${RUN}/${WORKDIR}$(echo $FAIL_LINES | awk '{print $2}')*/tmp_collated_versions.yml"
                cp ${VERSIONS_FILE_PATH} ${TMP_VERSIONS_FILE_PATH}
                # modify the file to remove problematic java warnings
                FLAG=0
                while read VERSION_LINE
                    do
                    if [ ${FLAG} -eq 0 ]
                        then
                        if [ -z "$(echo "${VERSION_LINE}" | grep 'warning')" ]
                        then 
                            echo "${VERSION_LINE}" >>${TMP_VERSIONS_FILE_PATH}
                        else
                            TMP_VERSION_LINE="${VERSION_LINE%%:*}: "
                            FLAG=1
                        fi
                    else
                        TMP_VERSION_LINE+="${VERSION_LINE}"
                        echo "${TMP_VERSION_LINE}" >>${TMP_VERSIONS_FILE_PATH}
                        FLAG=0
                    fi
                done < ${VERSIONS_FILE_PATH}
                mv ${TMP_VERSIONS_FILE_PATH} ${VERSIONS_FILE_PATH}
                # launch sarek again with resume
                eval ${COMMAND} -resume >"${SEQ_PATH}${RUNS_FILE}/${NF_STDOUT_FILE}"
                FAIL_LINES=$(grep FAILED "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}pipeline_info/execution_trace_*.txt" )
                if [ -z "$(echo "${FAIL_LINES}")" ]
                then
                    # Pipeline completed ! change run status to done
                    change_status 2
                else
                    #TODO: throw error with run info
            else
                #TODO: throw error with run info
            conda deactivate
            #TODO: remove temporary files (deal with other logs, may need CWD)
        done
done

