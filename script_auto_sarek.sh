#!/bin/bash

set -e

# load variables
source auto_sarek.conf  #TODO: modify path
HEADER="[HEADER]
A list of all the runs in this directory and their status regarding the sarek workflow
    0: to ignore
    1: to be treated
    2: done
    
[FILES]"

############################ FUNCTIONS

# create sarek_runs.txt file containing run names and status
make_runs_file () {
    # argument 1 runs ouput directory
    # argument 2 run list
    echo "${HEADER}" >"$1${RUNS_FILE}"
    for RUN in $2
        do
        echo "${RUN} 0" >>"$1${RUNS_FILE}"
        done
}

change_status () {
    # argument : new status
    sed -i "s/${RUN_NAME} .*/${RUN_NAME} $1/" "${SEQ_PATH}${RUNS_FILE}"
}

# get bed file from samplesheet. If does not exist, modify status to 0 (ignore)
get_bed_file () {
    if [ -f "${SEQ_PATH}${RUN}/${SAMPLESHEET_SEQ}" ]
    then 
        local DESCRIPTION=$(grep "Description.*\.bed" "${SEQ_PATH}${RUN}/${SAMPLESHEET_SEQ}")
        if [ -z "${DESCRIPTION}" ]
        then 
            change_status 0
            LINE=$(grep "${RUN}" "${SEQ_PATH}${RUNS_FILE}")
        else
            local tmp=${DESCRIPTION#*,}
            BEDFILE=${tmp%#*}
            #TODO: check if bed exists in valid beds list
        fi
    else
        change_status 0
    fi
}

create_samplesheet () {
    # header
    echo "patient,sample,lane,fastq_1,fastq_2" >"${SEQ_PATH}${RUN}/${SAMPLESHEET_SAREK}"
    # fastq path in the run directory
    case "${SEQ_PATH}" in
        "${MINISEQ}") local FASTQ_PATH="Alignment_*/${RUN}*/Fastq/" ;;
        "${MISEQ}") local FASTQ_PATH="Data/Intensities/BaseCalls/" ;;
        "${NEXTSEQ}") local FASTQ_PATH="FastQs/" ;;
    esac
    # loop on Samplesheet.csv on Data after header
    sed -n "/\[Data\]/ {n; :a; n; p; ba}" "${SEQ_PATH}${RUN}/${SAMPLESHEET_SEQ}" | \
    while read -r SAMPLESHEET_LINE
        do
        local SAMPLE_NAME=${SAMPLESHEET_LINE%%,*}
        local FASTQ1=$(ls ${SEQ_PATH}/${RUN}/${FASTQ_PATH} | grep "^${SAMPLE_NAME}_.*_R1_.*")
        local FASTQ2=$(ls ${SEQ_PATH}/${RUN}/${FASTQ_PATH}/ | grep "^${SAMPLE_NAME}_.*_R2_.*")
        echo "${SAMPLE_NAME},${SAMPLE_NAME},1,${FASTQ1},${FASTQ2}" >>"${SEQ_PATH}${RUN}/${SAMPLESHEET_SAREK}"
    done
}

launch_sarek () { 
    # get genome version
    VGENOME=$(echo "${BEDFILE}" | grep -oE 'hg[0-9]{2}')
    if [ "${VGENOME}" = "hg38" ] || [[ "${BEDLIST38}" =~ ${BEDFILE} ]]
    then 
        # string with genome-specific parameters
        PARAMS="--genome ${GENOME_38} \
                --igenomes_base ${IGENOMES_BASE_38} \
                --fasta ${FASTA_38} \
                --fasta_fai ${FASTA_FAI_38}"
    fi
    if [ "${VGENOME}" = "hg19" ] || [[ "${BEDLIST19}" =~ ${BEDFILE} ]]
    then 
        # string with genome-specific parameters
        PARAMS="--genome ${GENOME_19} \
                --igenomes_base ${IGENOMES_BASE_19} \
                --fasta ${FASTA_19} \
                --fasta_fai ${FASTA_FAI_19}"
    fi
    #TODO: if want to add an else to notify unknown bed, do cases ? or second in first's else
    
    COMMAND="nextflow run ${SAREK_PATH} -profile singularity \
        -c ${CONFIG_PATH} \
        -w "${SEQ_PATH}${RUN}/${WORKDIR}" \
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
    ${COMMAND} &> "${SEQ_PATH}${RUN}/${NF_STDOUT_FILE}"
}


############################ MAIN

for SEQ_PATH in ${MINISEQ_PATH} ${MISEQ_PATH} ${NEXTSEQ_PATH}
    do
    RUN_LIST=$(ls "${SEQ_PATH}" | grep -E '^[0-9]{6}_')

    # check runs file
    if [ ! -f "${SEQ_PATH}${RUNS_FILE}" ]
    then
        make_runs_file "${SEQ_PATH}" "${RUN_LIST}"
    fi

    for RUN in ${RUN_LIST}
        do
        set +e
        LINE=$(grep "${RUN}" "${SEQ_PATH}${RUNS_FILE}")
        set -e
        # add the new run to runs file
        if [ -z "${LINE}" ]
        then
              echo "${RUN} 1" >>"${SEQ_PATH}${RUNS_FILE}"
              LINE=$(grep "${RUN}" "${SEQ_PATH}${RUNS_FILE}")
        fi
        get_bed_file
        STATUS=${LINE##* }
        # check if the run is ready (status 1 and trigger file exists)
        if [[ ${STATUS} = 1 && -f ${SEQ_PATH}${RUN}/${TRIGGER_FILE} ]]
        then
            create_samplesheet
            conda activate NF-core && launch_sarek
            # get error
            FAIL_LINES=$(grep FAILED "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}pipeline_info/execution_trace_*.txt" )
            # change trace name
            mv "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}pipeline_info/execution_trace_*.txt" "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}pipeline_info/_execution_trace_$(date +%F).txt"
            # if error is from custom dump software versions, continue, else, throw error
            if [[ $(echo "${FAIL_LINES}" | wc -l) -eq 1 ]] && [ -n "$(echo "${FAIL_LINES}" | grep DSV)" ]
            then
                # get temporary version file
                VERSIONS_FILE=$(realpath "${SEQ_PATH}${RUN}/${WORKDIR}/tmp/*/*/collated_versions.yml")
                VERSIONS_FILE_TMP=tmp_collated_versions
                # modify the file to remove problematic java warnings
                FLAG=0
                while read -r VERSION_LINE
                    do
                    if [ ${FLAG} -eq 0 ]
                    then
                        if [ -z "$(echo "${VERSION_LINE}" | grep 'warning')" ]
                        then
                            if [ -z "$(echo "${VERSION_LINE}" | grep NFCORE_SAREK)" ]
                            then
                                echo "    ${VERSION_LINE}" >>${VERSIONS_FILE_TMP}
                            else    
                                echo "${VERSION_LINE}" >>${VERSIONS_FILE_TMP}
                            fi
                        else
                            TMP_VERSION_LINE="${VERSION_LINE%%:*}: "
                            FLAG=1
                        fi
                    else
                        TMP_VERSION_LINE+="${VERSION_LINE}"
                        echo "    ${TMP_VERSION_LINE}" >>${VERSIONS_FILE_TMP}
                        FLAG=0
                    fi
                done < "${VERSIONS_FILE}"
                # overwrite version file
                mv ${VERSIONS_FILE_TMP} "${VERSIONS_FILE}"
                # launch sarek again with resume
                eval "${COMMAND} -resume &> ${SEQ_PATH}${RUN}/${NF_STDOUT_FILE}"
                FAIL_LINES=$(grep FAILED "${SEQ_PATH}${RUN}/${OUTDIR_SAREK}pipeline_info/execution_trace_*.txt" )
                if [ -z "${FAIL_LINES}" ]
                then
                    # Pipeline completed ! change run status to done
                    change_status 2
                #else
                    #TODO: throw error with run info
                fi
            #else
                #TODO: throw error with run info
            fi
            conda deactivate
            #TODO: remove temporary files (deal with other logs, may need CWD)
        fi
    done
done

