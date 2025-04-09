abs_path () {
        case $1 in
                /*) echo ${1%/}
                ;;
                *)  echo ${PWD}/${1%/}
        esac
}

if [ -v SLURM_NPROCS ]; then
        CPUS=$SLURM_NPROCS

        # May be needed to make use of all threads
        export OMP_NUM_THREADS=$CPUS
else
        CPUS=1
fi

if [ -v SLURM_JOB_ID ]; then
        WORK=${PWD}/work/slurm/slurm_${SLURM_JOB_ID}
else
        WORK=${PWD}/work/$(basename $0)/$(date +%F_%H-%M-%S)
fi

if [ ! -z "$(which pigz)" ]; then
        GZIP=pigz
else
        GZIP=gzip
fi
