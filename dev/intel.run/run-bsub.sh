#!/bin/bash
#
#BSUB -n 1
#BSUB -q knlbq
#BSUB -R "{select[ekl] span[ptile=1]}"
#BSUB -J QE
#BSUB -o out-%J.txt
#BSUB -e out-%J.err
#BSUB -C 0
#
ROOT=${HOME}/espresso-5.4.0

WORKLOAD=${HOME}/espresso-data/AUSURF112/ausurf.in
#NUMACTL="numactl --preferred=1"

MAXNCORES=64
NRANKS=32
MAXNT=4
NT=1

SHIFT=$((((2*MAXNCORES+NRANKS-1)/NRANKS)/2))
XRANKS=$((MAXNCORES/SHIFT))
if [ "1" = "$((XRANKS<=NRANKS))" ]; then NRANKS=${XRANKS}; fi
NTHREADS=$((SHIFT*NT))
NUMNODES=$(cat ${PBS_NODEFILE} | wc -l)

#ARGS="-npool 2"
NPOOL=$(echo ${ARGS} | sed -n -e "s/.*-npool\s\s*\([0-9][0-9]*\).*/\1/p")
if [ "" = "${NPOOL}" ]; then NPOOL=1; fi
if [ "" = "$(echo ${ARGS} | grep '\-ndiag')" ]; then
  ARGS+=" -ndiag $((NUMNODES*NRANKS/NPOOL))"
fi
if [ "" = "$(echo ${ARGS} | grep '\-ntg')" ]; then
  NTG=$((NUMNODES*NRANKS/(NPOOL*2)))
  if [ "0" = "${NTG}" ]; then NTG=${NRANKS}; fi
  ARGS+=" -ntg ${NTG}"
fi

if [ "${MAXNT}" = "${NT}" ]; then
  AFFINITY=compact
else
  AFFINITY=scatter
fi

#export PSM2_MQ_RNDV_HFI_WINDOW=4194304
#export PSM2_MQ_EAGER_SDMA_SZ=65536
#export PSM2_MQ_RNDV_HFI_THRESH=200000
#export PSM2_IDENTIFY=1

#export I_MPI_HYDRA_PMI_CONNECT=alltoall
#export I_MPI_USE_DYNAMIC_CONNECTIONS=0
#export HFI_NO_CPUAFFINITY=1
#export IPATH_NO_CPUAFFINITY=1
#export I_MPI_SCALABLE_OPTIMIZATION=0

#export I_MPI_PIN_DOMAIN=node
#export I_MPI_PIN_MODE=lib

# I_MPI_FABRICS(I_MPI_DAPL_PROVIDER):
#   OPA: dapl(ofa-v2-hib0)
#        shm:tmi(psm2)
#export I_MPI_FABRICS=shm:dapl
export I_MPI_FABRICS=shm:tmi
export I_MPI_TMI_PROVIDER=psm2
RUN="mpirun -bootstrap ssh -genvall \
  -np $((NRANKS*NUMNODES)) -perhost ${NRANKS} -genv I_MPI_FALLBACK=0 -genv xI_MPI_DEBUG=4 \
  -genv I_MPI_PIN_DOMAIN=auto -genv I_MPI_PIN_ORDER=scatter \
  -genv KMP_AFFINITY=${AFFINITY},granularity=fine,1 \
  -genv OMP_NUM_THREADS=${NTHREADS} \
  ${NUMACTL} \
${ROOT}/bin/pw.x -i ${WORKLOAD} ${ARGS}"

source ${ROOT}/enable
#export LD_LIBRARY_PATH=${ROOT}/lib:${LD_LIBRARY_PATH}

mpdboot -n $(cat ${PBS_NODEFILE} | sort | uniq | wc -l) \
  -r ssh -f ${PBS_NODEFILE}
mpdtrace
cd $(dirname ${WORKLOAD})
echo "${RUN}"
eval ${RUN}
mpdallexit

