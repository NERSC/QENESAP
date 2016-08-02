#!/bin/sh
#
ROOT=${PWD}
#ROOT=$(cd ..; pwd -P)
MAXNCORES=64
MAXNT=4
NT=1

# consider the following in your .bashrc
# ulimit -s unlimited
# ulimit -c0

#PREFX="perf stat -e tlb:tlb_flush,irq_vectors:call_function_entry,syscalls:sys_enter_munmap,syscalls:sys_enter_madvise,syscalls:sys_enter_brk "
#VTUNE="-gtool 'amplxe-cl -r vtune -data-limit 0 -collect advanced-hotspots -knob collection-detail=stack-sampling:4=exclusive'"
#ADVXE="-gtool 'amplxe-cl -r vtune -data-limit 0 -collect advanced-hotspots -knob collection-detail=stack-sampling:4=exclusive'"
#NUMACTL="numactl --preferred=1"

if [ -f ${ROOT}/enable ]; then
  source ${ROOT}/enable
fi

HOSTS=$(${ROOT}/mynodes.sh 2> /dev/null | tr -s "\n" "," | tr -s " " "," | sed -e "s/^\(..*[^,]\),*$/\1/" | cut -d, -f1-${NUMNODES})
if [ "" = "${HOSTS}" ]; then HOSTS=localhost; fi
HOST=$(echo ${HOSTS} | cut -d, -f1)

if [ "" != "$1" ] && [ -f $1 ]; then
  WORKLOAD=$1
  shift
else
  WORKLOAD=${HOME}/espresso-data/AUSURF112/ausurf.in
fi
WORKLOAD=$(cd $(dirname ${WORKLOAD}); pwd -P)/$(basename ${WORKLOAD})

if [ "" != "$1" ]; then
  NUMNODES=$1
  shift
else
  NUMNODES=1
fi

if [ "" != "$1" ]; then
  NRANKS=$1
  shift
else
  NRANKS=16
fi

SHIFT=$((((2*MAXNCORES+NRANKS-1)/NRANKS)/2))
XRANKS=$((MAXNCORES/SHIFT))
if [ "1" = "$((XRANKS<=NRANKS))" ]; then NRANKS=${XRANKS}; fi
NTHREADS=$((SHIFT*NT))

ARGS=$*
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
#export PSM2_MQ_RNDV_HFI_THRESH=200000
#export PSM2_MQ_EAGER_SDMA_SZ=65536
#export PSM2_IDENTIFY=1

export I_MPI_FALLBACK=0
export I_MPI_SHM_LMT=shm
export I_MPI_HYDRA_PMI_CONNECT=alltoall
#export I_MPI_USE_DYNAMIC_CONNECTIONS=0
#export I_MPI_SCALABLE_OPTIMIZATION=0

#export IPATH_NO_CPUAFFINITY=1
#export HFI_NO_CPUAFFINITY=1

# I_MPI_FABRICS(I_MPI_DAPL_PROVIDER):
#   OPA: dapl(ofa-v2-hib0)
#        shm:tmi(psm2)
#export I_MPI_FABRICS=shm:dapl
#export I_MPI_DAPL_PROVIDER=ofa-v2-mlx5_0-1u
#export I_MPI_FABRICS=shm:tmi
#export I_MPI_TMI_PROVIDER=psm2
RUN="${PREFX}mpirun -bootstrap ssh -genvall -host ${HOSTS} \
  -np $((NRANKS*NUMNODES)) -perhost ${NRANKS} -genv xI_MPI_DEBUG=4 \
  -genv I_MPI_PIN_DOMAIN=auto -genv I_MPI_PIN_ORDER=scatter \
  -genv KMP_AFFINITY=${AFFINITY},granularity=fine,1 \
  -genv OMP_NUM_THREADS=${NTHREADS} \
  ${VTUNE} ${ADVXE} ${NUMACTL} \
${ROOT}/bin/pw.x -i ${WORKLOAD} ${ARGS}"

cd $(dirname ${WORKLOAD})
ssh ${HOST} "env | tr '\n' ' '"
echo; echo
ssh ${HOST} "numactl -H"
echo
echo "${RUN}"
eval ${RUN}
