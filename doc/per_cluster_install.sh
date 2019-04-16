#!/usr/bin/env bash

module purge
module load intel/2018a
module load Libint/2.4.2-intel-2018a
module load CMake/3.10.3-GCCcore-6.4.0
module load Eigen/3.3.4
module load Boost/1.66.0-intel-2018a

case ${VSC_INSTITUTE_CLUSTER} in
    "golett" )
	PPN=1
	MEM=3gb
	;;
    "swalot" )
	PPN=1
	MEM=3gb
	;;
    "phanpy" )
	PPN=1
	MEM=3gb
	;;
    # "delcatty" )
    # 	export LIBINT_ROOT=/apps/gent/CO7/sandybridge/software/Libint/2.4.2-intel-2018a
    # 	PPN=8
    # 	MEM=30gb
    # 	;;
    # "victini" )
    #     export LIBINT_ROOT=/apps/gent/CO7/skylake-ib/software/Libint/2.4.2-intel-2018a
    # 	PPN=18
    # 	MEM=30gb
    # 	;;
    * )
	echo "ERROR: Only the golett, swalot and phanpy clusters are supported."
	exit 1
	;;
esac

export EIGEN3_ROOT=$EBROOTEIGEN/include
export gqcp_DIR=${VSC_DATA}/apps/${VSC_INSTITUTE_CLUSTER}/gqcg/gqcp/
# 1. install
cd $HOME/drivers
rm -rf $VSC_INSTITUTE_CLUSTER
mkdir $VSC_INSTITUTE_CLUSTER
cd $VSC_INSTITUTE_CLUSTER
git clone https://github.com/tmhuysen/tmhuysen-gqcp-drivers.git
(cd tmhuysen-gqcp-drivers && rm -rf build && mkdir build && cd build && cmake .. -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc && make VERBOSE=1 -j ${PPN} && make test ARGS=-j${PPN})