# -*- coding: utf-8 -*-

import os
from remoteMachine import RemoteMachine

machines = {}


remote = {}
users=[u"ramet",u"lacoste"]
machinename=u"devel09.plafrim.cluster-%s"
for user in users:
    remote[machinename % user] = RemoteMachine(u"mygale.bordeaux.inria.fr",
                                               user,
                                               [u"ssh"],
                                               [u"scp"],
                                               # not in tmp as it is not shared 
                                               # between machines
                                               u'/lustre/%s/results' % user,
                                               u'qsub',
                                               qstat = u'qstat -f | grep __JOBID__ | wc -l',
                                               sleeptime = 5)


script = u"""
#PBS -N __JOBID__
#PBS -j oe
#PBS -l walltime=00:05:00 
#PBS -l nodes=__PROCNBR__:ppn=__THRDNBR__
### Positionne tous les noeud sur un meme switch, pour plus de 16 le job
### restera bloque indefiniment
### P_B_S -W x=\"NODESET:ONEOF:FEATURE:ib001:ib002:ib005:ib006\"

source /etc/profile.d/modules.sh

# Chargement des modules :
module purge
source /opt/intel/mkl/10.2.7.041/tools/environment/mklvarsem64t.sh
# Positionnement dans le bon repertoire qui va bien :
# cat $PBS_NODEFILE
# cat $PBS_NODEFILE | sort | uniq > hostfile
cd __DIR__
echo '/opt/cluster/mpi/openmpi/latest/gcc/latest/bin/mpiexec -n __PROCNBR__ __CMD__' > __LOG__ 2>&1
/opt/cluster/mpi/openmpi/latest/gcc/latest/bin/mpiexec -npernode 1 -n __PROCNBR__ __CMD__ >> __LOG__ 2>&1

iter=0
while [ ! -s __LOG__  -a $iter -le 10 ]
  do
  sleep 1
  iter=`expr $iter + 1`
  done

# On attend encore 2 secondes de plus au cas ou.
sleep 2
if [ ! -s __LOG__ ]; then echo "fichier de log absent..."; fi
"""

name   = u"plafrim-gnu-openmpi-mkl-%s"
libmkl = u'-L/opt/intel/mkl/10.2.7.041/lib/em64t/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core'
libgfortran = u"-lgfortran"

for user in users :
    machines[name % user] = remote[machinename % user].createConfig(name % user, 
                                                                    machinename % user, 
                                                                    32, 
                                                                    script,
                                                                    scotch_home = u"/opt/cluster/plafrim-dev/scotch-5.1.9/gcc-4.3.2/openmpi-1.4.2/",
                                                                    hwloc_home  = u"/opt/cluster/lib/hwloc/latest",
                                                                    cc          = u"/usr/bin/gcc", 
                                                                    fc          = u"/usr/bin/gfortran", 
                                                                    mpicc       = u'LD_LIBRARY_PATH=/opt/cluster/lib/mpc/0.8.2/lib:/opt/cluster/lib/mpfr/3.0.0/lib:/opt/cluster/lib/libelf/0.8.13/lib:/opt/cluster/lib/gmp/5.0.1/lib:$LD_LIBRARY_PATH /opt/cluster/mpi/openmpi/latest/gcc/latest/bin/mpicc', 
                                                                    mpifc       = u'LD_LIBRARY_PATH=/opt/cluster/lib/mpc/0.8.2/lib:/opt/cluster/lib/mpfr/3.0.0/lib:/opt/cluster/lib/libelf/0.8.13/lib:/opt/cluster/lib/gmp/5.0.1/lib:$LD_LIBRARY_PATH /opt/cluster/mpi/openmpi/latest/gcc/latest/bin/mpif90',
                                                                    fcflags     = u'', 
                                                                    ccflags     = u'', 
                                                                    ccopt       = u'-O3', 
                                                                    ccdeb       = u'-g3',
                                                                    libblas     = libmkl, 
                                                                    libfortran  = libgfortran,
                                                                    libmaths    = u'-lm', 
                                                                    librt       = u'-lrt',
                                                                    ar          = u'ar', 
                                                                    arflags     = u'-ruv', 
                                                                    make        = u'make', 
                                                                    make_j      = u'make -j 8')
