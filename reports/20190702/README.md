Output_meme and mast, comes from running meme on uniprot MG ligand, with uniref50, leading to ~9.5k seqs. Settings are:

/home/yincp/meme-5.0.5/bin/meme -text -w 30 -p 80 -nmotifs 50 /home/yincp/mg_50.fasta > /home/yincp/test.txt


Complete script is:

[yincp@santa ~]$ cat ibqueue_toMelvin.sh
#PBS -q batch
#PBS -l nodes=1:ppn=20
#PBS -l walltime=100:0:0

if [ ! -f $HOME/.mpd.conf ] ; then
  echo "FATAL! No $HOME/.mpd.conf"
  exit 1
fi

#-------------------- Defining the nodes and cpus ----------#
NODES=1
TOTAL_CPUS=20

source /cluster/apps/x86_64/packages/modules-tcl/init/bash
module load mpich2

#---- Parse $PBS_NODEFILE, and generate $MACHINE_FILE ---------#
MACHINE_FILE=$PBS_O_WORKDIR/machine_file.$$
sort < $PBS_NODEFILE | uniq > $MACHINE_FILE

#--- Startup MPD, do a quick test before we proceed ---------#

$MPI_DIR/bin/mpdboot -n $NODES -f $MACHINE_FILE
if [ $? -ne 0 ] ; then
  exit 1
fi
RESULT=`$MPI_DIR/bin/mpdtrace | wc -l`
if [ "$RESULT" != "$NODES" ] ; then
  $MPI_DIR/bin/mpdallexit
  exit 1
fi

cd $PBS_O_WORKDIR

/home/yincp/meme-5.0.5/bin/meme -text -w 30 -p 80 -nmotifs 50 /home/yincp/mg_50.fasta > /home/yincp/test.txt

Interpretation:
What we want (DXDXDG) is pretty far down the list and has rather few matches. Can consider matching without a fixed width, we might be able to get better matches that way. Raise to Igor this point. 
