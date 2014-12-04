#! /bin/sh
#$ -pe omp 2
#$ -l h_rt=24:00:00
#$ -N CONTOURS_OUT
#$ -V

export OMPI_MCA_btl=tcp,sm,self

/usr2/collab/ecowdery/EE509/Project/density.cluster.1.R

wait