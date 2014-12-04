#! /bin/sh
#$ -pe omp 2
#$ -l h_rt=24:00:00
#$ -N CONTOURS_OUT
#$ -V

export OMPI_MCA_btl=tcp,sm,self

/home/ecowdery/EE509/Project/contours.cluster.R

wait