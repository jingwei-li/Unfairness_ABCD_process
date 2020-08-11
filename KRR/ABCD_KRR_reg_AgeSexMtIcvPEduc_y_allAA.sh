#!/bin/sh
#
# Jingwei Li, 20200811

DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

########################
# setup for CIRC cluster
########################
curr_dir=$(pwd)
work_dir=${HOME}/cluster/

echo $curr_dir
echo $work_dir

if [ ! -d $work_dir ]; then
	mkdir -p $work_dir
fi

cd $work_dir


########################
# input setup
########################
proj_dir=/data/users/jingweil/storage/MyProject/fairAI/ABCD_race

