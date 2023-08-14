#!/bin/bash

condorLog="/gpfs/mnt/gpfs02/eic/pingwong/delphes_EIC2/condorLog"
echo create directory $condorLog

if [ ! -e "$condorLog" ]
then
    mkdir $condorLog
fi

cd $condorLog

if [ ! -e job ]
then
    mkdir job
fi

if [ ! -e "out" ]
then
    mkdir out
fi

if [ ! -e "log" ]
then
    mkdir log
fi

if [ ! -e "err" ]
then
    mkdir err
fi
