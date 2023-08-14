#!/bin/csh

#--------------------------------------------------------------------# 
# $ ./test directoryName JobNum   (they can be anything)
#--------------------------------------------------------------------# 

echo directoryName=$1
echo jobNum=$2

#--------------------------------------------------------------------#  
# set environment
#--------------------------------------------------------------------#  
setenv EIC_LEVEL EIC2022a
source /cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/etc/eic_cshrc.csh
setenv HAS_PYTHIA8 true

#setenv PYTHIA8 "/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a"
#echo PYTHIA8: $PYTHIA8
#echo $LD_LIBRARY_PATH
#--------------------------------------------------------------------#
# Setting up paths
#--------------------------------------------------------------------#
set DELPHESDIR="/gpfs/mnt/gpfs02/eic/pingwong/delphes_EIC2"
set OUTPUTDIR="/gpfs/mnt/gpfs02/eic/pingwong/delphes_EIC2/$1"
#--------------------------------------------------------------------#
# Run DELPHES with PYTHIA
#--------------------------------------------------------------------#
cd $DELPHESDIR
echo location: `pwd`

#$DELPHESDIR/DelphesPythia8 $DELPHESDIR/2ndDET_cards/delphes_card_ePIC.tcl $DELPHESDIR/2ndDET_pythia8cards/CC_DIS.cmnd $2.root

if ( -e $2.root ) then
    rm $2.root
endif

./DelphesPythia8 2ndDET_cards/delphes_card_ePIC.tcl 2ndDET_pythia8cards/CC_DIS.cmnd $2.root

if (! -e "${OUTPUTDIR}" ) then
   mkdir $OUTPUTDIR
endif

mv $2.root $OUTPUTDIR
