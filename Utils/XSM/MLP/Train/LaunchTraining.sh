#!/bin/bash

# Script to train MLPs
#@author BaL
#

echo "--------------------------"
echo "--- Run Training from MLP $1 to $2 ---"
echo "--------------------------"

#LigneDeDepart=$((0))
#NbreSimu=$((700))
LigneDeDepart=$1
LigneFinal=$2

echo LigneDeDepart $LigneDeDepart LigneFinal $LigneFinal exclu

for ((ligne=$LigneDeDepart; ligne<$LigneFinal; ligne=ligne+1 )) ;
do
		Train_XS $ligne
done
