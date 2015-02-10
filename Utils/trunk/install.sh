echo "********************************"
echo "******INSTALLING CLASS V4 ******"
echo "********************************"

mkdir -p lib
cd source/src
make clean; make -j
make install
echo "********************************"
cd $CLASS_PATH/gui
make clean; make -j
echo "********************************"
cd $CLASS_PATH/DATA_BASES/DECAY/ALL/
sed -i -e "s%/PATHTOBASE%`pwd`%" Decay.idx
echo "Decay base Configuration done!"
echo "********************************"
echo "**** INSTALLATION COMPLETED ****"
echo "********************************"