tname=$1
mkdir $tname
echo $tname
cp *.py $tname/
cp *.cpp $tname/
cp *.hpp $tname/
cd $tname/
g++ -o mh.o  mh.cpp  util.cpp `gsl-config --cflags --libs`
cp ../$tname.ssm.txt ssm_data.txt
cp ../$tname.cnv.txt cnv_data.txt
python evolve.py

