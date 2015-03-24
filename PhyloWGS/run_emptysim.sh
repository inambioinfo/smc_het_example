tname=$1
mkdir $tname.output
echo $tname
cp *.py $tname.output/
cp *.cpp $tname.output/
cp *.hpp $tname.output/
cd $tname.output/
g++ -o mh.o  mh.cpp  util.cpp `gsl-config --cflags --libs`
cp ../$tname ssm_data.txt
touch cnv_data.txt
python evolve.py

