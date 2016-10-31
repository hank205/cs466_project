mkdir data_sets
cd data_sets

mkdir default
cd default
python ../../step1.py 2 8 500 10
cd ..

mkdir a_icpc1
cd a_icpc1
python ../../step1.py 1 8 500 10
cd ..

mkdir a_icpc1.5
cd a_icpc1.5
python ../../step1.py 1.5 8 500 10
cd ..

mkdir b_ml6
cd b_ml6
python ../../step1.py 2 6 500 10
cd ..

mkdir b_ml7
cd b_ml7
python ../../step1.py 2 7 500 10
cd ..

mkdir c_sc5
cd c_sc5
python ../../step1.py 2 8 500 5
cd ..

mkdir c_sc20
cd c_sc20
python ../../step1.py 2 8 500 20
cd ..