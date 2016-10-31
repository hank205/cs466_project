mkdir outcomes
cd outcomes

mkdir default
mkdir a_icpc1
mkdir a_icpc1.5
mkdir b_ml6
mkdir b_ml7
mkdir c_sc5
mkdir c_sc20
cd ..

python step2.py default
mv 'predictedmotif.txt' outcomes/default/'predictedmotif.txt'
mv 'predictedsites.txt' outcomes/default/'predictedsites.txt'

python step2.py a_icpc1
mv 'predictedmotif.txt' outcomes/a_icpc1/'predictedmotif.txt'
mv 'predictedsites.txt' outcomes/a_icpc1/'predictedsites.txt'

python step2.py a_icpc1.5
mv 'predictedmotif.txt' outcomes/a_icpc1.5/'predictedmotif.txt'
mv 'predictedsites.txt' outcomes/a_icpc1.5/'predictedsites.txt'

python step2.py b_ml6
mv 'predictedmotif.txt' outcomes/b_ml6/'predictedmotif.txt'
mv 'predictedsites.txt' outcomes/b_ml6/'predictedsites.txt'

python step2.py b_ml7
mv 'predictedmotif.txt' outcomes/b_ml7/'predictedmotif.txt'
mv 'predictedsites.txt' outcomes/b_ml7/'predictedsites.txt'

python step2.py c_sc5
mv 'predictedmotif.txt' outcomes/c_sc5/'predictedmotif.txt'
mv 'predictedsites.txt' outcomes/c_sc5/'predictedsites.txt'

python step2.py c_sc20
mv 'predictedmotif.txt' outcomes/c_sc20/'predictedmotif.txt'
mv 'predictedsites.txt' outcomes/c_sc20/'predictedsites.txt'