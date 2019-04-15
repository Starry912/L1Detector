# Data Simulator for LINE-1 insertion

## require:
## python version '2.7.12 (default, Nov 12 2018, 14:36:49) [GCC 5.4.0 20160609]'
## python package: numpy

include:
1. L1create.py
Simulate1:create L1 sequence
2. L1insertion.py
Simulate2:insert L1seq into genome
3. NGS_sim.py
Simulate3:NGS sequence


## USING
`python NGS_sim.py coverage insert_number`
example `python ./simulation/NGS_sim.py 0.8 2000`


