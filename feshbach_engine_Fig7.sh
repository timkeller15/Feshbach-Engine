#!/bin/bash

# Find ground state for ramp starting point, can be omitted if already done for the data in Figs. 2 and 4. 
./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --Tf=1.4 --imaginarytime --atoms=10000 --gi=1.0 --rampmode=2 --input=groundstate --writewfc

cp data/feshbach_engine_ramp_durations_Fig7.txt data/ramp_durations.txt

for i in {1..9}
do
	# Find ground state for ramp target point
	./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --Tf=1.4 --imaginarytime --atoms=10000 --gi=$(bc <<<"scale=1; $i/10") --rampmode=2 --input=groundstate --writewfc
	
	# Perform the same shortcut ramp for a range of durations
	./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --atoms=10000 --gi=1.0 --gf=$(bc <<<"scale=1; $i/10") --rampmode=1 --input=groundstate --target=groundstate --batchmode --filename=stability
done

rm data/ramp_durations.txt