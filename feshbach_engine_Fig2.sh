#!/bin/bash

# Finding ground states of the four cycle end points
./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --Tf=1.4 --imaginarytime --atoms=10000 --gi=1.0 --rampmode=2 --input=groundstate --writewfc
./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --Tf=1.4 --imaginarytime --atoms=10000 --gi=0.8 --rampmode=2 --input=groundstate --writewfc
./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --Tf=1.4 --imaginarytime --atoms=8000 --gi=1.0 --rampmode=2 --input=groundstate --writewfc
./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --Tf=1.4 --imaginarytime --atoms=8000 --gi=0.8 --rampmode=2 --input=groundstate --writewfc

cp data/feshbach_engine_ramp_durations_Fig2.txt data/ramp_durations.txt

# Perform compression strokes using both the shortcut to adiabaticity and the time-rescaled adiabatic reference ramp for a range of durations
./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --atoms=10000 --gi=1.0 --gf=0.8 --rampmode=1 --input=groundstate --target=groundstate --batchmode --filename=compression_stroke
./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --atoms=10000 --gi=1.0 --gf=0.8 --rampmode=0 --input=groundstate --target=groundstate --batchmode --filename=compression_stroke

# Perform expansion strokes using both the shortcut to adiabaticity and the time-rescaled adiabatic reference ramp for a range of durations
./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --atoms=8000 --gi=0.8 --gf=1.0 --rampmode=1 --input=groundstate --target=groundstate --batchmode --filename=expansion_stroke
./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --atoms=8000 --gi=0.8 --gf=1.0 --rampmode=0 --input=groundstate --target=groundstate --batchmode --filename=expansion_stroke

rm data/ramp_durations.txt