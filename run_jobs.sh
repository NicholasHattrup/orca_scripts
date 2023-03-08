#!/bin/bash

# Set up variables
start=1
end=56
dist=1.2    
# Loop over range
for (( i=start; i<=end; i++ ))
do
    # Run ORCA job
    nohup $(which orca) job$i.inp > job$i.out &
    
    # Wait for job to finish
    wait %1
    ((dist=dist+.05))
    
    # Call Python script to generate next input file
    python write_job.py job$i.xyz job$((i+1)).inp $dist
done

