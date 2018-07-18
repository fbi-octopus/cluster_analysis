#!/bin/sh
#BSUB -n 1
#BSUB -o "/<path-to-data>/sample_data_roi/s07a02r21416-19926_756x3674_26green/output-%J.log"
#BSUB -e "/<path-to-data>/sample_data_roi/s07a02r21416-19926_756x3674_26green/error-%J.log"
#BSUB -W 72:00
#BSUB -M 64000000
#BSUB -J <some-job-name>
cd <path-to-the-R-scripts>
R -f run.R --args --target=/<path-to-data>/sample_data_roi/s07a02r21416-19926_756x3674_26green --config=/<path-to-data>/sample_data_roi/s07a02r21416-19926_756x3674_26green/config.txt
cd -
