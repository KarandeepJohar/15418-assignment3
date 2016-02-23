#!/bin/bash
make clean
make APP="grade" GRAPH="/home/15-418/asst3_graphs/random_500m.graph" jobs
qsub jobs/${USER}_${APP}_${GRAPH}.job