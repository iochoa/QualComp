#!/bin/bash

# Separate file in two, for QVs and no QVs
awk 'NR%4==0' $1 > $2
awk 'NR%4!=0' $1 > $3
