#!/bin/bash

# Copy output
mkdir -p /home/ovi/Desktop/hager6-vspqp-output/new
rsync -rtv --ignore-times --delete -e ssh \
      onechita@ls-hager6.math.ufl.edu:/home/onechita/codes/VSPQP/output/ \
      /home/ovi/Desktop/hager6-vspqp-output/new/