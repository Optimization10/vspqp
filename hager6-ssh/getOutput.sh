#!/bin/bash

# Copy output
rsync -rtv -e ssh \
      onechita@ls-hager6.math.ufl.edu:/home/onechita/codes/VSPQP/output/ \
      /home/ovi/Desktop/hager6-vspqp-output/