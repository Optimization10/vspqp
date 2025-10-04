#!/bin/bash

# Copy output
rsync -rtv -e ssh \
      onechita@ls-hager6.math.ufl.edu:/home/onechita/VSPQP/output/ \
      /home/ovi/Desktop/VSP/VSPQP/output/