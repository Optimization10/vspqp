#!/bin/bash

# Copy output
rsync -rtv -e ssh \
      onechita@ls-hager6.math.ufl.edu:/home/onechita/VSPQP/output/ \
      /home/ovi/Desktop/VSPQP/output/

# Copy outputs_archive                                      
rsync -rtv -e ssh \
       onechita@ls-hager6.math.ufl.edu:/home/onechita/VSPQP/outputs_archive/ \
       /home/ovi/Desktop/VSPQP/outputs_archive/

# Copy MultipleRunsResults.txt                                      
rsync --append -rtv -e ssh \
       onechita@ls-hager6.math.ufl.edu:/home/onechita/VSPQP/doc/MultipleRunsResults.txt \
       /home/ovi/Desktop/VSPQP/doc/MultipleRunsResults.txt