#!/bin/bash

# Update the code from local computer to remote computer
rsync -rtv -e "ssh" --delete \
      --exclude 'data' --exclude 'output' --exclude 'outputs_archive' \
      --exclude 'random_seed.txt' --exclude 'build' --exclude '.*' \
      --exclude 'hager6/' --exclude '*.zst' \
      /home/ovi/Desktop/VSP/VSPQP/ \
      onechita@ls-hager6.math.ufl.edu:/home/onechita/codes/VSPQP/