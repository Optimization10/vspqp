#!/bin/bash

# Update the code from local computer to remote computer
rsync -rtv -e "ssh" --delete \
      --exclude 'output' --exclude 'outputs_archive' --exclude '/arch/' \
      --exclude 'random_seed.txt' --exclude 'Build' --exclude '.*' \
      --exclude 'hager6-ssh/' --exclude '*.zst' \
      /home/ovi/vsp/vspqp/ \
      onechita@ls-hager6.math.ufl.edu:/home/onechita/vspqp/
