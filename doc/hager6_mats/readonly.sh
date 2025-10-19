#!/bin/bash

# Make all files read-only except the script
find ./doc/hager6_mats/ -type f ! -name "readonly.sh" -exec chmod a-w {} \;

# Make all folders read-only
find ./doc/hager6_mats/ -type d -exec chmod a-w {} \;
