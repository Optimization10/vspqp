#!/bin/bash

# Give write access to all files
find ./doc/hager6_mats/ -type f -exec chmod a+w {} \;

# Give write access to all folders
find ./doc/hager6_mats/ -type d -exec chmod a+w {} \;