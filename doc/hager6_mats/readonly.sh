#!/bin/bash

find ./doc/hager6_mats/ -type f ! -path "./doc/hager6_mats/readonly.sh" -exec chmod a-w {} \;