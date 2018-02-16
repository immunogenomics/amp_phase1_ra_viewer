#!/usr/bin/env bash

des="ig:/srv/shiny-server/ampviewer/"

rsync -rlzuvh \
   --progress \
   --exclude=*.png \
   --exclude=make_tar.sh \
   --exclude=launch.sh \
   --exclude=*.Rproj \
   --exclude=.* \
   "./" "$des"
