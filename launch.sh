#!/usr/bin/env bash

des="ig:/srv/shiny-server/ampviewer/"

#rsync -rlzuvh \
rsync -ahvz \
   --no-perms --no-owner --no-group --ignore-times --omit-dir-times \
   --progress \
   --exclude=*.png \
   --exclude=make_tar.sh \
   --exclude=launch.sh \
   --exclude=*.Rproj \
   --exclude=.* \
   "./" "$des"
