#!/usr/bin/env bash

src="/Users/slowikow/Dropbox/work/immunogenomics/ampviewer/"
des="ig:/srv/shiny-server/ampviewer/"

rsync -rlzuvh \
   --exclude=*.png \
   --exclude=make_tar.sh \
   --exclude=launch.sh \
   --exclude=*.Rproj \
   --exclude=data/*.rds \
   --exclude=.* \
   "$src" "$des"
