ampviewer
=========

Check it out here: https://amp.quickgene.net

Contact [Kamil] for the username and password.

[Kamil]: mailto:kslowikowski@gmail.com?subject=ampviewer

Overview
--------

ampviewer is an app for viewing AMP data.

<img height="300px" src="screenshot.png" />

Run the app on your laptop when you're developing new features.

When you're ready to launch the app on the public website, copy it to our Amazon EC2 server:

```bash
src="/PATH/TO/ampviewer/"
des="ig:/srv/shiny-server/ampviewer/"

rsync -avh \
    --exclude=*.sh \
    --exclude=*.Rproj \
    --exclude=data/*.rds \
    --exclude=.* \
    "$src" "$des"
```
