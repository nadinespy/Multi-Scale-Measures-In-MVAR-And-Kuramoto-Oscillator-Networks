#!/bin/bash

# local directory to upload
local_dir="/media/nadinespy/NewVolume1/work/current_projects/mec_experiments/mec_simulations/mec_simulations"

# connect to the SFTP server
sftp neocortex <<EOF

#rsync -av --exclude='*.png' --exclude='*.mat' --exclude='*.git/' -e ssh $local_dir ns508@allocortex:./

lcd $local_dir
put -r *
rm -r *.mat
rm -r *.png
rm -r *.git/

#! find . -type f -name '*.mat' -exec echo rm {} \;
#! find . -type f -name '*.png' -exec echo rm {} \;
#! find . -type f -name '*.git/' -exec echo rm {} \;
EOF

