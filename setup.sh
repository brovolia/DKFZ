#!/bin/bash


mkdir -p $HOME/Shiny_for_cell
cd $HOME/Shiny_for_cell
git clone https://github.com/brovolia/DKFZ.git


# download necessary R envoronment from the disk.yandex repository
wget https://disk.yandex.com/d/q12WxWEjkOaO0A


# install R libraries and dependencies
./install_required_packages.R
