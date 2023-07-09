#!/bin/bash

git submodule update --init --recursive
cd sfem_linear
git checkout -b furuhashi
git pull origin furuhashi
chmod +xr install_lib.sh
install_lib.sh
make

cd ..
cd generation
mkdir logs
mkdir Newton
mkdir inputfiles