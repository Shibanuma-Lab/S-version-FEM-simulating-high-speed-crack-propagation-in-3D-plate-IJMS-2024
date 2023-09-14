# Overview 
This is the code to perform the generation phase analysis for simulating the high-speed crack propagation in a 3D plate made by polymethyl methacrylate. 

# Description 
The global mesh is generated based on the geometry of the plate specimens and will not change with crack propagation. The local mesh is generated based on	 the shape of crack front, and update with the crack propagation. The history of crack velocity and crack shape is captured by the high-speed camera during the experiment. The experimental data is stored in “S-method-dynamic-crack-propgation-3D-plate/ generation/ experiments_data”. The input files like meshes and boundary conditions are generated using the Python and the all the input files will be moved to the SFEM solver which is developed based on Fortran to perform the analysis. The code to perform the postprocess is also developed based on Python and using the output files from the SFEM solver, the postprocess code can calculate the J-integral and stress intensity factor.

## Install Python 3.10 & Pipenv
```
sudo apt update
sudo apt upgrade
sudo apt install build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libsqlite3-dev libreadline-dev libffi-dev libbz2-dev
wget https://www.python.org/ftp/python/3.10.6/Python-3.10.6.tgz
tar -xf Python-3.10.6.tgz
cd Python-3.10.6
./configure --enable-optimizations
make -j$(nproc)
sudo make altinstall
```
**Verify Python installation: Check that Python 3.10 is installed correctly by running:**
```
python3.10 --version
```
**Install Pipenv:**
```
python3.10 -m pip install --user pipenv
```
**Open the .bashrc file in your home directory using a text editor, and add the following line at the end of the file:**
```
export PATH="$HOME/.local/bin:$PATH"
```
**Update the current terminal session: To apply the changes made to .bashrc, either open a new terminal window or run the following command:**
```
source ~/.bashrc
```
**Verify Pipenv installation: Confirm that Pipenv is installed correctly by running:**
```
pipenv --version
```

## Clone the code
```
git clone https://github.com/Shibanuma-Lab/S-method-dynamic-crack-propgation-3D-plate
cd S-method-dynamic-crack-propgation-3D-plate
```
**Activate the virtual environment:**
```
pipenv --python /usr/local/bin/python3.10
pipenv install
```
**Install modules:**
```
chmod +xr ./install_lib.sh
./install_lib.sh
```

## Run the code
Run main.py from the generation directory.
```
pipenv run python3 main.py --step_start 0 --step_end 300
```
The command is below.
```
--h, --help  Show help message and text.
--test_start  Specify first test number. (default: 1)
--test_end  Specify last test number. (default: 1)
--step_start  Specify first step number. (default: 0)
--step_end  Specify first test number. (default: 300)
--delete  If True, delete ALL files in /inputfiles and /results. (default: False)
--particular  If True, this program runs in the specified step. (default: False)
--is_jonly  If True, calculate only J-integral. (default: False)
--is_meshonly If True, calculate only meshs. (default: False)
```

## Demo

## Licence

[MIT](https://github.com/tcnksm/tool/blob/master/LICENCE)

## Reference
[1]User Name, 'Paper Titile' Conference Name pp.xx 20XX

[tcnksm](https://github.com/tcnksm)

cf. [how to write readme](https://deeeet.com/writing/2014/07/31/readme/)
