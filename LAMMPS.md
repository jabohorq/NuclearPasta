# LAMMPS in Vagrant
## Vagrant file (Vagrantfile)
Using Ubuntu 20.04
```
Vagrant.configure("2") do |config|
  config.vm.box = "bento/ubuntu-20.04"

  #config.vm.synced_folder "./data", "/vagrant_data"

end
```
## Start VB
```
vagrant up
vagrant ssh
vagrant up
vagrant ssh
sudo passwd root
(su -)
sudo apt-get update -y
```
## Pyhton 3 by default
```
sudo update-alternatives --install /usr/bin/python python  /usr/bin/python3 1
python --version
```

## Install dependencies
libssl-dev is the openssl dev library (need it for cmake)
```
sudo apt-get install -y wget build-essential ssh zlib1g-dev libfftw3-dev \
libopenblas-dev libopenmpi-dev libssl-dev python-dev 
sudo rm -rf /var/lib/apt/lists/*
```

## Install CMAKE
```
wget https://github.com/Kitware/CMake/releases/download/v3.19.6/cmake-3.19.6.tar.gz
tar -zxvf cmake-3.19.6.tar.gz
cd cmake-3.19.6
./bootstrap
make
make install
```

 # Install CCMAKE
 Edit cmake file to enable different lammps packages
 ```
 sudo apt-get install cmake-curses-gui
 ```

## Build LAMMPS
```
git clone -b release https://github.com/lammps/lammps.git
cd lammps
mkdir build
cd build
cmake -C ../cmake/presets/basic.cmake \
-D CMAKE_INSTALL_PREFIX=$HOME/.local \
-D BUILD_SHARED_LIBS=on -D LAMMPS_EXCEPTIONS=on -D PKG_PYTHON=on \
-D PKG_EXTRA-PAIR=on ../cmake
cmake --build .
cmake --install .
```
Executable (lmp) on $HOME/.local/bin
And libraries on $HOME/.local/lib
Export path
Edit ~/.bashrc
```
export PATH=$PATH:$HOME/.local/bin 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib 
```
Source bash file
```
source ~/.bashrc
```
