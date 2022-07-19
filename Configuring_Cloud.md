### CSHL 2022 Plant Science Course Materials
# Configuring your cloud instance for Scanpy 
-------

# Install Python via Anaconda
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod a+x Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
#echo 'export PATH=$HOME/miniconda3/bin:$PATH' >> ~/.bashrc
#source ~/.bashrc # alternatively, close and reopen a new terminal
#conda init
source ~/.bashrc # alternatively, close and reopen a new terminal
```

# Install Google Chrome
Because the pre-installed Firefox doesn't seem to work for some reason...
```bash
wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
sudo apt install ./google-chrome-stable_current_amd64.deb
```

# Create a new Conda environment for Scanpy analysis
```bash
# Download my configuration file and create a new Conda environment with it.
wget https://www.dropbox.com/s/c80cw6rykxw27w3/scanpy.yml
conda env create -f scanpy.yml

# Or do it from scratch:
conda create --name scanpy -c conda-forge python=3.8 scanpy python-igraph 
conda activate scanpy
conda install jupyter
conda install -c bioconda gtfparse harmony-pytorch
```

# Grab my Jupyter notebooks from Github:
```bash
wget https://www.dropbox.com/s/ln136h0hqcbn082/FTPS22_data_exploration.ipynb

# Launch the notebook
jupyter notebook FTPS22_data_exploration.ipynb
```
