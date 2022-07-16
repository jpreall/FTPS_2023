### CSHL 2022 Plant Science Course Materials
# Configuring your cloud instance for Scanpy 
-------

# Install Python via Anaconda
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
echo 'export PATH=$HOME/miniconda3/bin:$PATH' >> ~/.bashrc
source ~/.bashrc # alternatively, close and reopen a new terminal
conda init
source ~/.bashrc # alternatively, close and reopen a new terminal
```

# Install Google Chrome
Because the pre-installed Firefox doesn't seem to work for some reason...
```bash
wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
sudo apt install ./google-chrome-stable_current_amd64.deb
```

# Grab my Jupyter notebooks from Github:
```bash
mkdir FTPS22
cd FTPS22
git clone

conda create --name SingleCell python=3.7 --file environment.yml
```
