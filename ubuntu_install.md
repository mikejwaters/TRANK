# Step 0: to make our lives better, install virtual box guest additions

sudo apt-get install virtualbox-guest-utils  virtualbox-guest-x11 virtualbox-guest-dkms


# Step 1: install anaconda
## download their install script
cd ~/Downloads
bash Anaconda3-*-Linux-x86_64.sh 
## this name might need to be changed if they change their file names
## just use the defaults they give you, and make sure to prepend your path when asked. 
## Don't install microsoft visual studio, it's not needed


# Step 2:  install tmm
pip install tmm

# Step 3: install git
sudo apt-get install git

# Step 4: make personal python package directory and clone trank there
mkdir ~/my_python
cd ~/my_python
git clone https://github.com/mikejwaters/TRANK.git

# step 5: add your personal directory to your python path
echo 'export PYTHONPATH="${PYTHONPATH}:$HOME/my_python"' >> ~/.bashrc

# step 6: Open a new terminal to test if it's working and
## cd to the examples to run one, this will take some time and show you many plots
cd ~/my_python/TRANK/examples/example_1_TMM_and_back/
python 1_tmm_film_on_substrate.py && python 2_TRANK_fit_spectra.py  && python 3_compare_nk.py && python 4_compare_spectra.py


