#!/bin/bash

echo "Creating Conda environment for PanSSRAtor..."
conda create --name PanSSRAtor python=3.9 -y

echo "Activating environment..."
source activate PanSSRAtor

echo "Installing dependencies..."
conda install -c bioconda pysam primer3 intervaltree -y
conda install -c conda-forge numpy pandas tqdm -y

# Installing pyfastx separately due to compatibility issues
pip install pyfastx

# Installing alternative for 'tre' since it's problematic
echo "Installing regex handling libraries..."
pip install regex

# Optional: Check if installation is successful
echo "Installation complete. To activate, use: conda activate PanSSRAtor"

