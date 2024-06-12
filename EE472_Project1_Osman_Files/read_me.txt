#### Description of the EE472 Project
#### Power Flow Analysis with Newton-Raphson Iteration using Phyton

## Overview
This project implements a power flow solver using the Newton-Raphson iteration method in Python. 
The power flow solver calculates the steady-state voltages and phase angles in a power system network. 
It is commonly used in power system analysis to determine the operating conditions of the network, 
including power flows, voltage profiles, and losses.

## Features
- Power flow calculation using the Newton-Raphson iteration method
- Support for both radial and interconnected power system networks
- IEEE input format for bus and branch data
- Visualization of power flow results (optional)

## Requirements
- Python 3.x
- NumPy (for numerical computations)
- Matplotlib (for visualization)
- EE472_P1_OO_main.py file and ieee300cdf.txt should be in the same folder
- ieee300cdf.txt file must be IEEE common data format.

## Installation
1. Download the repository from [GitHub](https://github.com/odtu/project-1-OsmanOzutemiz).
2. Install Microsoft VSCode (recommended)
3. Install the required Python packages: pip install numpy matplotlib time


## Usage
1. Prepare input data:
- Create a text file containing bus data and branch data of the power system network.
- Ensure the data format complies with the IEEE specified input format.
2. Run the power flow solver:
Replace `input_data.txt` with the path to your input data file.
3. View the power flow results:
- The solver will display the calculated voltages and phase angles for each bus in the network.
- Optionally, visualize the power flow results using Matplotlib.

## Input Data Format
See the provided `cdf_readme.txt` file for a input data format explanation.

## Sample Input Data
See the provided `ieee300cdf.txt` file for a sample input data format.

## License
This project is licensed under the ODTU License. 

## Author
Osman Özütemiz
osmanozutemiz.blogspot.com
Email: osman.ozutemiz@metu.edu.tr
LinkedIn: [Osman Özütemiz](https://www.linkedin.com/in/osmanozutemiz/)

