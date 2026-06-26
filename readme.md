# Quantum Phase Interference (QPI) Simulation Toolkit

A comprehensive MATLAB toolkit for simulating and analyzing Quantum Phase Interference patterns using T-matrix formulism in squarenet materials, with particular focus on single and multiple defect systems. Used in Chapter 4.2 in this UBC thesis: [Defect-resolved scattering in quantum materials : a scanning tunneling microscope study with algorithmic multi-channel deconvolution](https://open.library.ubc.ca/soa/cIRcle/collections/ubctheses/24/items/1.0451233)

## Example result illustrations:
**Single defect simulation**
<img width="955" height="896" alt="Ch5_singe_defect" src="https://github.com/user-attachments/assets/25edf864-7e3e-4346-a575-6ae0393827b1" />

**Multi-defect simulation**
<img width="955" height="896" alt="Ch5_multi_defect" src="https://github.com/user-attachments/assets/08a97aaa-6908-4240-9027-f0b03ed55db9" />

**QPI lifetime analysis**
<img width="986" height="450" alt="Ch5_qp_lifetime" src="https://github.com/user-attachments/assets/0dc8edfd-ac30-47cc-8789-019bc0d1f605" />

## Features

- **Local Density of States (LDoS) Calculations**
  - Single defect LDoS computation
  - Multiple defect LDoS simulation
  - T-matrix calculations for various energy ranges
  - Bare Lattice Green's Function computations

- **Noise Analysis & Processing**
  - Multiple noise addition methods (amplitude, power, frequency-domain)
  - SNR and PSNR calculations
  - Wavelet-based and standard deviation noise estimation

- **Visualization Tools**
  - 2D and 3D grid displays of LDoS patterns

  - Energy dispersion visualization
  - QPI pattern analysis
  - Customizable colormap support

- **Analysis Tools**
  - Friedel oscillation analysis
  - Quasiparticle lifetime visualization
  - Convergence testing for numerical parameters

## Core Functions

### LDoS Computation
- `ComputeLDoS.m` - Calculates Local Density of States for single defect
- `computeLDoSWithMultipleDefects.m` - Handles multiple defect calculations
- `computeBLGF.m` - Computes Bare Lattice Green's Function
- `computeTMatrix.m` - Calculates T-matrix for defect scattering

### Data Processing
- `addNoise.m` - Adds controlled noise to LDoS data
- `estimate_noise.m` - Estimates noise levels in experimental data
- `assignDefectLocations.m` - Manages defect position assignment

### Visualization
- `gridDisplay.m` - 2D visualization of LDoS data
- `d3gridDisplay_QPISIM.m` - 3D visualization of LDoS patterns

## Scripts

### Main Simulation
- `QPI_simulation.m` - Primary simulation script for QPI patterns
- `TEST_iteration_num_QPI_simulation.m` - Convergence testing
- `Rahul_Sharma_etal_2021_SI_script.m` - Implementation of methods from referenced paper

### Analysis
- `convolved_vs_multi.m` - Comparison of convolution and multi-defect approaches
- `visualize_QP_lifetime.m` - Quasiparticle lifetime analysis

## Installation

1. Clone this repository
2. Add the `functions/` directory to your MATLAB path
3. Ensure required MATLAB toolboxes are installed:
   - Image Processing Toolbox
   - Signal Processing Toolbox
   - Wavelet Toolbox (optional)

## Usage

Basic example for single defect simulation:

## Requirements

- MATLAB R2020b or newer
- Image Processing Toolbox
- Parallel Computing Toolbox (recommended for performance)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

[Add your chosen license here]

## Citation

If you use this code in your research, please cite:
[Add relevant paper citations here]

## Authors

- Dong Chen

## Contact

[Add your contact information here]
