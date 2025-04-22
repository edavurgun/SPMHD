# SPMHD: Smoothed Particle Magnetohydrodynamics in Fortran

This repository contains a Fortran implementation of Smoothed Particle Magnetohydrodynamics (SPMHD), a method used to simulate magnetized fluid flows in astrophysical and plasma physics applications.

## Features

- Written in standard Fortran (compatible with modern Fortran compilers)
- Simulates magnetohydrodynamic interactions using the SPH (Smoothed Particle Hydrodynamics) approach
- Designed for research and educational purposes

## Requirements

- A Fortran compiler (e.g., `gfortran`)
- Basic familiarity with numerical simulation methods

## Usage

To compile the code:

```bash
gfortran -O3 -o spmhd SPMHD.f90
