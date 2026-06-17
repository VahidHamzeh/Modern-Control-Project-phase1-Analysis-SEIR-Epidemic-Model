![Control Theory](https://img.shields.io/badge/Domain-Modern_Control-blue?style=for-the-badge&logo=mathworks&logoColor=white)
![Model](https://img.shields.io/badge/Model-SEIR-red?style=for-the-badge&logo=microchip&logoColor=white)
![Tool](https://img.shields.io/badge/Tool-MATLAB-orange?style=for-the-badge&logo=matlab&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge&logo=github&logoColor=white)

# SEIR Epidemic Model Analysis & Control

This repository contains the Phase 1 project for the **Modern Control Engineering** course at **Sharif University of Technology (SUT)**. The project focuses on the mathematical modeling, linearization, structural analysis, and simulation of a fractional-order SEIR epidemic model, based on the research by *Soulaimani & Kaddar (2023)*.

## Project Overview
The primary goal of this project is to transition from theoretical control concepts to practical engineering applications. Key tasks performed include:
* **Modeling:** Extraction of non-linear differential equations and formulation into standard state-space form ($\dot{x}=f(x,u)$).
* **Linearization:** Jacobian-based linearization at the Disease-Free Equilibrium (DFE).
* **Structural Analysis:** Assessing controllability and observability using five distinct methods (Matrix Rank, PBH Test, Jordan Form, Gramian, and SVD).
* **Canonical Forms:** Conversion to Jordan, Controller, and Observer Canonical Forms.
* **Validation:** Comparative simulation between the non-linear and linearized models using Step and Sinusoidal inputs, with error analysis (RMSE & NRMSE).

## Repository Structure
* `/report`: Contains the final Phase 1 scientific report (`SEIR_Phase1_FINAL_v3.pdf`).
* `/matlab`: Contains all simulation scripts, numerical Jacobian implementations, and structural analysis functions.
    * `SEIR_Phase1_Final.m`: Main script for model definition and simulation.
    * `linearize_SEIR.m`: Function for numerical Jacobian calculation.

## Key Results
* **Stability:** The system is proven to be asymptotically stable under the control case ($v=4.0$, $R_0 < 1$).
* **Structural Analysis:** The system is found to be **partially controllable** ($rank=1/4$) and **partially observable** ($rank=2/4$), as confirmed by all five analytical methods.
* **Accuracy:** The linearized model provides high accuracy near the equilibrium point, with $NRMSE \approx 0.0012$.

## Tools & Frameworks
* **Modeling & Simulation:** MATLAB & Simulink
* **Documentation:** LaTeX (XeLaTeX)
* **Analysis:** Modern Control Theory (State-Space, Stability Analysis)

## Credits
* **Authors:** Vahid Hamzeh Garkani, Amirreza Farzaneh
* **Course:** Modern Control Engineering
* **Instructor:** Dr. Pourshamsi
* **Institution:** Sharif University of Technology (SUT)

---
*For more technical details, please refer to the [full report in the /Report folder](./Report/SEIR_Phase1_FINAL_v3.pdf).*
