
# SIR_MODEL

This repository contains a simulation and analysis of the **SIR (Susceptible-Infected-Recovered)** epidemiological model. The SIR model is a classic mathematical model used to understand the spread of infectious diseases in a population.

## Contents

- `sir_model.html` — An interactive HTML file visualizing the SIR model dynamics.
- `sir_model.py` — Python script used to generate the simulation.
- `README.md` — Project overview and usage instructions.

## About the SIR Model

The SIR model divides the population into three compartments:
- **S (Susceptible):** Individuals who can contract the disease.
- **I (Infected):** Individuals who have the disease and can transmit it.
- **R (Recovered):** Individuals who have recovered and are immune.

The model is governed by the following differential equations:

\[
\begin{aligned}
\frac{dS}{dt} &= -\beta \frac{SI}{N} \\
\frac{dI}{dt} &= \beta \frac{SI}{N} - \gamma I \\
\frac{dR}{dt} &= \gamma I
\end{aligned}
\]

Where:
- \( \beta \) is the transmission rate
- \( \gamma \) is the recovery rate
- \( N \) is the total population

## How to Use

1. Clone the repository:
   ```bash
   git clone https://github.com/CesarC314159/SIR_MODEL.git
   cd SIR_MODEL
   ```
