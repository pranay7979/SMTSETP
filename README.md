# Single Machine Total Stepwise Earliness and Tardiness Problem (SMTSETP)

This repository contains the source code and experimental framework for the **Single Machine Total Stepwise Earliness and Tardiness Problem (SMTSETP)**, a novel combinatorial optimization problem where both early and late job completions are penalized via discrete, stepwise cost functions.

The project implements and compares two evolutionary metaheuristics to minimize the total stepwise deviation cost across all jobs.

##  Problem Overview
The SMTSETP extends traditional single-machine scheduling frameworks to match real-world scenarios (such as cold-chain logistics, pharmaceutical manufacturing, and perishable goods production) where inventory and shipping costs escalate at distinct time thresholds. 

The objective notation of the problem is:
$$1||\sum_{i=1}^{n}\overline{f_{i}}(C_{i})$$

---

##  Core Algorithms Implemented
The codebase is implemented in **C** and evaluates two steady-state genetic algorithm (SSGA) variants:

1. **Approach I (SSGA):** Steady-State Genetic Algorithm featuring an initial population generated entirely at random, binary tournament selection, position-based crossover, a three-point swap mutation, and an **Adjacent Pair-wise Interchange (API)** local search technique.
2. **Approach II (HSSGA):** Heuristic-Based Steady-State Genetic Algorithm that utilizes identical genetic operators but improves convergence by seeding the initial population with high-quality constructive heuristics:
   * **Moore Heuristic** (Tardiness minimization)
   * **Reverse Moore Heuristic** (A novel approach designed to control excessive earliness penalties)
   * **NEH-based Insertion Heuristic**

Computational and statistical results (Wilcoxon signed-rank test) demonstrate that **HSSGA** significantly outperforms the purely random initialization of SSGA, particularly on larger instances.

---

## 🛠️ Getting Started

### Prerequisites
* A standard C compiler (e.g., `gcc`)
* Make or an IDE capable of handling multi-file C configurations

### Compilation
Compile the implementation directly using `gcc`:
```bash
gcc -O3 main.c algorithms.c utils.c -o smtsetp