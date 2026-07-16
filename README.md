# Gaussian Mixture Model (GMM) using EM Algorithm in R

## Project Overview

This project implements the **Expectation-Maximization (EM) Algorithm** for a **Gaussian Mixture Model (GMM)** from scratch using the R programming language.

Unlike built-in libraries, this implementation manually performs every step of the EM algorithm, making it useful for understanding how Gaussian Mixture Models work internally.

---

## Features

- Generates sample data from two Gaussian distributions.
- Implements the EM algorithm from scratch.
- Computes multivariate Gaussian probability density.
- Performs:
  - Expectation (E-Step)
  - Maximization (M-Step)
- Calculates Log-Likelihood for convergence.
- Assigns clusters to data points.
- Visualizes the final clusters.

---

## Algorithm

### 1. Data Generation

Two Gaussian datasets are generated:

- Cluster 1 → Mean = 2
- Cluster 2 → Mean = 7

These datasets are combined into a single dataset.

---

### 2. Parameter Initialization

The algorithm initializes:

- Mixture Weights
- Means
- Covariance Matrices

Random data points are selected as the initial cluster means.

---

### 3. Expectation Step (E-Step)

For every data point, the algorithm computes the probability that it belongs to each Gaussian component.

This probability is called the **responsibility (γ)**.

---

### 4. Maximization Step (M-Step)

Using the computed responsibilities, the algorithm updates:

- Mixture Weights
- Mean Vectors
- Covariance Matrices

---

### 5. Log-Likelihood

The log-likelihood is computed after each iteration.

The algorithm stops when the change in log-likelihood becomes smaller than the specified tolerance.

---

### 6. Cluster Assignment

Each data point is assigned to the Gaussian component with the highest responsibility.

---

## Technologies Used

- R Programming
- Base R Functions

No external libraries are required.

---

## Project Structure

```
.
├── gmm_em_algorithm.R
└── README.md
```

---

## How to Run

1. Install R or RStudio.
2. Clone this repository.

```bash
git clone https://github.com/yourusername/gmm-em-r.git
```

3. Open the R script.

4. Run the complete script.

---

## Output

The program displays:

- Iteration number
- Log-Likelihood
- Final Mixture Weights
- Final Means
- Final Covariance Matrices
- Cluster Assignments

Finally, it plots the clustered data with cluster centers.

---

## Example Output

```
Iteration: 1 LogLikelihood: -625.43

Iteration: 2 LogLikelihood: -542.12

...

Algorithm Converged!

Final Mixture Weights:
0.50
0.50
```

A scatter plot is generated showing:

- Colored clusters
- Black stars representing cluster centers

---

## Learning Objectives

This project helps understand:

- Gaussian Mixture Models (GMM)
- Expectation-Maximization (EM) Algorithm
- Probability Density Functions
- Maximum Likelihood Estimation
- Unsupervised Machine Learning
- Clustering Techniques

---

## Future Improvements

- Support multiple Gaussian components.
- Use real-world datasets.
- Add model selection using BIC/AIC.
- Compare results with the `mclust` package.
- Improve visualization with `ggplot2`.

---

## Visualization

<img width="800" height="600" alt="dataset" src="https://github.com/user-attachments/assets/c6db971e-52dd-4eae-805e-9da5a5a505ba" />
<img width="800" height="600" alt="final_clusters" src="https://github.com/user-attachments/assets/7e845af3-811f-4ea4-8283-df8cd2e428e2" />
<img width="800" height="600" alt="log_likelihood" src="https://github.com/user-attachments/assets/2981d15d-6208-46b7-b4b8-618a665e035f" />
<img width="1000" height="500" alt="true_vs_estimated" src="https://github.com/user-attachments/assets/4f76e099-a3ff-4bdb-9d33-3b080758d7c6" />



