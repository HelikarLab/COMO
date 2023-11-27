---
title: Obtaining a Gurobi License
permalink: /como_docker_gurobi_license.html
summary: This is how to obtain a Gurobi license
sidebar: como_sidebar
last_updated: August 29, 2023
---

## Obtaining a Gurobi License

From
the [MathWorks Website](https://www.mathworks.com/products/connections/product_detail/gurobi-optimizer.html#:~:text=The%20Gurobi%20Optimizer%20is%20a,implementations%20of%20the%20latest%20algorithms.):
> The Gurobi Optimizer is a state-of-the-art solver for mathematical programming. THe solvers in the Gurobi Optimizer
> were designed from the ground up to exploit modern architectures and multicore processors, using the most advanced
> implementations of the latest algorithms

This is to say, Gurobi will help when attempting to calculate linear programming solutions done during flux balance
analysis calculations

A Gurobi license is free for academic use. To obtain a license, perform the following:

1. Access the [Web License Manager](https://license.gurobi.com/manager/keys), and log in (or create an account)
2. Click on "API Keys" on the left sidebar
3. Click "Create API Key"
4. Enter the relevant information. The "Application Name" can be anything; I have mine set to "COMO", with an empty
   description
5. Click "Create"
6. You **MUST** download the license file now, as it is not available later. Creating new licenses is easily done,
   however.
7. Save this file to a location on your computer, and note the location. This will be used later.
