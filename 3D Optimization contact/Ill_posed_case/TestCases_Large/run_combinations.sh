#!/bin/bash

# Define the values for each parameter
methods=("BFGS" "SLSQP" "L-BFGS-B")
actions=("d" "f" "df")

# Loop through each combination
for method in "${methods[@]}"; do
  for action in "${actions[@]}"; do
    # Run Python script with the current combination of parameters
    /usr/bin/python3 "/home/diego/fem/3D Optimization contact/Ill_posed_case/TestCases_Large/IllConditioned.py" "$method" "$action"
  done
done
