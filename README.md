# FEM - Computational Contact Mechanics and Machine Learning

This repository contains the code used for my PhD thesis, which focuses on two main topics:

1. **Minimization-Based Contact Mechanics**: Implementations of different minimization techniques (BFGS, L-BFGS, and Trust Region methods) to solve contact problems in quasi-static simulations.
2. **Neural Network-Based Contact Detection**: Multi-task neural network models for accelerating contact detection in rigid-body simulations.

## Repository Structure

```
FEM/
│── 1_Minimization_solvers/         # Code for minimization-based contact solvers
│── 2_Accelerated_Contact_Detection/ # Code for ANN-based contact detection
│── PyClasses/                       # Python classes used in the main scripts
│── Meshes/                           # Mesh data used for simulations
│── Other_work/                       # Additional or exploratory work
│── Paper1: Remedy_to_Newton_v0/      # Material related to a paper submission
│── .gitignore                        # Files to be ignored by git
│── environment.yml                    # Conda environment setup
│── fem-env.yml                        # Additional environment configuration
│── requirements.txt                    # Python dependencies
```

## Getting Started

### 1. Clone the Repository
```bash
git clone https://github.com/dhurtadocath/fem.git
cd fem
```

### 2. Set Up the Environment
#### Using Conda:
```bash
conda env create -f environment.yml
conda activate fem-env
```
#### Using pip:
```bash
pip install -r requirements.txt
```

### 3. Running the Code
Each of the main sections (`1_Minimization_solvers/` and `2_Accelerated_Contact_Detection/`) contains 2D and 3D examples that require specific input parameters.

#### **Minimization-Based Contact Mechanics**
- **2D Case:**
```bash
python pseudo2d.py --min_method BFGS --mesh 5 --plastic 0
```
  - Available meshes: `5` (5x5x5), `10`, `15`
  - Minimization methods: `BFGS`, `LBFGSnn` (nn is the number of iteration differences stored), `TR`, `TR-icho` (Trust regions with incomplete Cholesky decomposition as a preconditioner)

- **3D Case:**
```bash
python ContactPotato_Ex1.py --min_method BFGS --mesh 10 --plastic 0
```

#### **Neural Network-Based Contact Detection**
The scripts in `2_Accelerated_Contact_Detection/` follow the same execution principle but include an additional flag to determine if the Multi-Task ANN is used.

- **Example (3D with ANN disabled):**
```bash
python ContactPotato_Ex2.py --min_method BFGS --mesh 15 --plastic 0 --ann 0
```

## Compilation Notice
The material solver in `PyClasses/FEAssembly.py` uses compiled functions to compute the material contributions to the system. **Only the first time you run a script, the code will need to compile this, which may take several minutes.** After that, the compiled files remain in `PyClasses/` and will be reused in subsequent runs, significantly reducing execution time.

## Future Improvements
- Adding **unit tests** for core functionalities.
- Implementing **GitHub Actions** for automated testing.
- Improving documentation with detailed explanations of models and methods.

## License
TBD

