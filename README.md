Also, let's mention that the material solver in PyClasses/FEAssembly.py uses compiled functions to compute the material contributions to the system so, ONLY the first time we run a script, the code will have to compile this, which may take several minutes. After that the compiled files stay in PyClasses and those will be used any other time that the code is run.

Also, each part ("1_..." and "2_...") contain 2d and 3d examples. Running this examples require inputs so they have to be run like:

part 1, 2d case: "python pseudo2d.py --min_method BFGS --mesh 5 --plastic 0" # here the meshes are 5 (5x5x5), 10 or 15. Minimization methods BFGS, LBFGSnn  (nn is the number of iteration differences stored), TR, TR-icho (Trust regions incomplete-cholesky decomposition as preconditioner)

part1, 3d case: "python ContactPotato_Ex1.py --min_method BFGS --mesh 10 --plastic 0" # there's a similar principle

For part 2 (In "2_Accel.../"), there are the same scripts but in addition they include whether they are solved using the MultiTask ANN or not.

For example (3D): python ContactPotato_Ex2.py --min_method BFGS --mesh 15 --plastic 0 --ann 0


