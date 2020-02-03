This repository contains R implementation of the algorithms proposed in "SpMKL: Sparse multiple kernel learning via optimal kernel selection with application to cancer biology" (manuscript under review).

* 1_run_SpMKL_Step1.R : shows how to solve the first step SpMKL problem based on the TCGA cohorts and Hallmark pathways to produce significance values for the second step.
* 2_produce_kernel_significance_values_from_step1.R : produces the significance matrix by combining the results obtained in step 1
* 3_run_SpMKL_Step2.R : shows how to perform the final step of SpMKL on the TCGA cohorts using the significance matrix produced in the previous step.

SpMKL methods
------------
* classification_helper.R => helper functions
* solve_classification_models_cplex.R => support vector machine classification solver, and the cutting plane model using CPLEX optimization software
* solve_classification_models_mosek.R => support vector machine classification solver, and the cutting plane model using Mosek optimization software
* spmkl_train.R => training procedures used in different steps and phases of SpMKL
* spmkl_test.R => test procedure for SpMKL


