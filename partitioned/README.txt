==================================================================================================

Directory Structure: 

	fortran 	- Contains implementations of EPBM, EAB, ESDC, IMEXSDC, and ETDRK4 for
			  unpartitioned problems. fortran/equations contains relevant equation
			  files.

	results		- Code outputs results to specified path located in this directory.

==================================================================================================	

Compiling Instructions:	

	A Makefile are located in both the fortran directory. The makefile accepts the option
	OMP=TRUE to activate -fopenmp flag The makefile has the following targets:

	sisc-experiments.exe	- compiles program experiments.f90 into which reproduces one of 
				  the four partitioned numerical experiments included in the paper. 

	Assuming Dependencies are satisfied (See makefile for details), typical usage might be:

		make OMP=TRUE
		./sisc-experiments 	  

==================================================================================================

Solving Different Equations

	To change the equation being solved, modify line 9 of sisc-experiments.f90
	and change equationname_mod.f90 to desired equation module. To generate new equation 
	modules see README.txt inside equation folder.
	
==================================================================================================

Figure Data:

	Performance data for numerical schemes is collected in Fortran and saved in text files in
	the folder results.