Contents:

'dbnsPotentialFoam' : Most up-to-date code; compiles as a library instead of
  modifying the main foam-extend source tree
'extrapolatedGradient' : code from internet for a gradient matched BC
  (didn't end up using this, but could be useful)
'in_source_reconstruction' : older version of code, which modifies the
  foam-extend source files
'new_solver_instructions.txt' : quick guide to how the new solver files were
 initially created for modification


Compiling:

***** All steps assume foam-extend is properly installed *****

1. Navigate to the 'dbnsPotentialFoam/balancedReconstruction' directory
2. Run command 'wmake libso' to generate the balancedFlux code library
3. Navigate up to the 'dbnsPotentialFoam' directory
4. Run command 'wmake' to generate the solver executable
5. Use the command 'dbnsPotentialFoam' to invoke the new solver
6. See cases for examples of 'fvScheme' and 'fvSolution' files

***** swak4foam is also useful for groovyBC and other capabilities *****