Run the program for evaluation?
=====================
1. Simply run make - This should create five executables, one each for Serial, OpenMP, MPI, MPI+OpenMP,and CUDA.
2. Run the executable you wish to test. This will generate 7 text files.
    2a. If you want to test the MPI executables then you will have to submit a batch. Use the ofrun batch file to submit the case.
3. Run the Readfiles.m matlab file from the same folder as the text files.
4. Run the matlab file to visualize the vortex velocity advection test case.
    4a. If you face an error. Just check the path provided in lines 7-9 for the text file. If you are working on a windows machine,
        you will have to change line 8 with a backward slash.
5. The test case set up in the files currently is different from the case that was used to report the timing results in the report.
   This change was made because the test case used for timing analysis was much larger and takes longer to complete.
6. Important Note: Running the code generates a file "details.txt". Currently, there is an issue with this file from MPI codes.
   If you wish to check the MPI code, run the serial code first. Copy this details.txt somewhere else and then run the MPI code.
   Use the text files generated from this code with the originally generated details file with the Matlab code to visualize.



How to run the program for other test cases?
=====================

For the serial program (GALS_Advection_Serial)
1. Check the Constants.h file - See if all the parameters for the level set and velocity field are correctly defined.
2. Check the ******Velocity.h file - See if the correct velocity field is implemented
3. Check the InitializeLevelSet.h file - See if the correct level set field is initialized
4. In GALS_Advection.cpp - Modify the user setting parameters that are in the beginning of the main program.
5. Change in the include file name in the main program GALS_Advection.cpp
6. Run the Makefile to generate the executable GALS_Advection_Serial


For the MPI+OpenMP program (GALS_Advection_MPI_OMP)
1. Do the same steps 1-3 as serial program
2. Make changes to Setting.h to change the program settings
3. Change in the include file name in the main program GALS_Advection_MPI_OpenMP.cpp
4. Run the Makefile to generate the executable GALS_Advection_MPI_OMP


For the MPI program (GALS_Advection_MPI)
1. Same instructions as the MPI+OpenMP program, and the executable is GALS_Advection_MPI


For the OpenMP program (GALS_Advection_OMP)
1. Same instructions as the MPI+OpenMP program, and the executable is GALS_Advection_OMP


For the CUDA program (GALS_Advection_CUDA)
The modularity of this program is worse because of problems with scoping of certain functions and variables.
1. Go to the CUDA folder
2. Update the Initialize Level Set and Velocity files as in previous implementations
3. Change in the include file name in the main program GALS_Advection.cu and in Allocation.h
4. Update the first function of Allocation.h with the new constant values as in InitializeLevelSet.h
5. Now do Make and generate the executable GALS_Advection_CUDA
