This is a C++ reimplementation of the IPDDP (Interior Point Differential Dynamic Programming) algorithm originally written in MATLAB by @xapavlov (https://github.com/xapavlov/ipddp).
The code is written as a header-only library using class templates. The problem description (system dynamics, cost function, constraints and gradiens) must be supplied as a class following certain requirements.

Example problem descriptions are automatically generated from the original MATLAB examples. A mex function wrapper and a MATLAB script that automatically compiles the mex function for a given problem description are also included in this repo.

A full demonstration comparing the result of the original code with the C++ implementation is given in the ``test_ipddp`` script.

I tested the code on linux using g++.
