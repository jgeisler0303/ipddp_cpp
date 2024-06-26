function mex_ipddp_problem(problem_name, problem_func)
ipddp_problem_cpp(problem_name, problem_func)

mex(['-DPROBLEM_NAME=' problem_name], 'CXXFLAGS="-std=c++17 -fPIC"', '-output', ['ipddp_' problem_name '_mex'], 'ipddp_mex.cpp')