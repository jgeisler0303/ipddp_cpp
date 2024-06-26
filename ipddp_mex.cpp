#define makeSFunName(c) ipddp_ ## c ## _mex
#define exp_SFunName(c) makeSFunName(c)

#define S_FUNCTION_NAME exp_makeSFunName(PROBLEM_NAME)
#define S_FUNCTION_LEVEL 2

#include <ctime>
#include <stdlib.h>
#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

#define printf mexPrintf

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define FILE_NAME(x) ipddp_problem_ ## x
#define EXPAND_FILE_NAME(x) FILE_NAME(x)
#define FILE_NAME_EXT EXPAND_FILE_NAME(PROBLEM_NAME).hpp

#include TOSTRING(FILE_NAME_EXT)

#include "ipddp.hpp"

#define makeProblemClass(c) Problem ## c
#define exp_makeProblemClass(c) makeProblemClass(c)
#define ProblemClass exp_makeProblemClass(PROBLEM_NAME)


enum input_idx {
    in_idx_x0= 0,
    in_idx_u0,
    in_idx_opt,
    in_idx_last
};
const input_idx in_idx_optional= in_idx_opt;

enum output_idx {
    out_idx_success= 0,
    out_idx_x,
    out_idx_u,
    out_idx_last
};



typedef ipddp<ProblemClass, n_hor> Solver;

bool tryGetOption(double *value, const char *name, const mxArray *mxOptions, int m__=1, int n__=1);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(!(nrhs>=in_idx_optional && nrhs<=in_idx_last)) { mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Wrong number of arguments. Expecting (x0, u0, {opt})"); return; }
    
    if(nlhs!=out_idx_last) { mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Wrong number of return values. Expecting [success, x, u]"); return; }
    
    if(!mxIsDouble(prhs[in_idx_x0]) || mxGetM(prhs[in_idx_x0])!=nx || !(mxGetN(prhs[in_idx_x0])==1 || mxGetN(prhs[in_idx_x0])==n_hor)) { mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Wrong number of elements in 'x0' (%d x 1|%d expected)", nx, n_hor); return; }
    
    if(!mxIsDouble(prhs[in_idx_u0]) || mxGetM(prhs[in_idx_u0])!=nu || mxGetN(prhs[in_idx_u0])!=(n_hor-1)) { mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Wrong number of elements in 'u0' (%d x %d expected)", nu, n_hor-1); return; }

    if(nrhs>in_idx_opt) {
        const mxArray *mxOptions= prhs[in_idx_opt];
        if(!mxIsStruct(mxOptions) || mxGetNumberOfElements(mxOptions)!=1) {
            mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Input options must be a scalar struct.\n");
            return;
        }
    }
    
    
    Solver s;

    if(nrhs>=in_idx_opt) {
        const mxArray *mxOptions= prhs[in_idx_opt];
        double value;
        
        if(tryGetOption(&value, "step_size_factor", mxOptions)) {
            if(value<=0.0 || value>=1.0) {
                mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Option step_size_factor must be >=0 and <=1.\n");
                return;
            }
            s.stepsize_factor= value;
        }
        if(tryGetOption(&value, "tol", mxOptions)) {
            if(value<=0.0) {
                mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Option tol must be positive.\n");
                return;
            }
            s.tol= value;
        }
        if(tryGetOption(&value, "maxiter", mxOptions)) {
            if(value<1.0) {
                mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Option maxiter must be >=1.\n");
                return;
            }
            s.maxiter= value;
        }
        if(tryGetOption(&value, "max_reg_exp", mxOptions)) {
            if(value<1.0) {
                mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Option max_reg_exp must be >=1.\n");
                return;
            }
            s.max_reg_exp= value;
        }
        if(tryGetOption(&value, "infeas", mxOptions)) {
            if(value!=0.0 && value!=1.0) {
                mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Option infeas must 0 or 1.\n");
                return;
            }
            s.infeas= value;
        }
        if(tryGetOption(&value, "max_step_reductions", mxOptions)) {
            if(value<1.0) {
                mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Option max_step_reductions must be >=1.\n");
                return;
            }
            s.max_stepsize_reductions= value;
        }
        if(tryGetOption(&value, "show_progress", mxOptions)) {
            if(value!=0.0 && value!=1.0) {
                mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Option show_progress must 0 or 1.\n");
                return;
            }
            s.show_progress= value;
        }
        if(tryGetOption(&value, "mu", mxOptions)) {
            if(value<0.0) {
                mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Option mu must be positive or 0.\n");
                return;
            }
            s.mu= value;
        }
    }
    
    double *u0= mxGetPr(prhs[in_idx_u0]);
    for(int i= 0; i<n_hor-1; ++i) {
        for(int j= 0; j<nu; ++j) {
            s.get_nominal()[i].u(j)= u0[j + nu*i];
        }
    }
        
    double *x0= mxGetPr(prhs[in_idx_x0]);
    if(mxGetN(prhs[in_idx_x0])==1) {
        for(int j= 0; j<nx; ++j) {
            s.get_nominal()[0].x(j)= x0[j];
        }
        s.initialroll();
    } else {
        for(int i= 0; i<n_hor; ++i) {
            for(int j= 0; j<nx; ++j) {
                s.get_nominal()[i].x(j)= x0[j + nx*i];
            }
        }
    }

    std::clock_t startcputime = std::clock();
    int success= s.solve();
    double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
    
    plhs[out_idx_success]= mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(plhs[out_idx_success])[0]= success;
    
    plhs[out_idx_x]= mxCreateDoubleMatrix(nx, n_hor, mxREAL);
    for(int i= 0; i<n_hor; ++i) {
        for(int j= 0; j<nx; ++j) {
            mxGetPr(plhs[out_idx_x])[j + nx*i]= s.get_nominal()[i].x(j);
        }
    }

    plhs[out_idx_u]= mxCreateDoubleMatrix(nx, n_hor-1, mxREAL);
    for(int i= 0; i<n_hor-1; ++i) {
        for(int j= 0; j<nu; ++j) {
            mxGetPr(plhs[out_idx_u])[j + nu*i]= s.get_nominal()[i].u(j);
        }
    }

}

bool tryGetOption(double *value, const char *name, const mxArray *mxOptions, int m__, int n__) {
    const mxArray *mxOption;
    if((mxOption= mxGetField(mxOptions, 0, name))!=NULL) {
        int m_= mxGetM(mxOption);
        int n_= mxGetN(mxOption);
        if(mxIsSparse(mxOption) || !mxIsDouble(mxOption) || (m_!=m__ && n_!=n__)) {
            mexErrMsgIdAndTxt("IPDDP:InvalidArgument", "Option name '%s' must be a scalar.\n", name);
            return false;
        }
        
        for(int i= 0; i<m__; i++)
            for(int j= 0; j<n__; j++)
                value[i + j*m__]= mxGetPr(mxOption)[i + j*m__];
        
        return true;
    }
    return false;
}
