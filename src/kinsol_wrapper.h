#pragma once

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */

#include <vector>

////////////////////////////////////////////////////////////////////////////////

int FuncKinsol(N_Vector yvec, N_Vector fvec, void *user_data);

////////////////////////////////////////////////////////////////////////////////

class kinsol_wrapper {


//////////////////// DATA MEMBERS //////////////////////

public: 

    void           *user_data;
    size_t          nvar;

    SUNContext      sun;
    N_Vector        vars;
    N_Vector        scales_v;
    N_Vector        scales_f;
    N_Vector        constraints;
    void           *kmem;
    SUNMatrix       J;
    SUNLinearSolver LS;
    int             rv;
    int             solver_type;                // KIN_NONE, KIN_LINESEARCH, KIN_FP, KIN_PICARD
    int             exact_or_modified_newton;   // 1 or 0, respectively

//////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public: 

kinsol_wrapper(
        void                      *user_data_, 
        const size_t               nvar_,
        const std::vector<double> &scales_v_,
        const std::vector<double> &scales_f_,
        const std::vector<double> &constraints_,
        const double               ftol = 1E-5,
        const double               stol = 1E-5) :
    user_data(user_data_),
    nvar(nvar_) {

        solver_type              = KIN_LINESEARCH;
        exact_or_modified_newton = 1;

        rv = SUNContext_Create(NULL, &sun);

        vars        = N_VNew_Serial(nvar, sun);
        scales_v    = N_VNew_Serial(nvar, sun);
        scales_f    = N_VNew_Serial(nvar, sun);
        constraints = N_VNew_Serial(nvar, sun);

        J  = SUNDenseMatrix(nvar, nvar, sun);   // linear solver matrix J

        for(size_t k=0; k<nvar; k++) {
            NV_Ith_S(scales_v, k)    = scales_v_[k];
            NV_Ith_S(scales_f, k)    = scales_f_[k];
            NV_Ith_S(constraints, k) = constraints_[k];
        }

        kmem = KINCreate(sun);

        rv = KINSetUserData(     kmem, user_data);
        rv = KINSetConstraints(  kmem, constraints);
        rv = KINSetScaledStepTol(kmem, stol);
        rv = KINSetFuncNormTol(  kmem, ftol);
        rv = KINSetMaxSetupCalls(kmem, exact_or_modified_newton);

}

//--------------

    ~kinsol_wrapper(){
        N_VDestroy(vars);
        N_VDestroy(scales_v);
        N_VDestroy(scales_f);
        N_VDestroy(constraints);
        KINFree(&kmem);
        SUNMatDestroy(J);
        SUNLinSolFree(LS);
        SUNContext_Free(&sun);
    }

//////////////////// MEMBER FUNCTIONS /////////////////

public:

int solve(std::vector<double> &y) {

    //---------

    for(size_t k=0; k<nvar; k++)
        NV_Ith_S(vars, k) = y[k];

    rv = KINInit(kmem, FuncKinsol, vars);
    LS = SUNLinSol_Dense(vars, J, sun);  // set linear solver
    rv = KINSetLinearSolver(kmem, LS, J);                 // associate matrix J with solver LS

    //---------

    rv = KINSol(kmem, vars, solver_type, scales_v, scales_f);

    //---------

    for(size_t k=0; k<nvar; k++)
        y[k] = NV_Ith_S(vars,k);

    return rv;
}

};
