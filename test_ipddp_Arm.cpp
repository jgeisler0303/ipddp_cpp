#include "ipddp_problem_Arm.hpp"
#include "ipddp.hpp"

typedef ipddp<ProblemArm, n_hor> sol;

void init(sol &s) {
    s.get_nominal()[0].x(0)= -M_PI/2.0; // [-pi/2;0;0];
    s.get_nominal()[0].x(1)= 0.0;
    s.get_nominal()[0].x(2)= 0.0;

    for(int i= 0;i<s.get_nominal().size()-1; ++i) {
        s.get_nominal()[i].u(0)= 0.0; // 0.02*random<double>(-0.01, 0.01);
        s.get_nominal()[i].u(1)= 0.0; // 0.02*random<double>(-0.01, 0.01);
        s.get_nominal()[i].u(2)= 0.0; // 0.02*random<double>(-0.01, 0.01);        
    }
}

int main(int argc, char* argv[]) {
    sol solver;
    
    init(solver);
    solver.initialroll();
    solver.solve();    
}
