#include <array>
#include <random>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cassert>

#ifdef MATLAB_MEX_FILE
    #define printf mexPrintf
#endif

template<typename Numeric, typename Generator = std::mt19937>
Numeric random(Numeric from, Numeric to)
{
    thread_local static Generator gen(std::random_device{}());

    using dist_type = typename std::conditional
    <
        std::is_integral<Numeric>::value
        , std::uniform_int_distribution<Numeric>
        , std::uniform_real_distribution<Numeric>
    >::type;

    thread_local static dist_type dist;

    return dist(gen, typename dist_type::param_type{from, to});
}


template <class ipddp_problem, int n_hor, int max_filter_len=11>
class ipddp {
public:
    typedef std::array<ipddp_problem, n_hor> trajectoriesy_t;

    ipddp() {
        filter_len= 1;
        filter_logcost[0]= std::numeric_limits<double>::infinity();
        filter_err[0]= 0.0; 
        
        for(int i= 0;i<trajectories[index_nominal].size()-1; ++i) {
            for(int j= 0;j<nc; ++j) {
                trajectories[index_nominal][i].y(j)= 0.01;
                trajectories[index_nominal][i].s(j)= 0.1;
//               trajectories[index_nominal][i].mu(j)= trajectories[index_nominal][i].y(j)*trajectories[index_nominal][i].s(j);
            }
        }
    }
    
    int solve() {
        if(mu==0.0) {
            mu= cost / (trajectories[index_nominal].size()-1.0) / trajectories[index_nominal][0].s.size();
        }
        
        resetfilter();
        resetreg();
        
        if(show_progress)
            printf("%-12s%-12s%-12s%-12s%-12s%-12s\n","Iteration","mu","Cost","Opt. error","Reg. power","Stepsize");
        for(int iter= 0; iter<maxiter; ++iter)  {
            trajectories[index_nominal][trajectories[index_nominal].size()-1].calc_final_grads(); // px, pxx
            for(int i= trajectories[index_nominal].size()-2; i>=0; --i)
                trajectories[index_nominal][i].calc_grads();
            
            // backward pass with increasing regularisation
            while(true) {
                backwardpass();
                if(!bp_failed || bp_reg>max_reg_exp) break;
            }
            if(bp_failed && bp_reg>=max_reg_exp) {
                if(show_progress)
                    printf("Failed: maximum regularisation reached\n");
                return -1;
            }
            forwardpass();
            if(show_progress)
                printf("%-12d%-12.4g%-12.4g%-12.4g%-12d%-12.3f\n", iter, mu, cost, bp_opterr, bp_reg, stepsize);
            
            if(std::max(bp_opterr, mu)<=tol) {
                if(show_progress)
                    printf("Optimality reached\n");
                return 1;
            }
            
            if(bp_opterr<=0.2*mu) {
                mu= std::max(tol/10.0, std::min(0.2*mu, pow(mu, 1.2)));
                resetfilter();
                resetreg();
            }
        }
        return 0;
    }
    
    void backwardpass() {
        trajectoriesy_t &nominal= trajectories[index_nominal];
        dV1= 0.0;
        dV2= 0.0;
        double c_err = 0.0;
        double mu_err = 0.0;
        double Qu_err = 0.0;

        if(fp_failed || bp_failed) {
            bp_reg++;
        } else if(step == 1) {
            bp_reg--;
        } else if(step <= 4) {
            bp_reg = bp_reg;
        } else {
            bp_reg++;
        }

        if(bp_reg < 0) {
            bp_reg = 0;
        } else if(bp_reg > max_reg_exp) {
            bp_reg = max_reg_exp;
        }

        double V = nominal[nominal.size()-1].p;
        VecX Vx = nominal[nominal.size()-1].px;
        MatXX Vxx = nominal[nominal.size()-1].pxx;
        bp_failed = false;
        for(int i= nominal.size()-2; i>=0; --i) {
            VecX Qx= nominal[i].qx + nominal[i].cx.transpose()*nominal[i].s + nominal[i].fx.transpose()*Vx;
            VecU Qu= nominal[i].qu + nominal[i].cu.transpose()*nominal[i].s + nominal[i].fu.transpose()*Vx; 
            
            MatXX Qxx= nominal[i].qxx + nominal[i].fx.transpose()*Vxx*nominal[i].fx;
            MatXU Qxu= nominal[i].qxu + nominal[i].fx.transpose()*Vxx*nominal[i].fu;
            MatUU Quu= nominal[i].quu + nominal[i].fu.transpose()*Vxx*nominal[i].fu;
            if constexpr(has_fxx) {
                Qxx+= tensdot<VecX, TensXXX, MatXX>(Vx, nominal[i].fxx);
                Qxu+= tensdot<VecX, TensXXU, MatXU>(Vx, nominal[i].fxu);
                Quu+= tensdot<VecX, TensXUU, MatUU>(Vx, nominal[i].fuu);
            }
            if constexpr(has_cxx) {
                Qxx+= tensdot<VecC, TensCXX, MatXX>(nominal[i].s, nominal[i].cxx);
                Qxu+= tensdot<VecC, TensCXU, MatXU>(nominal[i].s, nominal[i].cxu);
                Quu+= tensdot<VecC, TensCUU, MatUU>(nominal[i].s, nominal[i].cuu);
            }
            
            MatCC S= nominal[i].s.asDiagonal();
  
            MatUU Quu_reg= Quu + nominal[i].quu*(pow(1.6, bp_reg)-1.0);
            
            VecC r;
            if(infeas) {
                r= nominal[i].s.cwiseProduct(nominal[i].y);
                r.array()-= mu;
                VecC rhat= nominal[i].s.cwiseProduct(nominal[i].c+nominal[i].y) - r;
                VecC yinv= nominal[i].y.cwiseInverse();
                MatCC SYinv = nominal[i].s.cwiseProduct(yinv).asDiagonal();
        
                Eigen::LLT<MatUU, Eigen::Upper> llt;
                llt.compute(Quu_reg + nominal[i].cu.transpose() * SYinv * nominal[i].cu);
                if(llt.info()!=Eigen::Success) {
                    bp_failed= true;
                    break;
                }
                
                nominal[i].ku= -llt.solve(Qu + nominal[i].cu.transpose()*yinv.cwiseProduct(rhat));
                nominal[i].Ku= -llt.solve(Qxu.transpose() + nominal[i].cu.transpose() * SYinv * nominal[i].cx);
            
                nominal[i].ks = yinv.cwiseProduct(rhat + S * nominal[i].cu * nominal[i].ku);
                nominal[i].Ks = SYinv * (nominal[i].cx + nominal[i].cu * nominal[i].Ku);
                
                nominal[i].ky = -(nominal[i].c + nominal[i].y) - nominal[i].cu * nominal[i].ku;
                nominal[i].Ky = -nominal[i].cx - nominal[i].cu * nominal[i].Ku;

                Quu = Quu + nominal[i].cu.transpose() * SYinv * nominal[i].cu;
                Qxu = Qxu + nominal[i].cx.transpose() * SYinv * nominal[i].cu;
                Qxx = Qxx + nominal[i].cx.transpose() * SYinv * nominal[i].cx;

                Qu = Qu + nominal[i].cu.transpose() * (yinv.cwiseProduct(rhat));
                Qx = Qx + nominal[i].cx.transpose() * (yinv.cwiseProduct(rhat));
            } else {
                r = S * nominal[i].c;
                r.array()+= mu;
                VecC cinv = nominal[i].c.cwiseInverse();
                MatCC SCinv = nominal[i].s.cwiseProduct(cinv).asDiagonal();

                Eigen::LLT<MatUU, Eigen::Upper> llt;
                MatUU H= Quu_reg - nominal[i].cu.transpose() * SCinv * nominal[i].cu;
                llt.compute(H);
                if(llt.info()!=Eigen::Success) {
//                     printf("llt.info(): %d\n", llt.info());
                    bp_failed= true;
                    break;
                }

                nominal[i].ku = -llt.solve(Qu - nominal[i].cu.transpose() * (cinv.cwiseProduct(r)));
                nominal[i].Ku = -llt.solve(Qxu.transpose() - nominal[i].cu.transpose() * SCinv * nominal[i].cx);
                
                nominal[i].ks = -cinv.cwiseProduct(r + S * nominal[i].cu * nominal[i].ku);
                nominal[i].Ks = -SCinv * (nominal[i].cx + nominal[i].cu * nominal[i].Ku);

                Quu = Quu - nominal[i].cu.transpose() * SCinv * nominal[i].cu;
                Qxu = Qxu - nominal[i].cx.transpose() * SCinv * nominal[i].cu;
                Qxx = Qxx - nominal[i].cx.transpose() * SCinv * nominal[i].cx;

                Qu = Qu - nominal[i].cu.transpose() * (cinv.cwiseProduct(r));
                Qx = Qx - nominal[i].cx.transpose() * (cinv.cwiseProduct(r));
            }
            
            dV1+= nominal[i].ku.transpose() * Qu;
            dV2+= (0.5 * nominal[i].ku.transpose() * Quu * nominal[i].ku)(0);
            
            Vx = Qx + nominal[i].Ku.transpose() * Qu + nominal[i].Ku.transpose() * Quu * nominal[i].ku + Qxu * nominal[i].ku;
            Vxx = Qxx + nominal[i].Ku.transpose() * Qxu.transpose() + Qxu * nominal[i].Ku + nominal[i].Ku.transpose() * Quu * nominal[i].Ku;


            // Optimality error
            Qu_err = std::max(Qu_err, Qu.cwiseAbs().maxCoeff());
            mu_err = std::max(mu_err, r.cwiseAbs().maxCoeff());
            
            if(infeas) {
                c_err = std::max(c_err, (nominal[i].c+nominal[i].y).cwiseAbs().maxCoeff());
            }
        } // end loop over nominal

        bp_opterr= std::max({Qu_err, c_err, mu_err}); // std::max(std::max(Qu_err, c_err), mu_err);
    }
    
    void forwardpass() {
        trajectoriesy_t &nominal= trajectories[index_nominal];
        trajectoriesy_t &candidate= trajectories[index_candidate];
        double tau= std::max(0.99, 1.0-mu);

        int failed;
        double pass_cost= 0.0;
        double pass_logcost= 0.0;
        double pass_err= 0.0;
        
        stepsize= 1.0;
        for(step= 1; step<=max_stepsize_reductions; ++step, stepsize*= stepsize_factor) {
            failed= 0;
            
            candidate[0].x= nominal[0].x;
            
            if(infeas) {
                for(int i= 0; i<candidate.size()-1; ++i) {
                    VecX dx= candidate[i].x-nominal[i].x;
                    candidate[i].y = nominal[i].y + stepsize*nominal[i].ky + nominal[i].Ky*dx;
                    candidate[i].s = nominal[i].s + stepsize*nominal[i].ks + nominal[i].Ks*dx;
                    if((candidate[i].y.array()<(1-tau)*nominal[i].y.array()).any() || (candidate[i].s.array()<(1-tau)*nominal[i].s.array()).any()) {
                        failed= 1;
                        break;
                    }
                    
                    candidate[i].u = nominal[i].u + stepsize*nominal[i].ku + nominal[i].Ku*dx;
                    candidate[i+1].x = candidate[i].step();
                }
            } else {
                for(int i= 0; i<candidate.size()-1; ++i) {
                    VecX dx= candidate[i].x-nominal[i].x;
                    candidate[i].s = nominal[i].s + stepsize*nominal[i].ks + nominal[i].Ks*dx;
                    candidate[i].u = nominal[i].u + stepsize*nominal[i].ku + nominal[i].Ku*dx;
                    
                    candidate[i].calc_c();
                    if((candidate[i].c.array()>(1-tau)*nominal[i].c.array()).any() || (candidate[i].s.array()<(1-tau)*nominal[i].s.array()).any()) {
                        failed= 1;
                        break;
                    }
                    candidate[i+1].x = candidate[i].step();
                }
            }
    
            if(failed)
                continue;
            else {
                pass_cost= 0.0;
                for(int i= 0;i<candidate.size()-1; ++i) {
                    pass_cost+= candidate[i].calc_q();
                }
                pass_cost+= candidate[candidate.size()-1].calc_p();
      
                if(infeas) {
                    pass_logcost= 0.0;
                    for(int i= 0; i<candidate.size()-1; ++i)
                        for(int j= 0; j<nc; ++j) 
                            pass_logcost-= std::log(candidate[i].y(j));
                    pass_logcost*= mu;
                    pass_logcost+= pass_cost;
                    
                    pass_err= tol;
                    for(int i= 0; i<candidate.size()-1; ++i) {
                        candidate[i].calc_c();
                        pass_err= std::max(pass_err, (candidate[i].c+candidate[i].y).cwiseAbs().sum());
                    }
                } else {
                    pass_logcost= 0.0;
                    for(int i= 0; i<candidate.size()-1; ++i)
                        for(int j= 0; j<nc; ++j) 
                            pass_logcost-= std::log(-candidate[i].c(j));
                    pass_logcost*= mu;
                    pass_logcost+= pass_cost;

                    pass_err= 0;
                }
      
                for(int i= 0; i<filter_len; ++i)
                    if(pass_logcost>=filter_logcost[i] && pass_err>=filter_err[i]) {
                        failed=2;
                        break;
                    }
                if(failed)
                    continue;
                    
                for(int i= filter_len-1; i>=0; --i)
                    if(pass_logcost<=filter_logcost[i] && pass_err<=filter_err[i]) {
                        // delete values from filter and close gap
                        for(int j= i+1; j<filter_len; ++j) {
                            filter_logcost[j-1]= filter_logcost[j];
                            filter_err[j-1]= filter_err[j];
                        }
                        filter_len--;
                    }
                assert(filter_len<filter_logcost.size());
                filter_logcost[filter_len]= pass_logcost;
                filter_err[filter_len]= pass_err;
                filter_len++;
                
                break;
            }      
        }
        
        if(failed) {
            fp_failed= failed;
            stepsize= 0;
        } else {
            cost= pass_cost;
            logcost= pass_logcost;
            err= pass_err;
            fp_failed= 0;
            
            std::swap(index_nominal, index_candidate);
        }
    }

    template<class type_a, class type_b, class type_res>type_res tensdot(const type_a &a, const type_b &b) {
        type_res res= a(0)*b[0];
        for(int i= 1; i<a.size(); ++i)
            res+= a(i)*b[i];
        return res;
    }
    
    void initialroll() {
        trajectoriesy_t &nominal= trajectories[index_nominal];
        cost= 0.0;
        for(int i= 0;i<nominal.size()-1; ++i) {
            nominal[i+1].x= nominal[i].step();
            nominal[i].calc_c();
            cost+= nominal[i].calc_q();
        }
        cost+= nominal[nominal.size()-1].calc_p();
    }
    
    void resetfilter() {
        trajectoriesy_t &nominal= trajectories[index_nominal];
        if(infeas) {
            logcost= 0.0;
            for(int i= 0; i<nominal.size()-1; ++i)
                for(int j= 0; j<nc; ++j) 
                    logcost-= std::log(nominal[i].y(j));
            logcost*= mu;
            logcost+= cost;
            if(err<tol)
                err= 0.0;
        } else {
            logcost= 0.0;
            for(int i= 0; i<nominal.size()-1; ++i)
                for(int j= 0; j<nc; ++j) 
                    logcost-= std::log(-nominal[i].c(j));
            logcost*= mu;
            logcost+= cost;
            err= 0.0;
        }
        
        filter_len= 1;
        filter_logcost[0]= logcost;
        filter_err[0]= err;
        
        step= 0;
        fp_failed= false;
    }
    
    void resetreg() {
        bp_reg= 0;
        bp_failed= false;
        recovery= 0.0;
    }

    trajectoriesy_t& get_nominal() {
        return trajectories[index_nominal];
    }

    // Options
    double tol= 1e-7;
    int maxiter= 1000;
    int max_reg_exp= 24;
    bool infeas= false;
    int max_stepsize_reductions= 11;
    double stepsize_factor= 0.5;
#ifdef MATLAB_MEX_FILE
    bool show_progress= false;
#else
    bool show_progress= true;
#endif
    
    double cost;
    double logcost;
    double mu= 0;
    double err;
    double bp_opterr;
    int filter_len;
    std::array<double, max_filter_len> filter_logcost;
    std::array<double, max_filter_len> filter_err;
    int step= 0;
    double stepsize;
    bool fp_failed= false;
    bool bp_failed= false;
    int bp_reg= 0;
    double recovery;
    double dV1, dV2;
     
    std::array<trajectoriesy_t, 2> trajectories;
    int index_nominal= 0;
    int index_candidate= 1;
};
