/*
 * solver.cpp
 *
 *  Created on: Dec 1, 2020
 *      Author: d-w-h
 */

#include <stdio.h>
#include <stdlib.h>
#include "main.hpp"
#include "user_types.hpp"

double dra_dCa(double Ca) {
    double dCa = 0.00001;
    return (ra(Ca + dCa) - ra(Ca)) / dCa;
}

void solver(g_params grid_params, p_params physical_params, t_data time_data, s_data* solver_data) {
    double H, W, L, Ug, Ul, kg, kl, K, Cag0, Cal0, D;
    double to, tf;
    int N, Nt;

    /* Parameters */
    N = grid_params.N;
    H = grid_params.H;
    W = grid_params.W;
    L = grid_params.L;

    Ug = physical_params.Ug;
    Ul = physical_params.Ul;
    kg = physical_params.kg;
    kl = physical_params.kl;
    K = physical_params.K;
    Cag0 = physical_params.Cag0;
    Cal0 = physical_params.Cal0;
    D = physical_params.D;

    Nt = time_data.Nt;
    to = time_data.to;
    tf = time_data.tf;

    /* Allocate memory */
    double* Cag = new double[N];
    double* Cago = new double[N];
    double* Cag_prev = new double[N];
    double* Cal = new double[N];
    double* Calo = new double[N];
    double* Cal_prev = new double[N];
    double* x_c = new double[N];

    /* Start calculations */
    double del_x = L / N;
    double del_t = (tf - to) / Nt;
    double B = kg + kl / K;
    double t = to;
    double Calin = 0.0;
    int max_outer_iter = 300;
    int max_inner_iter = 300;

    /* Initialize data */
    for(int i = 0; i < N; ++i) {
        Cag[i] = 0.0;
        Cago[i] = 0.0;
        Cal[i] = 0.0;
        Calo[i] = 0.0;
        Cag_prev[i] = 0.0;
        Cal_prev[i] = 0.0;
        x_c[i] = i*del_x + 0.5*del_x;
    }

    /* Start simulation */
    while(t < tf) {
        /* Start Gauss-Seidel iterations */
        int outer_iter = 0;
        while(outer_iter < max_outer_iter) {
            int inner_iter = 0;
            for(int i = 0; i < N; ++i) {
                Cag_prev[i] = Cag[i];
                Cal_prev[i] = Cal[i];
            }
            /* Inner iterations */
            while(inner_iter < max_inner_iter) {
                /*First node gas phase*/
                Cag[0] = (Ug*H*W*Cag0 + kg*W*del_x*kl*Cal[0]/B + H*W*del_x*Cago[0]/del_t) /
                         (Ug*H*W + kg*W*del_x + H*W*del_x/del_t - kg*W*del_x*kg/(B+1e-20));

                /* 2nd node to last node gas phase */
                for(int i = 1; i < N; ++i) {
                    Cag[i] = (Ug*H*W*Cag[i-1] + kg*W*del_x*kl*Cal[i]/B + H*W*del_x*Cago[i]/del_t) /
                             (Ug*H*W + kg*W*del_x + H*W*del_x/del_t - kg*W*del_x*kg/(B+1e-20));
                }

                /* Inlet node liquid phase */
                Calin = (Ul*Cal0 + D*Cal[0]/(0.5*del_x)) / (Ul + D/(0.5*del_x));

                /* First node liquid phase */
                Cal[0] = (D*H*W*Calin/(0.5*del_x) + D*H*W*Cal[1]/del_x + Ul*H*W*Calin + kl*W*del_x*kg*Cag[0]/(B*K+1e-20) +
                          ra(Cal_prev[0])*H*W*del_x - dra_dCa(Cal_prev[0])*Cal_prev[0]*H*W*del_x + H*W*del_x*Calo[0]/del_t) /
                         (D*H*W/(0.5*del_x) + D*H*W/del_x + Ul*H*W + kl*W*del_x + H*W*del_x/del_t - kl*W*del_x*kl/(B*K+1e-20) - dra_dCa(Cal_prev[0])*H*W*del_x);

                /* Central nodes liquid phase */
                for(int i = 1; i < N - 1; ++i) {
                    Cal[i] = (D*H*W*Cal[i-1]/del_x + D*H*W*Cal[i+1]/del_x + Ul*H*W*Cal[i-1] + kl*W*del_x*kg*Cag[i]/(B*K+1e-20) +
                              ra(Cal_prev[i])*H*W*del_x - dra_dCa(Cal_prev[i])*Cal_prev[i]*H*W*del_x + H*W*del_x*Calo[i]/del_t) /
                             (2*D*H*W/del_x + Ul*H*W + kl*W*del_x + H*W*del_x/del_t - kl*W*del_x*kl/(B*K+1e-20) - dra_dCa(Cal_prev[i])*H*W*del_x);
                }

                /* Last node liquid phase */
                Cal[N-1] = (D*H*W*Cal[N-2]/del_x + Ul*H*W*Cal[N-2] + kl*W*del_x*kg*Cag[N-1]/(B*K+1e-20) +
                            ra(Cal_prev[N-1])*H*W*del_x - dra_dCa(Cal_prev[N-1])*Cal_prev[N-1]*H*W*del_x + H*W*del_x*Calo[N-1]/del_t) /
                           (D*H*W/del_x + Ul*H*W + kl*W*del_x + H*W*del_x/del_t - kl*W*del_x*kl/(B*K+1e-20) - dra_dCa(Cal_prev[N-1])*H*W*del_x);

                inner_iter++;
            }

            outer_iter++;
        }

        /* Set 'old' timestep values */
        for(int i = 0; i < N; ++i) {
            Cago[i] = Cag[i];
            Calo[i] = Cal[i];
        }

        t += del_t;
    }

    /* Set solver results */
    for(int i = 0; i < N; ++i) {
        solver_data->Cag[i] = Cag[i];
        solver_data->Cal = Cal;
        solver_data->x_c[i] = x_c[i];
    }

    /* Deallocate memory */
    delete [] Cag;
    delete [] Cago;
    delete [] Cag_prev;
    delete [] Cal;
    delete [] Calo;
    delete [] Cal_prev;
    delete [] x_c;

}
