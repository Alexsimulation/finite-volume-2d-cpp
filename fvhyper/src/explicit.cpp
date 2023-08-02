/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Explicit solver sources
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#include <fvhyper/explicit.h>
#include <fvhyper/post.h>
#include <array>
#include <chrono>



namespace fvhyper {


namespace solver {
    const double limiter_k_value = 10.;
}



void smooth_residuals(
    std::vector<double>& qt_,
    std::vector<double>& smoother_qt,
    std::vector<double>& smoother,
    mesh& m
) {
    /*
        Smooth the residuals in r implicitly using jacobi iteration
    */
    const uint iters = 2;
    const double epsilon = 0.6;

    for (int i=0; i<smoother_qt.size(); ++i) {
        smoother_qt[i] = qt_[i];
    }
    for (int jacobi=0; jacobi<iters; ++jacobi) {
        for (int i=0; i<smoother.size(); ++i) {
            smoother[i] = 0.;
        }
        for (int e=0; e<m.edgesCells.cols(); ++e) {
            const uint& i = m.edgesCells(e, 0);
            const uint& j = m.edgesCells(e, 1);
            for (uint k=0; k<vars; ++k) {
                smoother[vars*i+k] -= qt_[vars*j+k] * epsilon;
                smoother[vars*j+k] -= qt_[vars*i+k] * epsilon;
            }
        }
        for (int i=0; i<m.nRealCells; ++i) {
            const double ne = m.cellsIsTriangle[i] ? 3 : 4;
            for (uint k=0; k<vars; ++k) {
                qt_[vars*i+k] = (smoother_qt[vars*i+k] - smoother[vars*i+k])/(1. + ne * epsilon);
            }
        }
    }
}



void gradient_for_diffusion(
    double* gradx, double* grady, 
    const double* gxi, const double* gyi,
    const double* gxj, const double* gyj,
    const double* qi, const double* qj, 
    const double* ci, const double* cj
) {
    // Evaluate gradient for diffusive fluxes

    // Distance from i to j
    double tij[2];
    tij[0] = ci[0] - cj[0];
    tij[1] = ci[1] - cj[1];

    const double lij = sqrt(tij[0]*tij[0] + tij[1]*tij[1]);

    tij[0] = tij[0] / lij;
    tij[1] = tij[1] / lij;

    // Directional derivative
    double grad_dir[vars];
    for (uint i=0; i<vars; ++i) {
        grad_dir[i] = (qi[i] - qj[i])/lij;
    }

    // Arithmetic average gradient
    double gradx_bar[vars];
    double grady_bar[vars];
    for (uint i=0; i<vars; ++i) {
        gradx_bar[i] = (gxi[i] + gxj[i])*0.5;
        grady_bar[i] = (gyi[i] + gyj[i])*0.5;
    }

    // Corrected gradient
    for (uint i=0; i<vars; ++i) {
        const double grad_dot = gradx_bar[i]*tij[0] + grady_bar[i]*tij[1];
        gradx[i] = gradx_bar[i] - (grad_dot - grad_dir[i])*tij[0];
        grady[i] = grady_bar[i] - (grad_dot - grad_dir[i])*tij[1];
    }
}


void calc_gradients(
    std::vector<double>& gx,
    std::vector<double>& gy,
    const std::vector<double>& q,
    mesh& m
) {
    // reset gradients to be null
    for (uint i=0; i<gx.size(); ++i) {
        gx[i] = 0.;
        gy[i] = 0.;
    }
    
    // Update gradients using green gauss cell based
    for (uint e=0; e<m.edgesNodes.cols(); ++e) {
        const auto& i = m.edgesCells(e, 0);
        const auto& j = m.edgesCells(e, 1);
        const auto& nx = m.edgesNormalsX[e];
        const auto& ny = m.edgesNormalsY[e];
        const auto& le = m.edgesLengths[e];

        const double dxif = m.edgesCentersX[e] - m.cellsCentersX[i];
        const double dyif = m.edgesCentersY[e] - m.cellsCentersY[i];
        const double dif = sqrt(dxif*dxif + dyif*dyif);

        const double dxij = m.cellsCentersX[i] - m.cellsCentersX[j];
        const double dyij = m.cellsCentersY[i] - m.cellsCentersY[j];
        const double dij = sqrt(dxij*dxij + dyij*dyij);

        const double geom_factor = dif / dij;

        if (i != j) {
            double f[vars];
            for (uint k=0; k<vars; ++k) {
                f[k] = (q[vars*i+k]*(1.0 - geom_factor) + q[vars*j+k] * geom_factor) * le;
            }
            for (uint k=0; k<vars; ++k) {
                gx[vars*i+k] += f[k] * nx;
                gy[vars*i+k] += f[k] * ny;

                gx[vars*j+k] -= f[k] * nx;
                gy[vars*j+k] -= f[k] * ny;
            }
        }
    }
    // normalize by cell areas
    for (uint i=0; i<m.nRealCells; ++i) {
        const double invA = 1./m.cellsAreas[i];
        for (uint k=0; k<vars; ++k) {
            gx[vars*i+k] *= invA;
            gy[vars*i+k] *= invA;
        }
    }
    // Set boundary gradients to zero
    for (uint i=m.nRealCells; i<m.cellsAreas.size(); ++i) {
        for (uint k=0; k<vars; ++k) {
            gx[vars*i+k] = 0.;
            gy[vars*i+k] = 0.;
        }
    }
}


void calc_limiters(
    std::vector<double>& limiters,
    std::vector<double>& qmin,
    std::vector<double>& qmax,
    const std::vector<double>& q,
    const std::vector<double>& gx,
    const std::vector<double>& gy,
    mesh& m
) {
    // Reset limiters to two
    for (uint i=0; i<limiters.size(); ++i) {
        limiters[i] = 1.;
    }
    // Set qmin and qmax as q
    for (uint i=0; i<q.size(); ++i) {
        qmin[i] = q[i];
        qmax[i] = q[i];
    }
    // Compute qmin and qmax
    for (uint e=0; e<m.edgesLengths.size(); ++e) {
        const auto& i = m.edgesCells(e, 0);
        const auto& j = m.edgesCells(e, 1);
        
        for (uint k=0; k<vars; ++k) {
            qmin[vars*i+k] = std::min(qmin[vars*i+k], q[vars*j+k]);
            qmin[vars*j+k] = std::min(qmin[vars*j+k], q[vars*i+k]);

            qmax[vars*i+k] = std::max(qmax[vars*i+k], q[vars*j+k]);
            qmax[vars*j+k] = std::max(qmax[vars*j+k], q[vars*i+k]);
        }
    }
    // Compute limiters
    const double tol = 1e-15;
    std::array<uint, 2> ids;
    for (uint e=0; e<m.edgesLengths.size(); ++e) {
        const auto& i = m.edgesCells(e, 0);
        const auto& j = m.edgesCells(e, 1);
        ids[0] = i;
        ids[1] = j;

        for (auto& id : ids) {
            if ((id < m.nRealCells)&(!m.cellsIsGhost[id])) {
                const double dx = m.edgesCentersX[e] - m.cellsCentersX[id];
                const double dy = m.edgesCentersY[e] - m.cellsCentersY[id];
                const double sqrt_area = sqrt(m.cellsAreas[id]);

                for (uint k=0; k<vars; ++k) {
                    double dqg = gx[vars*id+k]*dx + gy[vars*id+k]*dy;
                    
                    double delta_max = qmax[vars*id+k] - q[vars*id+k];
                    double delta_min = qmin[vars*id+k] - q[vars*id+k];

                    const double Ka = solver::limiter_k_value * sqrt_area;
                    const double K3a = Ka * Ka * Ka;
                    const double dMaxMin2 = (delta_max - delta_min)*(delta_max - delta_min); 

                    double sig;
                    if (dMaxMin2 <= K3a) {
                        sig = 1.;
                    } else if (dMaxMin2 <= 2*K3a) {
                        double y = (dMaxMin2/K3a - 1.0);
                        sig = 2.0*y*y*y - 3.0*y*y + 1.0;
                    } else {
                        sig = 0.;
                    }
                    
                    double lim = 1.0;
                    if (sig < 1.0) {
                        if (dqg > tol) {
                            lim = limiter_func(delta_max/dqg);
                        } else if (dqg < -tol) {
                            lim = limiter_func(delta_min/dqg);
                        } else {
                            lim = 1.0;
                        }
                    }

                    lim = sig + (1.0 - sig)*lim;

                    limiters[vars*id+k] = std::min(limiters[vars*id+k], lim);
                }
            }
        }
    }
}


void calc_time_derivatives(
    std::vector<double>& qt,
    const std::vector<double>& q,
    const std::vector<double>& gx,
    const std::vector<double>& gy,
    const std::vector<double>& limiters,
    mesh& m
) {
    // reset qt to be null
    for (uint i=0; i<qt.size(); ++i) {
        qt[i] = 0.;
    }
    // Compute time derivatives qt of q
    for (uint e=0; e<m.edgesNodes.cols(); ++e) {
        double n[2];
        double di[2];
        double dj[2];
        double ci[2];
        double cj[2];

        const uint i = m.edgesCells(e, 0);
        const uint j = m.edgesCells(e, 1);
        const double le = m.edgesLengths[e];
        
        n[0] = m.edgesNormalsX[e];
        n[1] = m.edgesNormalsY[e];

        const double cx = m.edgesCentersX[e];
        const double cy = m.edgesCentersY[e];

        ci[0] = m.cellsCentersX[i];
        ci[1] = m.cellsCentersY[i];

        cj[0] = m.cellsCentersX[j];
        cj[1] = m.cellsCentersY[j];

        di[0] = cx - ci[0];
        di[1] = cy - ci[1];

        dj[0] = cx - cj[0];
        dj[1] = cy - cj[1];

        // Compute edge center values
        double qi[vars];
        double qj[vars];

        if (solver::linear_interpolate) {
            for (uint k=0; k<vars; ++k) {
                const uint ki = vars*i+k;
                const uint kj = vars*j+k;
                qi[k] = q[ki] + (gx[ki]*di[0] + gy[ki]*di[1])*limiters[ki]*vars_limiters[k];
                qj[k] = q[kj] + (gx[kj]*dj[0] + gy[kj]*dj[1])*limiters[kj]*vars_limiters[k];
            }
        } else {
            for (uint k=0; k<vars; ++k) {
                qi[k] = q[vars*i+k];
                qj[k] = q[vars*j+k];
            }
        }

        // Compute viscous fluxes
        double gxv[vars];
        double gyv[vars];

        if (solver::diffusive_gradients) {
            gradient_for_diffusion(
                gxv, gyv,
                &gx[vars*i], &gy[vars*i],
                &gx[vars*j], &gy[vars*j],
                &q[vars*i], &q[vars*j],
                ci, cj
            );
        } else {
            for (uint k=0; k<vars; ++k) {
                gxv[k] = 0.;
                gyv[k] = 0.;
            }
        }

        // Compute fluxes
        double f[vars];
        calc_flux(
            f, qi, qj, gxv, gyv, n
        );

        // Update qt
        for (uint k=0; k<vars; ++k) {
            qt[vars*i+k] -= f[k] * le / m.cellsAreas[i];
            qt[vars*j+k] += f[k] * le / m.cellsAreas[j];
        }
    }
    if (solver::source_term) {
        for (uint i=0; i<m.nRealCells; ++i) {

            double source[vars];
            calc_source(source, &q[vars*i], &gx[vars*i], &gy[vars*i], m.wall_dist[i]);
            for (uint k=0; k<vars; ++k) {
                qt[vars*i+k] += source[k] * m.cellsAreas[i];
            }
        }
    }
    // if not a real cell, qt = 0
    for (uint i=0; i<m.cellsAreas.size(); ++i) {
        if (i >= m.nRealCells) {
            for (uint k=0; k<vars; ++k) {
                qt[vars*i+k] = 0.;
            }
        } else if (m.cellsIsGhost[i]) {
            for (uint k=0; k<vars; ++k) {
                qt[vars*i+k] = 0.;
            }
        }
    }
}



void update_cells(
    std::vector<double>& q,
    std::vector<double>& ql,
    const std::vector<double>& qt,
    const std::vector<double>& dt,
    const double v
) {
    for (uint i=0; i<q.size(); ++i) {
        q[i] = ql[i] + qt[i] * dt[i] * v;
    }
}

void update_bounds(
    std::vector<double>& q,
    std::vector<double>& gx,
    std::vector<double>& gy,
    std::vector<double>& limiters,
    mesh& m
) {
    // Update the ghost cells with boundary conditions
    for (uint b=0; b<m.boundaryEdges.size(); ++b) {
        const uint e = m.boundaryEdges[b];
        double n[2];
        n[0] = m.edgesNormalsX[e];
        n[1] = m.edgesNormalsY[e];
        const uint id_internal = m.edgesCells(e, 0);
        const uint id_bound = m.edgesCells(e, 1);

        double q_int[vars];
        if (solver::linear_interpolate) {
            double di[2];
            di[0] = m.edgesCentersX[e] - m.cellsCentersX[id_internal];
            di[1] = m.edgesCentersY[e] - m.cellsCentersY[id_internal];
            for (uint i=0; i<vars; ++i) {
                const uint k = vars*id_internal+i;
                q_int[i] = q[k] + (gx[k]*di[0] + gy[k]*di[1])*limiters[k]*vars_limiters[i];
            }
        } else {
            for (uint k=0; k<vars; ++k) {
                q_int[k] = q[vars*id_internal+k];
            }
        }

        m.boundaryFuncs[b](
            &q[vars*id_bound],
            q_int,
            n
        );
    }
}


void update_comms(
    std::vector<double>& q,
    mesh& m
) {

    std::vector<MPI_Request> reqs(m.comms.size());
    uint k = 0;
    for (auto& comm : m.comms) {

        uint iter = 0;
        for (const auto& i : comm.snd_indices) {
            for (uint j=0; j<vars; ++j) {
                comm.snd_q[vars*iter + j] = q[vars*i + j];
            }
            iter += 1;
        }

        // Send values
        MPI_Isend(
        /* data         = */ &comm.snd_q[0], 
        /* count        = */ comm.snd_q.size(), 
        /* datatype     = */ MPI_DOUBLE, 
        /* destination  = */ comm.out_rank, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD,
        /* request      = */ &reqs[k]
        );
        k += 1;
    }

    // Recieve values from all communicating cells
    for (auto& comm : m.comms) {
        // Recieve values
        MPI_Recv(
        /* data         = */ &comm.rec_q[0], 
        /* count        = */ comm.rec_q.size(), 
        /* datatype     = */ MPI_DOUBLE, 
        /* source       = */ comm.out_rank, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD,
        /* status       = */ MPI_STATUS_IGNORE
        );
        uint iter = 0;
        for (const auto& i : comm.rec_indices) {
            for (uint j=0; j<vars; ++j) {
                q[vars*i + j] = comm.rec_q[vars*iter + j];
            }
            iter += 1;
        }
    }

    MPI_Barrier( MPI_COMM_WORLD );

    // Free requests
    for (uint i=0; i<reqs.size(); ++i) {
        MPI_Request_free(&reqs[i]);
    }
}



void calc_residuals(
    double* R,
    std::vector<double>& qt,
    mesh& m,
    mpi_wrapper& pool
) {
    for (uint i=0; i<vars; ++i) {
        R[i] = 0.;
    }
    for (uint i=0; i<m.nRealCells; ++i) {
        if (!m.cellsIsGhost[i]) {
            for (uint j=0; j<vars; ++j) {
                R[j] += qt[vars*i+j]*qt[vars*i+j] * m.cellsAreas[i];
            }
        }
    }

    if (pool.rank != 0) {
        MPI_Send(
        /* data         = */ R, 
        /* count        = */ vars, 
        /* datatype     = */ MPI_DOUBLE, 
        /* destination  = */ 0, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD
        );
    } else {
        double R_other[vars];
        for (uint i=1; i<pool.size; ++i) {
            MPI_Recv(
            /* data         = */ R_other, 
            /* count        = */ vars, 
            /* datatype     = */ MPI_DOUBLE, 
            /* source       = */ i, 
            /* tag          = */ 0,
            /* communicator = */ MPI_COMM_WORLD,
            /* status       = */ MPI_STATUS_IGNORE
            );
            for (uint j=0; j<vars; ++j) {
                R[j] += R_other[j];
            }
        }
    }

    for (uint i=0; i<vars; ++i) {
        R[i] = sqrt(R[i]);
    }

    if (pool.rank == 0) {
        for (uint i=1; i<pool.size; ++i) {
            MPI_Send(
            /* data         = */ R, 
            /* count        = */ vars, 
            /* datatype     = */ MPI_DOUBLE, 
            /* destination  = */ i, 
            /* tag          = */ 0,
            /* communicator = */ MPI_COMM_WORLD
            );
        }
    } else {
        MPI_Recv(
        /* data         = */ R, 
        /* count        = */ vars, 
        /* datatype     = */ MPI_DOUBLE, 
        /* source       = */ 0, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD,
        /* status       = */ MPI_STATUS_IGNORE
        );
    }
    
}


void min_dt(std::vector<double>& dt, mesh& m) {
    // Minimize dt
    double min_dt = dt[0];
    for (uint i=0; i<dt.size(); ++i) {
        min_dt = std::min(min_dt, dt[i]);
    }
    for (uint i=0; i<dt.size(); ++i) {
        dt[i] = min_dt;
    }
}


void validate_dt(std::vector<double>& dt, mpi_wrapper& pool) {
    // Send dt to node 0
    if (pool.rank != 0) {
        MPI_Send(
        /* data         = */ &dt[0], 
        /* count        = */ 1, 
        /* datatype     = */ MPI_DOUBLE, 
        /* destination  = */ 0, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD
        );
    } else {
        double dti;
        for (uint i=1; i<pool.size; ++i) {
            MPI_Recv(
            /* data         = */ &dti, 
            /* count        = */ 1, 
            /* datatype     = */ MPI_DOUBLE, 
            /* source       = */ i, 
            /* tag          = */ 0,
            /* communicator = */ MPI_COMM_WORLD,
            /* status       = */ MPI_STATUS_IGNORE
            );
        }
        dt[0] = std::min(dt[0], dti);
    }
    // Send dt to all other nodes
    if (pool.rank == 0) {
        for (uint i=1; i<pool.size; ++i) {
            MPI_Send(
            /* data         = */ &dt[0], 
            /* count        = */ 1, 
            /* datatype     = */ MPI_DOUBLE, 
            /* destination  = */ i, 
            /* tag          = */ 0,
            /* communicator = */ MPI_COMM_WORLD
            );
        }
    } else {
        MPI_Recv(
        /* data         = */ &dt[0], 
        /* count        = */ 1, 
        /* datatype     = */ MPI_DOUBLE, 
        /* source       = */ 0, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD,
        /* status       = */ MPI_STATUS_IGNORE
        );

        // Dispatch this dt to all other dts
        for (uint i=1; i<dt.size(); ++i) {
            dt[i] = dt[0];
        }
    }
}



void complete_calc_qt(
    std::vector<double>& qt,
    std::vector<double>& q,
    std::vector<double>& gx,
    std::vector<double>& gy,
    std::vector<double>& qmin,
    std::vector<double>& qmax,
    std::vector<double>& limiters,
    mesh& m,
    mpi_wrapper& pool
) {
    // Compute gradients
    if (solver::do_calc_gradients) {
        calc_gradients(gx, gy, q, m);
        if (pool.size > 1) {
            update_comms(gx, m);
            update_comms(gy, m);
        }
    }

    // Compute limiters
    if (solver::do_calc_limiters) {
        calc_limiters(limiters, qmin, qmax, q, gx, gy, m);
        if (pool.size > 1) update_comms(limiters, m);
    }

    // Compute time derivative
    calc_time_derivatives(qt, q, gx, gy, limiters, m);
}



void run(
    const std::string name,
    std::vector<double>& q,
    mpi_wrapper& pool,
    mesh& m,
    solverOptions& opt
) {

    q.resize(vars*m.cellsAreas.size());
    generate_initial_solution(q, m);

    std::vector<double> qk(q.size());

    std::vector<double> qt(q.size());
    std::vector<double> gx(q.size());
    std::vector<double> gy(q.size());
    std::vector<double> limiters(q.size());
    std::vector<double> qmin(q.size());
    std::vector<double> qmax(q.size());
    std::vector<double> dt(q.size());
    std::vector<double> q_smooth0(q.size());
    std::vector<double> q_smooth1(q.size());

    // RK5 stage coefficients
    std::vector<double> alpha = {
        0.0695, 
        0.1602, 
        0.2898, 
        0.506, 
        1.
    };

    bool running = true;

    uint step = 0;
    double time = 0;

    double R0[vars];
    double R[vars];
    for (uint i=0; i<vars; ++i) {R[i] = 1.0;}

    if ((opt.verbose)&(pool.rank == 0)) {
        std::cout << "Step, Time, RealTime, ";
        for (uint i=0; i<vars; ++i) {
            std::cout << "R(q[" << i << "])";
            if (i < vars-1) {std::cout << ", ";}
        }
        std::cout << std::endl;
    }


    double save_time = opt.time_series_interval;
    uint time_step = 0;

    if (solver::do_calc_gradients) {
        calc_gradients(gx, gy, q, m);
    }

    // Init the ghost cells with boundary conditions
    update_bounds(q, gx, gy, limiters, m);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    while (running) {

        // Convergence check
        double Rmax = 0.0;
        if (step > 0) {
            for (uint i=0; i<vars; ++i) {Rmax = std::max(R[i], Rmax);}
        } else {Rmax = 1.0;}

        if ((step >= opt.max_step)|(time >= opt.max_time)|(Rmax < opt.tolerance)) {
            running = false;
            break;
        }

        // Compute time step and update comms with dt
        calc_dt(dt, q, gx, gy, m);
        if (pool.size > 1) update_comms(dt, m);
        if (solver::global_dt) {
            min_dt(dt, m);
            if (pool.size > 1) validate_dt(dt, pool);
        }
        
        // Runge kutta iterations

        // Store q in qk
        for (uint i=0; i<q.size(); ++i) qk[i] = q[i];

        for (const double& a : alpha) {
            complete_calc_qt(qt, qk, gx, gy, qmin, qmax, limiters, m, pool);
            if (solver::smooth_residuals) smooth_residuals(qt, q_smooth0, q_smooth1, m);
            update_cells(qk, q, qt, dt, a);
            update_bounds(qk, gx, gy, limiters, m);
            if (pool.size > 1) update_comms(qk, m);
        }
        // Get back qk values into q
        for (uint i=0; i<q.size(); ++i) q[i] = qk[i];

        // Compute residuals
        if (step == 0) {
            calc_residuals(R0, qt, m, pool);
            for (uint i=0; i<vars; ++i) {R[i] = R0[i];}

            if ((pool.rank == 0) & (opt.verbose)) {
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                int microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                double seconds = ((double) microseconds) / 1e6;

                std::cout << step << ", " << time << ", " << seconds << ", ";
                for (uint i=0; i<vars; ++i) {
                    std::cout << "1.0";
                    if (i < vars-1) {std::cout << ", ";}
                }
                std::cout << std::endl;
            }
        } else if ((step % opt.print_interval == 0)|(opt.tolerance > 1.01e-16)) {
            
            calc_residuals(R, qt, m, pool);
            for (uint i=0; i<vars; ++i) {R[i] = R[i]/R0[i];}

            if ((step % opt.print_interval == 0) & (opt.verbose) & (pool.rank == 0)) {
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                int microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
                double seconds = ((double) microseconds) / 1e6;

                std::cout << step << ", " << time << ", " << seconds << ", ";
                for (uint i=0; i<vars; ++i) {
                    std::cout << R[i];
                    if (i < vars-1) {std::cout << ", ";}
                }
                std::cout << std::endl;
            }
        }

        // If save time series, save time series
        if (opt.save_time_series) {
            if (time > save_time) {
                // Save file to ./times/ folder
                std::string timename = std::to_string(time_step);
                writeVtk(name, q, m, pool.rank, pool.size, timename);
                save_time += opt.time_series_interval;
                time_step += 1;
            }
        }

        // Edit step and time
        step += 1;
        time += dt[0];
    }

    if ((opt.verbose)&(pool.rank == 0)) {
        // End prints
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        int microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        double seconds = ((double) microseconds) / 1e6;

        std::cout << step << ", " << time << ", " << seconds << ", ";
        for (uint i=0; i<vars; ++i) {
            std::cout << R[i];
            if (i < vars-1) {std::cout << ", ";}
        }
        std::cout << std::endl;
    }

}


}