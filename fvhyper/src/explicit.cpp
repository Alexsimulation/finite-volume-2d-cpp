#include <fvhyper/explicit.h>
#include <array>



namespace fvhyper {

void calc_gradients(
    std::vector<double>& g,
    const std::vector<double>& q,
    mesh& m
) {
    // reset gradients to be null
    for (uint i=0; i<g.size(); ++i) {
        g[i] = 0.;
    }
    double f[vars];
    // Update gradients using green gauss cell based
    for (uint e=0; e<m.edgesNodes.cols(); ++e) {
        const auto& i = m.edgesCells(e, 0);
        const auto& j = m.edgesCells(e, 1);
        const auto& nx = m.edgesNormalsX[e];
        const auto& ny = m.edgesNormalsY[e];
        const auto& le = m.edgesLengths[e];

        for (uint k=0; k<vars; ++k) {
            f[k] = (q[vars*i+k] + q[vars*j+k]) * 0.5 * le;
        }
        for (uint k=0; k<vars; ++k) {
            g[2*vars*i+2*k] += f[k] * nx;
            g[2*vars*i+2*k+1] += f[k] * ny;

            g[2*vars*j+2*k] -= f[k] * nx;
            g[2*vars*j+2*k+1] -= f[k] * ny;
        }
    }
    // normalize by cell areas
    for (uint i=0; i<m.cellsAreas.size(); ++i) {
        const double invA = 1./m.cellsAreas[i];
        for (uint k=0; k<vars; ++k) {
            g[2*vars*i+2*k] *= invA;
            g[2*vars*i+2*k+1] *= invA;
        }
    }
}


void calc_limiters(
    std::vector<double>& limiters,
    std::vector<double>& qmin,
    std::vector<double>& qmax,
    const std::vector<double>& q,
    const std::vector<double>& g,
    mesh& m
) {
    // Reset limiters to two
    for (uint i=0; i<limiters.size(); ++i) {
        limiters[i] = 2.;
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
    // Compute limiters r variable
    const double tol = 1e-15;
    std::array<uint, 2> ids;
    for (uint e=0; e<m.edgesLengths.size(); ++e) {
        const auto& i = m.edgesCells(e, 0);
        const auto& j = m.edgesCells(e, 1);
        const double dx_i = m.edgesCentersX[e] - m.cellsCentersX[i];
        const double dy_i = m.edgesCentersY[e] - m.cellsCentersY[i];
        const double dx_j = m.edgesCentersX[e] - m.cellsCentersX[j];
        const double dy_j = m.edgesCentersY[e] - m.cellsCentersY[j];
        ids[0] = i;
        ids[1] = j;
        for (auto& id : ids) {
            for (uint k=0; k<vars; ++k) {
                double dqg = g[2*vars*id+2*k]*dx_i + g[2*vars*id+2*k+1]*dy_i;
                if (dqg > std::max(qmax[vars*i+k] - q[vars*i+k], tol)) {
                    limiters[vars*id+k] = std::min(limiters[vars*id+k], (qmax[vars*id+k] - q[vars*id+k])/dqg);
                } else if (dqg < std::min(qmin[vars*i+k] - q[vars*id+k], -tol)) {
                    limiters[vars*id+k] = std::min(limiters[vars*id+k], (qmin[vars*id+k] - q[vars*id+k])/dqg);
                } else {
                    limiters[vars*id+k] = std::min(limiters[vars*id+k], 1.);
                }
            }
        }
    }
    // Compute limiters
    for (uint i=0; i<m.cellsAreas.size(); ++i) {
        double r[vars];
        for (uint j=0; j<vars; ++j) {
            r[j] = limiters[vars*i + j];
        }
        limiter_func(&limiters[i], r);
    }
}


void calc_time_derivatives(
    std::vector<double>& qt,
    const std::vector<double>& q,
    const std::vector<double>& g,
    const std::vector<double>& limiters,
    mesh& m
) {
    // reset qt to be null
    for (uint i=0; i<qt.size(); ++i) {
        qt[i] = 0.;
    }
    // Compute time derivatives qt of q
    double f[vars];
    double n[2];
    double di[2];
    double dj[2];
    for (uint e=0; e<m.edgesNodes.cols(); ++e) {
        const auto& i = m.edgesCells(e, 0);
        const auto& j = m.edgesCells(e, 1);
        const auto& le = m.edgesLengths[e];
        
        n[0] = m.edgesNormalsX[e];
        n[1] = m.edgesNormalsY[e];

        const double cx = (m.nodesX[m.edgesNodes(e, 0)] + m.nodesX[m.edgesNodes(e, 1)])*0.5;
        const double cy = (m.nodesY[m.edgesNodes(e, 0)] + m.nodesY[m.edgesNodes(e, 1)])*0.5;

        di[0] = cx - m.cellsCentersX[i];
        di[1] = cy - m.cellsCentersY[i];

        dj[0] = cx - m.cellsCentersX[j];
        dj[1] = cy - m.cellsCentersY[j];

        calc_flux(
            f, &q[vars*i], &q[vars*j], 
            &g[2*vars*i], &g[2*vars*j], &limiters[i],
            n, di, dj, m.cellsAreas[i], m.edgesLengths[e]
        );
        for (uint k=0; k<vars; ++k) {
            qt[vars*i+k] -= f[k] * le;
            qt[vars*j+k] += f[k] * le;
        }
    }
    // normalize by cell areas, and if not a real cell, qt = 0
    for (uint i=0; i<m.cellsAreas.size(); ++i) {
        if (i >= m.nRealCells) {
            for (uint k=0; k<vars; ++k) {
                qt[vars*i+k] = 0.;
            }
        } else if (m.cellsIsGhost[i]) {
            for (uint k=0; k<vars; ++k) {
                qt[vars*i+k] = 0.;
            }
        } else {
            const double invA = 1./m.cellsAreas[i];
            for (uint k=0; k<vars; ++k) {
                qt[vars*i+k] *= invA;
            }
        }
    }
}



void update_cells(
    std::vector<double>& q,
    const std::vector<double>& qt,
    const double dt
) {
    for (uint i=0; i<q.size(); ++i) {
        q[i] += qt[i] * dt;
    }
}


void update_comms(
    std::vector<double>& q,
    mesh& m
) {
    for (auto& comm : m.comms) {
        uint iter = 0;
        for (const auto& i : comm.snd_indices) {
            for (uint j=0; j<vars; ++j) {
                comm.snd_q[vars*iter + j] = q[vars*i + j];
            }
            iter += 1;
        }
        // Send values
        MPI_Send(
        /* data         = */ &comm.snd_q[0], 
        /* count        = */ comm.snd_q.size(), 
        /* datatype     = */ MPI_DOUBLE, 
        /* destination  = */ comm.out_rank, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD
        );
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
                R[j] += (qt[vars*i+j]*qt[vars*i+j]) * m.cellsAreas[i];
            }
        }
    }
    for (uint i=0; i<vars; ++i) {
        R[i] = sqrt(R[i]);
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
    
}


void validate_dt(double& dt, mpi_wrapper& pool) {
    // Send dt to node 0
    if (pool.rank != 0) {
        MPI_Send(
        /* data         = */ &dt, 
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
        dt = std::min(dt, dti);
    }
    // Send dt to all other nodes
    if (pool.rank == 0) {
        for (uint i=1; i<pool.size; ++i) {
            MPI_Send(
            /* data         = */ &dt, 
            /* count        = */ 1, 
            /* datatype     = */ MPI_DOUBLE, 
            /* destination  = */ i, 
            /* tag          = */ 0,
            /* communicator = */ MPI_COMM_WORLD
            );
        }
    } else {
        MPI_Recv(
        /* data         = */ &dt, 
        /* count        = */ 1, 
        /* datatype     = */ MPI_DOUBLE, 
        /* source       = */ 0, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD,
        /* status       = */ MPI_STATUS_IGNORE
        );
    }
}



void run(
    std::vector<double>& q,
    mpi_wrapper& pool,
    mesh& m,
    solverOptions& opt
) {

    q.resize(vars*m.cellsAreas.size());
    generate_initial_solution(q, m);

    std::vector<double> qt(q.size());
    std::vector<double> g(2*q.size());
    std::vector<double> limiters(q.size());
    std::vector<double> qmin(q.size());
    std::vector<double> qmax(q.size());

    bool running = true;
    double dt = 0;

    uint step = 0;
    double time = 0;

    double R0[vars];
    while (running) {

        if ((step >= opt.max_step)|(time >= opt.max_time)) {
            running = false;
            break;
        }

        // Update the ghost cells with boundary conditions
        for (uint i=0; i<m.boundaryEdges.size(); ++i) {
            double n[2];
            uint e = m.boundaryEdges[i];
            n[0] = m.edgesNormalsX[e];
            n[1] = m.edgesNormalsY[e];
            m.boundaryFuncs[i](
                &q[vars*m.edgesCells(e, 1)],
                &q[vars*m.edgesCells(e, 0)],
                n
            );
        }

        // Compute stuff
        calc_dt(dt, q, m);
        if (pool.size > 1) validate_dt(dt, pool);

        // Compute gradients
        //calc_gradients(g, q, m);

        // Compute limiters
        //calc_limiters(limiters, qmin, qmax, q, g, m);

        // Compute time derivative
        calc_time_derivatives(qt, q, g, limiters, m);

        // Euler explicit update cells
        update_cells(q, qt, dt);

        // Link with MPI
        // Compute update to send to neighboor nodes

        if (pool.size > 1) update_comms(q, m);

        if (step == 0) {
            calc_residuals(R0, qt, m, pool);
        } else if (step % opt.print_interval == 0) {
            
            double R[vars];
            calc_residuals(R, qt, m, pool);

            if (pool.rank == 0) {
                std::cout << "Step = " << step << ", ";
                std::cout << "Time = " << time << "\n";
                for (uint i=0; i<vars; ++i) {
                    std::cout << " - R(q[" << i << "]) = " << R[i]/R0[i] << "\n";
                }
                std::cout << std::endl;
            }
        }

        // Edit stop criteria
        step += 1;
        time += dt;
    }

}


}