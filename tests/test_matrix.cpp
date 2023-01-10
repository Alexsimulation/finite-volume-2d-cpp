#include <fvhyper/matrix.h>
#include <string>

/*
Test of the sparse matrix-vector multiplication
Multiply :
    m = (
        1.9   1.8     0     0
          0   1.6     0     0
          0     0   1.5   1.2
          0     0   1.1     0
    )
with vector :
    v = (1, 0, 3, 4)

Compile with : mpic++ test_matrix.cpp ../fvhyper/src/parallel.cpp -o test_matrix -I../
*/


int main() {
    fvhyper::mpi_wrapper pool;


    fvhyper::matrix<double> m;
    m.set_pool(&pool);

    if (pool.rank == 0) {
        m.place_new(1.9, 0, 0);
        m.place_new(1.8, 0, 1);
        m.place_new(1.6, 1, 1);
    } else if (pool.rank == 1) {
        m.place_new(1.5, 2, 2);
        m.place_new(1.2, 2, 3);
        m.place_new(1.1, 3, 2);
    }

    fvhyper::vector<double> v;
    if (pool.rank == 0) {
        v.set(0) = 1;
    } else if (pool.rank == 1) {
        v.set(2) = 3;
        v.set(3) = 4;
    }

    fvhyper::vector<double> ans;
    if (pool.rank == 0) {
        ans.set(0) = 0;
        ans.set(1) = 0;
    } else if (pool.rank == 1) {
        ans.set(2) = 0;
        ans.set(3) = 0;
    }

    m.mult(ans, v);

    ans.print();

    return pool.exit();
}
