#include "fem_2.h"

int main() {
   CGM matrix;
   matrix.read_data();
   matrix.read_boundary();
   matrix.making_grid();
   matrix.making_time_grid();
   matrix.glob_G();
   matrix.glob_M();
   for (int i = 2; i < matrix.n_times; i++) {
      double dt = matrix.times[i] - matrix.times[i - 1];
      matrix.glob_F(matrix.times[i - 1]);
      matrix.time_scheme(dt);
      matrix.first_cond(matrix.times[i]);
      matrix.CGM_precond_ILU();
      matrix.q_init_0 = matrix.q_init_1;
      matrix.q_init_1 = matrix.q;
      matrix.print_result(matrix.times[i]);
   }

   return 0;
}
