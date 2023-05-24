#include "fem_2.h"

int main() {
   CGM matrix;
   matrix.read_data();
   matrix.read_boundary();
   matrix.making_grid();
   matrix.read_time_nonform();
   //matrix.making_time_grid();
   matrix.init_cond();
   matrix.glob_G();
   matrix.glob_M();
   //matrix.time_scheme_A();
   for (int i = 2; i < matrix.n_times; i++) {
      //matrix.time_scheme_d(i);
      matrix.time_scheme_nonform(i);
      matrix.CGM_precond_ILU();
      matrix.print_result(i);
   }
   return 0;
}
