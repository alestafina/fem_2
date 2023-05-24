#include "fem_2.h"

double Functions::func(double x, double y, double t) {
   return 0.1 * (x + y);
}

double Functions::u_g(double x, double y, double t) {
   return (x + y) * t;
}

double Functions::u_betta(double x, double y, double t) {
   return 1;
}

double Functions::theta(double x, double y, double t) {
   return 1;
}

double Functions::u_real(double x, double y, double t) {
   return (x + y) * t;
}

FEM::FEM() {
   A = matrix();
   Msigma = matrix();
   Mchi = matrix();
   G = matrix();
   loc_M = vector<vector<double>>();
   loc_G = vector<vector<double>>();
   loc_f = vector<double>();
   F = vector<double>();
   q = vector<double>();
   q_real = vector<double>();

   d = vector<double>();
   times = vector<double>();
   q_init_2 = vector<double>();
   q_init_1 = vector<double>();

   x_nodes = vector<double>();
   y_nodes = vector<double>();
   nodes = vector<int>();
   W = vector<vector<int>>();

   first_bc = vector<int>();
   second_bc = vector<int>();
   thrid_bc = vector<int>();

   hx = 0, hy = 0, ht = 0;
   N = 0, M = 0, nx = 0, ny = 0;
   betta = 0, lambda = 0, chi = 0, sigma = 0;
   n_times = 0;
}

FEM::~FEM() {
   A.~matrix();
   loc_M.~vector();
   loc_G.~vector();
   loc_f.~vector();
   F.~vector();
   q.~vector();
   q_real.~vector();
   W.~vector();
   x_nodes.~vector();
   y_nodes.~vector();
   nodes.~vector();
   first_bc.~vector();
   second_bc.~vector();
   thrid_bc.~vector();
   d.~vector();
   times.~vector();
   q_init_2.~vector();
   q_init_1.~vector();
}

void FEM::read_data() {
   ifstream grid("grid.txt");

   grid >> lambda >> betta >> chi >> sigma;

   grid >> nx;
   x_nodes.resize(nx);
   grid >> x_nodes[0] >> x_nodes[nx - 1];

   grid >> ny;
   y_nodes.resize(ny);
   grid >> y_nodes[0] >> y_nodes[ny - 1];

   grid >> n_times;
   times.resize(n_times);
   grid >> times[0] >> times[n_times - 1];

   grid.close();
}

void FEM::making_grid() {
   M = nx * ny;
   N = (nx - 1) * (ny - 1);

   hx = (x_nodes[nx - 1] - x_nodes[0]) / (nx - 1.0);
   hy = (y_nodes[ny - 1] - y_nodes[0]) / (ny - 1.0);

   for (int i = 1; i < nx; i++) x_nodes[i] = x_nodes[i - 1] + hx;
   for (int i = 1; i < ny; i++) y_nodes[i] = y_nodes[i - 1] + hy;

   W.resize(N);
   for (int i = 0; i < N; i++) W[i].resize(4);

   for (int i = 0; i < ny - 1; i++) {
      for (int j = 0; j < nx - 1; j++) {
         W[i * (nx - 1) + j][0] = j;
         W[i * (nx - 1) + j][1] = j + 1;
         W[i * (nx - 1) + j][2] = i;
         W[i * (nx - 1) + j][3] = i + 1;
      }
   }

   A.di.resize(M);
   A.ig.resize(M + 1);

   Msigma.di.resize(M);
   Msigma.ig.resize(M + 1);

   Mchi.di.resize(M);
   Mchi.ig.resize(M + 1);

   G.di.resize(M);
   G.ig.resize(M + 1);

   loc_G.resize(4);
   for (int i = 0; i < 4; i++) loc_G[i].resize(4);

   loc_M.resize(4);
   for (int i = 0; i < 4; i++) loc_M[i].resize(4);

   loc_f.resize(4);
   F.resize(M);

   d.resize(M);
   q.resize(M);
}

void FEM::making_time_grid() {
   ht = (times[n_times - 1] - times[0]) / (n_times - 1.0);
   for (int i = 1; i < n_times; i++) times[i] = times[i - 1] + ht;
}

void FEM::init_cond() {
   Functions init;
   q_init_2.resize(M);
   q_init_1.resize(M);
   for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
         q_init_2[i + j * nx] = init.u_real(x_nodes[i], y_nodes[j], times[0]);
         q_init_1[i + j * nx] = init.u_real(x_nodes[i], y_nodes[j], times[1]);
      }
   }
}

void FEM::read_time_nonform() {
   ifstream gridT("time.txt");
   times.clear();
   n_times = 0;

   gridT >> n_times;

   times.resize(n_times);
   for (int i = 0; i < n_times; i++) gridT >> times[i];
}

void FEM::build_portrait(matrix &A) {
   int n = 0;
   vector<int> tmp(4);

   for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
         n = nx * j + i;
         if (n > 1) {
            if (n <= nx && n > 1) A.ig[n] = A.ig[n - 1] + 1;
            else if (n == j * nx + 1) A.ig[n] = A.ig[n - 1] + 2;
            else if (n == j * nx) A.ig[n] = A.ig[n - 1] + 3;
            else A.ig[n] = A.ig[n - 1] + 4;
         }
      }
   }
   A.ig[n + 1] = A.ig[n] + 3;

   for (int j = 1; j < M; j++) {
      tmp[0] = j - nx - 1;
      tmp[1] = j - nx;
      tmp[2] = j - nx + 1;
      tmp[3] = j - 1;
      for (int i = 0; i < A.ig[j + 1] - A.ig[j]; i++) {
         if (j < nx) A.jg.push_back(tmp[3]);
         else if (j % nx == 0) A.jg.push_back(tmp[i + 1]);
         else if (j % nx + 1 == nx) A.jg.push_back(tmp[i == 2 ? i + 1 : i]);
         else A.jg.push_back(tmp[i]);
      }
   }

   A.au.resize(A.ig.back());
   A.al.resize(A.ig.back());
}

void FEM::glob_M() {
   vector<double> g_nodes(4);
   double el1, el2, el3, el4;

   build_portrait(Msigma);

   for (int i = 0; i < N; i++) {
      el1 = sigma *hx * hy / 9.0;
      el2 = el3 = sigma * hx * hy / 18.0;
      el4 = sigma * hx * hy / 36.0;

      loc_M[0][0] = el1, loc_M[0][1] = el2, loc_M[0][2] = el3, loc_M[0][3] = el4;
      loc_M[1][0] = el2, loc_M[1][1] = el1, loc_M[1][2] = el4, loc_M[1][3] = el3;
      loc_M[2][0] = el3, loc_M[2][1] = el4, loc_M[2][2] = el1, loc_M[2][3] = el2;
      loc_M[3][0] = el4, loc_M[3][1] = el3, loc_M[3][2] = el2, loc_M[3][3] = el1;

      g_nodes[0] = nx * W[i][2] + W[i][0];
      g_nodes[1] = nx * W[i][2] + W[i][1];
      g_nodes[2] = nx * W[i][3] + W[i][0];
      g_nodes[3] = nx * W[i][3] + W[i][1];
      
      for (int j = 0; j < 4; j++) {
         for (int l = 0; l < 4; l++) {
            add_elem(g_nodes[j], g_nodes[l], loc_M[j][l], Msigma);
         }
      }
   }
   
   build_portrait(Mchi);

   for (int i = 0; i < N; i++) {
      el1 = chi * hx * hy / 9.0;
      el2 = el3 = chi * hx * hy / 18.0;
      el4 = chi * hx * hy / 36.0;

      loc_M[0][0] = el1, loc_M[0][1] = el2, loc_M[0][2] = el3, loc_M[0][3] = el4;
      loc_M[1][0] = el2, loc_M[1][1] = el1, loc_M[1][2] = el4, loc_M[1][3] = el3;
      loc_M[2][0] = el3, loc_M[2][1] = el4, loc_M[2][2] = el1, loc_M[2][3] = el2;
      loc_M[3][0] = el4, loc_M[3][1] = el3, loc_M[3][2] = el2, loc_M[3][3] = el1;

      g_nodes[0] = nx * W[i][2] + W[i][0];
      g_nodes[1] = nx * W[i][2] + W[i][1];
      g_nodes[2] = nx * W[i][3] + W[i][0];
      g_nodes[3] = nx * W[i][3] + W[i][1];

      for (int j = 0; j < 4; j++) {
         for (int l = 0; l < 4; l++) {
            add_elem(g_nodes[j], g_nodes[l], loc_M[j][l], Mchi);
         }
      }
   }
}

void FEM::glob_G() {
   vector<double> g_nodes(4);
   double el1, el2, el3, el4;

   build_portrait(G);
   build_portrait(A);

   for (int i = 0; i < N; i++) {
      el1 = lambda / 6.0 * 2.0 * (hy / hx + hx / hy);
      el2 = lambda / 6.0 * ((-2.0 * hy) / hx + hx / hy);
      el3 = lambda / 6.0 * (hy / hx + (-2.0 * hx) / hy);
      el4 = -lambda / 6.0 * (hy / hx + hx / hy);

      loc_G[0][0] = el1, loc_G[0][1] = el2, loc_G[0][2] = el3, loc_G[0][3] = el4;
      loc_G[1][0] = el2, loc_G[1][1] = el1, loc_G[1][2] = el4, loc_G[1][3] = el3;
      loc_G[2][0] = el3, loc_G[2][1] = el4, loc_G[2][2] = el1, loc_G[2][3] = el2;
      loc_G[3][0] = el4, loc_G[3][1] = el3, loc_G[3][2] = el2, loc_G[3][3] = el1;

      g_nodes[0] = nx * W[i][2] + W[i][0];
      g_nodes[1] = nx * W[i][2] + W[i][1];
      g_nodes[2] = nx * W[i][3] + W[i][0];
      g_nodes[3] = nx * W[i][3] + W[i][1];

      for (int j = 0; j < 4; j++) {
         for (int l = 0; l < 4; l++) {
            add_elem(g_nodes[j], g_nodes[l], loc_G[j][l], G);
         }
      }
   }
}

void FEM::glob_F(double t) {
   vector<double> g_nodes(4);
   double el1, el2, el3, el4;

   Functions func;

   for (int i = 0; i < N; i++) {
      el1 = func.func(x_nodes[W[i][0]], y_nodes[W[i][2]], t);
      el2 = func.func(x_nodes[W[i][1]], y_nodes[W[i][2]], t);
      el3 = func.func(x_nodes[W[i][0]], y_nodes[W[i][3]], t);
      el4 = func.func(x_nodes[W[i][1]], y_nodes[W[i][3]], t);

      loc_f[0] = hx * hy / 36.0 * (4.0 * el1 + 2.0 * el2 + 2.0 * el3 + el4);
      loc_f[1] = hx * hy / 36.0 * (2.0 * el1 + 4.0 * el2 + el3 + 2.0 * el4);
      loc_f[2] = hx * hy / 36.0 * (2.0 * el1 + el2 + 4.0 * el3 + 2.0 * el4);
      loc_f[3] = hx * hy / 36.0 * (el1 + 2.0 * el2 + 2.0 * el3 + 4.0 * el4);

      g_nodes[0] = nx * W[i][2] + W[i][0];
      g_nodes[1] = nx * W[i][2] + W[i][1];
      g_nodes[2] = nx * W[i][3] + W[i][0];
      g_nodes[3] = nx * W[i][3] + W[i][1];

      for (int j = 0; j < 4; j++) {
         F[g_nodes[j]] += loc_f[j];
      }
   }
}

vector<double> FEM::matr_vec_mult(vector<double> &x, matrix &A) {
   vector<double> result(x.size()), L = A.al, U = A.au;
   for (int i = 0; i < x.size(); ++i) {
      result[i] = A.di[i] * x[i];
      for (int j = A.ig[i]; j < A.ig[i + 1]; j++) {
         result[i] += L[j] * x[A.jg[j]];
         result[A.jg[j]] += U[j] * x[i];
      }
   }
   return result;
}

void FEM::add_elem(int i, int j, double elem, matrix &A) {
   int k = 0;
   if (i == j) {
      A.di[i] += elem;
   }
   else {
      if (i > j) {
         for (k = A.ig[i]; k < A.ig[i + 1]; k++) {
            if (A.jg[k] == j) A.al[k] += elem;
         }
      }
      else {
         for (k = A.ig[j]; k < A.ig[j + 1]; k++) {
            if (A.jg[k] == i) A.au[k] += elem;
         }
      }
   }
}

void FEM::time_scheme_A() {
   // считаем матрицу —Ћј” на i-том слое
   for (int j = 0; j < Mchi.al.size(); j++) {
      A.al[j] = (Mchi.al[j] / (ht * ht)) + (Msigma.al[j] / (2.0 * ht));
      A.au[j] = (Mchi.au[j] / (ht * ht)) + (Msigma.au[j] / (2.0 * ht));
   }
   for (int j = 0; j < Mchi.di.size(); j++) {
      A.di[j] = (Mchi.di[j] / (ht * ht)) + (Msigma.di[j] / (2.0 * ht));
   }
}

void FEM::time_scheme_d(int i) {
   matrix tmpMatrix = Mchi;
   vector<double> tmpVect(M);

   d.clear();
   F.clear();
   d.resize(M);
   F.resize(M);

   // считаем вектор правой части на i-том слое
   glob_F(times[i - 1]);
   d = F;

   for (int j = 0; j < Mchi.al.size(); j++) {
      tmpMatrix.al[j] = (2.0 * Mchi.al[j]) / (ht * ht);
      tmpMatrix.au[j] = (2.0 * Mchi.au[j]) / (ht * ht);
   }
   for (int j = 0; j < Mchi.di.size(); j++) {
      tmpMatrix.di[j] = (2.0 * Mchi.di[j]) / (ht * ht);
   }
   tmpVect = matr_vec_mult(q_init_1, tmpMatrix);
   for (int j = 0; j < M; j++) d[j] += tmpVect[j];


   for (int j = 0; j < Mchi.al.size(); j++) {
      tmpMatrix.al[j] = (Mchi.al[j]) / (ht * ht);
      tmpMatrix.au[j] = (Mchi.au[j]) / (ht * ht);
   }
   for (int j = 0; j < Mchi.di.size(); j++) {
      tmpMatrix.di[j] = (Mchi.di[j]) / (ht * ht);
   }
   tmpVect = matr_vec_mult(q_init_2, tmpMatrix);
   for (int j = 0; j < M; j++) d[j] -= tmpVect[j];


   for (int j = 0; j < Msigma.al.size(); j++) {
      tmpMatrix.al[j] = (Msigma.al[j]) / (2.0 * ht);
      tmpMatrix.au[j] = (Msigma.au[j]) / (2.0 * ht);
   }
   for (int j = 0; j < Msigma.di.size(); j++) {
      tmpMatrix.di[j] = (Msigma.di[j]) / (2.0 * ht);
   }
   tmpVect = matr_vec_mult(q_init_2, tmpMatrix);
   for (int j = 0; j < M; j++) d[j] += tmpVect[j];


   tmpVect = matr_vec_mult(q_init_1, G);
   for (int j = 0; j < M; j++) d[j] -= tmpVect[j];

   first_cond(times[i]);
}

void FEM::time_scheme_nonform(int i) {
   double dt = times[i] - times[i - 2];
   double dt0 = times[i] - times[i - 1];
   double dt1 = times[i - 1] - times[i - 2];

   matrix tmpMatrix = Mchi;
   vector<double> tmpVec(M);

   A.di.clear();
   A.al.clear();
   A.au.clear();
   d.clear();
   F.clear();
   A.di.resize(M);
   A.al.resize(A.ig.back());
   A.au.resize(A.ig.back());
   d.resize(M);
   F.resize(M);

   // matrix A
   for (int k = 0; k < A.al.size(); k++) {
      A.al[k] = ((2.0 * Mchi.al[k]) / (dt * dt0)) + ((dt1 * Msigma.al[k]) / (dt * dt0));
      A.au[k] = ((2.0 * Mchi.au[k]) / (dt * dt0)) + ((dt1 * Msigma.au[k]) / (dt * dt0));
   }
   for (int k = 0; k < A.di.size(); k++) {
      A.di[k] = ((2.0 * Mchi.di[k]) / (dt * dt0)) + ((dt1 * Msigma.di[k]) / (dt * dt0));
   }

   // vector d
   glob_F(times[i - 1]);
   d = F;
   for (int k = 0; k < Mchi.al.size(); k++) {
      tmpMatrix.al[k] = (2.0 * Mchi.al[k]) / (dt * dt1);
      tmpMatrix.au[k] = (2.0 * Mchi.au[k]) / (dt * dt1);
   }
   for (int k = 0; k < Mchi.di.size(); k++) {
      tmpMatrix.di[k] = (2.0 * Mchi.di[k]) / (dt * dt1);
   }
   tmpVec = matr_vec_mult(q_init_2, tmpMatrix);
   for (int k = 0; k < M; k++) d[k] -= tmpVec[k];


   for (int k = 0; k < Mchi.al.size(); k++) {
      tmpMatrix.al[k] = (2.0 * Mchi.al[k]) / (dt0 * dt1);
      tmpMatrix.au[k] = (2.0 * Mchi.au[k]) / (dt0 * dt1);
   }
   for (int k = 0; k < Mchi.di.size(); k++) {
      tmpMatrix.di[k] = (2.0 * Mchi.di[k]) / (dt0 * dt1);
   }
   tmpVec = matr_vec_mult(q_init_1, tmpMatrix);
   for (int k = 0; k < M; k++) d[k] += tmpVec[k];


   for (int k = 0; k < Msigma.al.size(); k++) {
      tmpMatrix.al[k] = (dt0 * Msigma.al[k]) / (dt * dt1);
      tmpMatrix.au[k] = (dt0 * Msigma.au[k]) / (dt * dt1);
   }
   for (int k = 0; k < Mchi.di.size(); k++) {
      tmpMatrix.di[k] = (dt0 * Msigma.di[k]) / (dt * dt1);
   }
   tmpVec = matr_vec_mult(q_init_2, tmpMatrix);
   for (int k = 0; k < M; k++) d[k] += tmpVec[k];


   for (int k = 0; k < Msigma.al.size(); k++) {
      tmpMatrix.al[k] = ((dt0 - dt1) * Msigma.al[k]) / (dt0 * dt1);
      tmpMatrix.au[k] = ((dt0 - dt1) * Msigma.au[k]) / (dt0 * dt1);
   }
   for (int k = 0; k < Mchi.di.size(); k++) {
      tmpMatrix.di[k] = ((dt0 - dt1) * Msigma.di[k]) / (dt0 * dt1);
   }
   tmpVec = matr_vec_mult(q_init_1, tmpMatrix);
   for (int k = 0; k < M; k++) d[k] += tmpVec[k];


   tmpVec = matr_vec_mult(q_init_1, G);
   for (int k = 0; k < M; k++) d[k] += tmpVec[k];

   first_cond(times[i]);
}

void FEM::read_boundary() {
   int num, count;
   ifstream bc("boundary.txt");

   bc >> num >> count;
   if (num == 1) {
      first_bc.resize(count);
      for (int i = 0; i < count; i++) {
         bc >> first_bc[i];
      }
      bc >> num >> count;
   }

   if (num == 2) {
      second_bc.resize(2 * count);
      for (int i = 0; i < 2 * count; i++) {
         bc >> second_bc[i];
      }
      bc >> num >> count;
   }

   if (num == 3) {
      thrid_bc.resize(2 * count);
      for (int i = 0; i < 2 * count; i++) {
         bc >> thrid_bc[i];
      }
   }

   bc.close();
}

void FEM::thrid_cond(double t) {
   Functions func;
   vector<vector<double>> loc_AS3;
   vector<double> loc_fS3(2);
   double el1A = 0, el2A = 0, el1f = 0, el2f = 0;
   double h = 0;

   loc_AS3.resize(2, vector<double>(2));

   for (int k = 0; k < thrid_bc.size() - 1; k += 2) {
      int i = thrid_bc[k] % nx;
      int j = thrid_bc[k] / nx;
      if (i == thrid_bc[k + 1] % nx) { // граница вертикальна€
         h = y_nodes[j + 1] - y_nodes[j];

         el1A = 2.0 * betta * h / 6.0;
         el2A = betta * h / 6.0;

         el1f = betta * h * (2.0 * func.u_betta(x_nodes[i], y_nodes[j], t) + func.u_betta(x_nodes[i], y_nodes[j + 1], t)) / 6.0;
         el2f = betta * h * (func.u_betta(x_nodes[i], y_nodes[j], t) + 2.0 * func.u_betta(x_nodes[i], y_nodes[j + 1], t)) / 6.0;

      }
      else if (j == thrid_bc[k + 1] / nx) { // граница горизонтальна€ 
         h = x_nodes[i + 1] - x_nodes[i];

         el1A = 2.0 * betta * h / 6.0;
         el2A = betta * h / 6.0;

         el1f = betta * h * (2.0 * func.u_betta(x_nodes[i], y_nodes[j], t) + func.u_betta(x_nodes[i + 1], y_nodes[j], t)) / 6.0;
         el2f = betta * h * (func.u_betta(x_nodes[i], y_nodes[j], t) + 2.0 * func.u_betta(x_nodes[i + 1], y_nodes[j], t)) / 6.0;
      }

      loc_AS3[0][0] = loc_AS3[1][1] = el1A;
      loc_AS3[0][1] = loc_AS3[1][0] = el2A;

      loc_fS3[0] = el1f;
      loc_fS3[1] = el2f;

      d[thrid_bc[k]] += loc_fS3[0];
      d[thrid_bc[k + 1]] += loc_fS3[1];

      add_elem(thrid_bc[k], thrid_bc[k], loc_AS3[0][0], A);
      add_elem(thrid_bc[k], thrid_bc[k + 1], loc_AS3[0][1], A);
      add_elem(thrid_bc[k + 1], thrid_bc[k], loc_AS3[1][0], A);
      add_elem(thrid_bc[k + 1], thrid_bc[k + 1], loc_AS3[1][1], A);
   }
}

void FEM::second_cond(double t) {
   Functions func;
   vector<double> loc_fS2(2);
   double el1 = 0, el2 = 0;
   double h = 0;

   for (int k = 0; k < second_bc.size(); k += 2) {
      int i = second_bc[k] % nx;
      int j = second_bc[k] / nx;
      if (i == second_bc[k + 1] % nx) { // граница вертикальна€
         h = y_nodes[j + 1] - y_nodes[j];

         el1 = h * (2.0 * func.theta(x_nodes[i], y_nodes[j], t) + func.theta(x_nodes[i], y_nodes[j + 1], t)) / 6.0;
         el2 = h * (func.theta(x_nodes[i], y_nodes[j], t) + 2.0 * func.theta(x_nodes[i], y_nodes[j + 1], t)) / 6.0;
      }
      else if (j == second_bc[k + 1] / nx) { // граница горизонтальна€ 
         h = x_nodes[i + 1] - x_nodes[i];

         el1 = h * (2.0 * func.theta(x_nodes[i], y_nodes[j], t) + func.theta(x_nodes[i + 1], y_nodes[j], t)) / 6.0;
         el2 = h * (func.theta(x_nodes[i], y_nodes[j], t) + 2.0 * func.theta(x_nodes[i + 1], y_nodes[j], t)) / 6.0;
      }

      loc_fS2[0] = el1;
      loc_fS2[1] = el2;

      d[second_bc[k]] += loc_fS2[0];
      d[second_bc[k + 1]] += loc_fS2[1];
   }
}

void FEM::first_cond(double t) {
   Functions func;

   for (int k = 0; k < first_bc.size(); k++) {
      int i = first_bc[k] % nx;
      int j = first_bc[k] / nx;

      A.di[first_bc[k]] = 1.0;
      d[first_bc[k]] = func.u_g(x_nodes[i], y_nodes[j], t);

      for (int h = A.ig[first_bc[k]]; h < A.ig[first_bc[k] + 1]; h++) {
         A.al[h] = 0.0;
      }
   }

   for (int i = 0; i < M; i++) {
      for (int j = A.ig[i]; j < A.ig[i + 1]; j++) {
         for (int k = 0; k < first_bc.size(); k++) {
            if (A.jg[j] == first_bc[k]) A.au[j] = 0.0;
         }
      }
   }
}

void FEM::print_result(int l) {
   double error = 0.0;
   Functions func;

   ofstream ans("answer.csv", ios::app);

   ans << "time" << ";" << "answer" << ";" << "real" << ";" << "error" << ";" << endl;

   ans.precision(20);

   q_real.resize(M);
   for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
         q_real[j * nx + i] = func.u_real(x_nodes[i], y_nodes[j], times[l]);
         error = abs(q_real[j * nx + i] - q[j * nx + i]);
         ans << scientific << times[l] << ";" << q[j * nx + i] << ";" << q_real[j * nx + i] << ";" << error << ";" << endl;
      }
   }

   ans.close();
}
