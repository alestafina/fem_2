#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

struct matrix {
   vector<int> ig, jg;
   vector<double> al, au, di;
};

class FEM {
protected:
   matrix Mchi, Msigma, G;            // матрицы массы и жесткости для 
   matrix A;                          // итоговая матрица
   vector<vector<double>> loc_M;      // локальная матрица масс
   vector<vector<double>> loc_G;      // локальная матрица жесткости
   vector<double> loc_f;              // локальный вектор правой части
   vector<double> F;                  // вектор правой части
   vector<double> q_real;             // ожидаемые значения весов
                                      
   vector<double> q;                  // вектор весов (ответ)
   vector<double> d;                  // вектор правой части нестац. задачи
   vector<double> q_init_2, q_init_1; // начальные условия
   vector<double> times;              // временная сетка

   vector<double> x_nodes, y_nodes;   // сетка по x и y
   vector<int> nodes;                 // глобальные номера узлов 
   vector<vector<int>> W;             // хранит "адреса" границ каждого КЭ
                                      
   vector<int> first_bc, second_bc, thrid_bc; // узлы, в которых выполнены ку

   int N, M, nx, ny;                  // кол-во конечных элементов, кол-во узлов
   double hx, hy, ht;                 // шаги сетки
   double betta, lambda, chi, sigma;
public:
   FEM();
   ~FEM();

   int n_times;                       // количество временных слоев 
   
   void read_data();                  // ввод данных
   void read_time_nonform();          // чтение неравномерной сетки времени
   void making_grid();                // составление стеки
   void making_time_grid();           // составление временной сетки
   void init_cond();                  // подсчет начальных условий
   
   void build_portrait(matrix &A);    // составление портрета матрицы по сетке
   void glob_M();                     // вычисление глобальной матрицы масс
   void glob_G();                     // вычисление глобальной матрицы жесткости
   void glob_F(double t);             // вычисление глобального вектора правой части
   void time_scheme_d(int i);         // вычисление векторов правой части 
   void time_scheme_A();              // вычисление матрицы А
   void time_scheme_nonform(int i);   // аппроксмация по неравномерн. сетке времени
   void add_elem(int i, int j, double elem, matrix &A); // добавление элемента в глобальную матрицу
   vector<double> matr_vec_mult(vector<double> &x, matrix &A); // умножение матрицы на вектор

   void read_boundary();              // ввод краевых условий
   void thrid_cond(double t);         // учет третьего краевого условия
   void second_cond(double t);        // учет второго краевого условия
   void first_cond(double t);         // учет первого краевого условия

   void print_result(int l);
};

class Functions {
public:
   double func(double x, double y, double t);
   double u_g(double x, double y, double t);
   double u_betta(double x, double y, double t);
   double theta(double x, double y, double t);
   double u_real(double x, double y, double t);
};

class CGM : public FEM {
private:
   vector<double> r, z, buf1, buf2;
   vector<double> al_LU, au_LU, di_LU;
   int maxiter, k;
   double eps;
   double residual, norm;
public:
   CGM();
   ~CGM();

   void ILU_precond(); // неполная LU факторизация

   void set_vector();
   vector<double> matr_vec_mult(vector<double> &x, bool flag);
   double dot_product(vector<double> &a, vector<double> &b);

   vector<double> direct(vector<double> &L, vector<double> &b);
   vector<double> direct(vector<double> &L, vector<double> &D, vector<double> &b);
   vector<double> reverse(vector<double> &U, vector<double> &b);
   vector<double> reverse(vector<double> &U, vector<double> &D, vector<double> &b);

   void CGM_precond_ILU();
};

