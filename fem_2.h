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
   matrix Mm, G;                      // ������� ����� � ��������� ��� 
   matrix A;                          // �������� �������
   vector<vector<double>> loc_M;      // ��������� ������� ����
   vector<vector<double>> loc_G;      // ��������� ������� ���������
   vector<double> loc_f;              // ��������� ������ ������ �����
   vector<double> F;                  // ������ ������ �����
   vector<double> q_real;             // ��������� �������� �����
                                      
   vector<double> d;                  // ������ ������ ����� ������. ������

   vector<double> x_nodes, y_nodes;   // ����� �� x � y
   vector<int> nodes;                 // ���������� ������ ����� 
   vector<vector<int>> W;             // ������ "������" ������ ������� ��
                                      
   vector<int> first_bc, second_bc, thrid_bc; // ����, � ������� ��������� ��

   int N, M, nx, ny;                  // ���-�� �������� ���������, ���-�� �����
   int hx, hy, ht;                    // ���� �����
   double betta, lambda, chi, sigma;
public:
   FEM();
   ~FEM();

   int n_times;                       // ���������� ��������� ����� 
   vector<double> q_init_0, q_init_1; // ��������� �������
   vector<double> times;              // ��������� �����
   vector<double> q;                  // ������ ����� (�����)
   
   void read_data();                  // ���� ������
   void making_grid();                // ����������� �����
   void making_time_grid();           // ����������� ��������� �����
   void build_portrait(matrix &A);     // ����������� �������� ������� �� �����
   void glob_M();                     // ���������� ���������� ������� ����
   void glob_G();                     // ���������� ���������� ������� ���������
   void glob_F(double t);             // ���������� ����������� ������� ������ �����
   void time_scheme(double dt);       // ���������� �������� ������ ����� � ������ �
// void glob_matrix();                // ������ ���������� �������
   void add_elem(int i, int j, double elem, matrix &A); // ���������� �������� � ���������� �������
   vector<double> matr_vec_mult(vector<double> &x, matrix &A); // ��������� ������� �� ������

   void read_boundary();              // ���� ������� �������
   void thrid_cond(double t);         // ���� �������� �������� �������
   void second_cond(double t);        // ���� ������� �������� �������
   void first_cond(double t);         // ���� ������� �������� �������

   void print_result(double t);
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

   void ILU_precond(); // �������� LU ������������

   void set_vector();
   vector<double> matr_vec_mult(vector<double> &x, bool flag);
   double dot_product(vector<double> &a, vector<double> &b);

   vector<double> direct(vector<double> &L, vector<double> &b);
   vector<double> direct(vector<double> &L, vector<double> &D, vector<double> &b);
   vector<double> reverse(vector<double> &U, vector<double> &b);
   vector<double> reverse(vector<double> &U, vector<double> &D, vector<double> &b);

   void CGM_precond_ILU();
};

