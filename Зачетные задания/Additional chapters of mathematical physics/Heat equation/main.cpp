#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

// Параметры расчета
const int T_ = 1;  // Время
const double TAU = 0.01; // Шаг по времени
const int L = 1; // Длина
const double H = 0.01; // Шаг по длине
const double SIGMA = 0.7; // Вес
const int A_ = 1; // Температуропроводность
const int SIZE_X = int(L / H) + 1, SIZE_T = int(T_ / TAU) + 1;

// Начальные условия
double initial(double x) {
    return pow(x, 2) + 1.0 / 2;
}

// Левое граничное условие
std::tuple<double, double, double> left_bord_cond(double t) {
    std::tuple<double, double, double> result(1, 0, pow(t, 2) + 1.0 / 2);
    return result;
}

// Правое граничное условие
std::tuple<double, double, double> right_bord_cond(double t) {
    std::tuple<double, double, double> result(3, 1, 5 - 8 * t + 3 * pow(t, 3) + (3 + t * tan(t)) / (2 * cos(t)));
    return result;
}

// Неоднородность
double inhomogeneity(double x, double t) {
    return -2 * (1 + x + t) + (x * tan(x * t) - pow(t, 2) * (1 + 2 * pow(tan(x * t), 2))) / (2 * cos(x * t));
}

int main() {
    // Инициализация вектора координаты
    std::vector<double> X(SIZE_X, 0);
    for (int i = 1; i < SIZE_X; ++i) {
        X[i] = X[i - 1] + H;
    }

    // Инициализация вектора времени
    std::vector<double> T(SIZE_T, 0);
    for (int i = 1; i < SIZE_T; ++i) {
        T[i] = T[i - 1] + TAU;
    }

    // Инициализация поля температур
    std::vector<double> U(SIZE_X, 0), U_prev(SIZE_X, 0);

    // Заполняем согласно начальным условиям
    for (int i = 0; i < SIZE_X; ++i) {
        U[i] = initial(X[i]);
        for (int j = 0; i < SIZE_T; ++i) {
            // Сохраняем предыдущий шаг
            U_prev = U;
            // Заполняем граничные условия
            std::vector<double> d(SIZE_X, 0);

            std::tuple<double, double, double> left = left_bord_cond(T[j]);
            double alpha_0 = get<0>(left);
            double beta_0 = get<1>(left);
            d[0] = get<2>(left);

            std::tuple<double, double, double> right = right_bord_cond(T[j]);
            double alpha_l = get<0>(right);
            double beta_l = get<1>(right);
            d[SIZE_X - 1] = get<2>(right);

            double A_0 = alpha_0 - beta_0 / H;
            double B_0 = beta_0 / H;
            double A_l = alpha_l - beta_l / H;
            double B_l = beta_l / H;
            double A = - A_ * SIGMA / pow(H, 2);
            double C = A;
            double B = 2 * A_ * SIGMA / pow(H, 2) + 1 / TAU;

            std::vector<double> a(SIZE_X - 1, 0);
            std::vector<double> b(L - 1, 0);
            a[0] = - B_0 / A_0;
            b[0] = d[0] / A_0;
            // Прямой прогон
            for (int k = 1; k < SIZE_X - 1; ++k) {
                d[k] = A_ * (1 - SIGMA) / pow(H, 2) * (U_prev[k + 1] - 2 * U_prev[k] + U_prev[k - 1]) +
                       U_prev[k] / TAU + inhomogeneity(X[k], T[j]);
                a[k] = -C / (A * a[k - 1] + B);
                b[k] = (d[k] - A * b[k - 1]) / (A * a[k - 1] + B);
            }

            U[SIZE_X - 1] = (d[SIZE_X - 1] - A_l * b[SIZE_X - 2]) / (a[SIZE_X - 2] * A_l + B_l);

            // Обратный прогон
            for (int k = SIZE_X - 1; k < 1; --k) {
                U[k-1] = a[k-1] * U[k] + b[k-1];
            }
        }
    }
    std::vector<double> U_0(L, 0);
    for (int i = 0; i < SIZE_X; ++i){
        U_0[i] = pow((X[i] - T[SIZE_T - 1]), 2) + 1 / (2 * cos(X[i] * T[SIZE_T - 1]));
    }

    // Вывод в файл
    std::ofstream fout;
    fout.open(R"(..\\..\\Heat equation\\solve.csv)");
    for (auto &x: U) {
        fout << std::setprecision(8) << x << " ";
    }
    fout.close();

    fout.open(R"(..\\..\\Heat equation\\analytical.csv)");
    for (auto &x: U_0) {
        fout << std::setprecision(8) << x << " ";
    }
    fout.close();
    return 0;
}