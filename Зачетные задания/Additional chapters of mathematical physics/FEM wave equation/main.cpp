#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

// Параметры задачи
const double L = 1, T = 1; // Пределы по x и t
const double H = 0.05, TAU = 0.05; // Шаги по пространству и по времени
const double A = 0.5;
const int SIZE_X = int(L / H) + 1, SIZE_T = int(T / TAU) + 1;
const double alpha0 = 2, alphaL = 0, beta0 = -1, betaL = 1;

// Начальное условие на координату
double phi(double x) {
    return 5.0 / 2 * tanh(x);
}

// Начальное условие на импульс
double psi(double x) {
    return -5.0 / (2 * pow(cosh(x), 2));
}

// Граничное условие в нуле
double gamma0(double t) {
    return -5 * tanh(t) - 5 / (2 * pow(cosh(t), 2));
}

// Граничное условие в L
double gammaL(double t) {
    return 5 / (2 * pow(cosh(1 - t), 2));
}

// Неоднородность уравнения
double f(double x, double t) {
    return 5.0 / 2 * tanh(t - x) / pow(cosh(t - x), 2);
}

// Подстановка начальных условий
void initial(std::vector<std::vector<double>> &u) {
    for (int i = 0; i < SIZE_X; ++i) {
        u[0][i] = phi(i * H);
        u[1][i] = u[0][i] + TAU * psi(i * H);
    }
}

int main() {
    std::vector<std::vector<double>> u(SIZE_T, std::vector<double>(SIZE_X, 0));
    initial(u);
    for (int n = 2; n < SIZE_T; n++) {
        for (int i = 1; i < SIZE_X; i++) {
            if (i != 0 and i * H != L and n > 1)
                // Шаг схемой "Крест"
                u[n][i] = 2 * u[n - 1][i] - u[n - 2][i] +
                          A * pow(TAU / H, 2) * (u[n - 1][i + 1] - 2 * u[n - 1][i] + u[n - 1][i - 1]) +
                          pow(TAU, 2) * f(i * H, (n - 1) * TAU);

        }
        // Подстановка граничных условий
//        u[n][0] = (u[n][1] / H + gamma0(n * TAU)) / (2 + 1 / H);
//        u[n][SIZE_X - 1] = (u[n][SIZE_X - 2] / H + gammaL(n * TAU)) / (1 / H);
        u[n][0] = (gamma0(n * TAU) * H - beta0 * u[n][1]) / (alpha0 * H - beta0);
        u[n][SIZE_X - 1] = (gammaL(n * TAU) * H + betaL * u[n][SIZE_X - 2]) / (alphaL * H + betaL);
//        u[n][0] = (2 * H * gamma0(n * TAU) - 4 * beta0 * u[n][1] + beta0 * u[n][2]) / (2 * H * alpha0 - 3 * beta0);
//        u[n][SIZE_X - 1] = (2 * H * gammaL(n * TAU) + 4 * betaL * u[n][SIZE_X - 2] - betaL * u[n][SIZE_X - 3]) / (2 * H * alphaL + 3 * betaL);
    }

    // Вывод в файл
    std::ofstream fout;
    fout.open(R"(..\\..\\FEM wave equation\\solve.csv)");
    for (auto &t: u) {
        for (double x: t) {
            fout << std::setprecision(8) << x << " ";
        }
        fout << std::endl;
    }
    fout.close();
    return 0;
}
