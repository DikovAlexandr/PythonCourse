#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

double phi(double x) {
    return 5.0 / 2 * tan(x);
}

double psi(double x) {
    return -5 / (2 * pow(cosh(x), 2));
}

double gamma0(double t) {
    return -5 * tan(t) - 5 / (2 * pow(cosh(t), 2));
}

double gammaL(double t) {
    return 5 / (2 * pow(cosh(1 - t), 2));
}

double f(double x, double t) {
    return 5 * tanh(t - x) / pow(cosh(t - x), 2);
}

int main() {
    double L = 1, T = 1; // пределы по x и t
    double h = 0.05, tau = 0.05; // шаги по пространству и по времени
    double a = 0.5;
    std::vector<std::vector<double>> u(int(L / h) + 1, std::vector<double> (int(T / tau) + 1, 0));

    for (int n = 0; n * tau <= T; n++) {
        for (int i = 0; i * h <= L; i++) {
            //std::cout << i << ' ' << n << std::endl;
            if (n == 0)
                u[i][n] = phi(i * h);
            if (n == 1)
                u[i][n] = u[i][0] + tau * psi(i * h);
            if (i != 0 and i * h != L and n * tau < T and n > 1)
                u[i][n + 1] =
                        2 * u[i][n] - u[i][n - 1] + pow(a * tau / h, 2) * (u[i + 1][n] - 2 * u[i][n] + u[i - 1][n]) +
                        pow(tau, 2) * f(i * h, n * tau);
        }
    }

    for (int n = 0; n * tau <= T; n++) {
        for (int i = 0; i * h <= L; i++) {
            if (i == 0)
                u[i][n] = (u[i + 1][n] + gamma0(n * tau)) / (2 - 1 / h);
            if (i == L)
                u[i][n] = (u[i - 1][n] + gammaL(n * tau)) / (1 / h);
        }
    }

    std::ofstream fout;
    fout.open(R"(..\\..\\FEM wave equation\\solve.csv)");

    for (auto &i: u) {
        for (double j: i) {
            fout << std::setprecision(3) << j << " ";
        }
        fout << std::endl;
    }
    fout.close();
    return 0;
}
