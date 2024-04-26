#include <vector>
#include <functional>
#include <cmath>


// Simpson's rule for integrating a 1D function
double integrate1D(std::function<double(double, const std::vector<double>&)> f, double a, double b, int n, const std::vector<double>& params) {
    double h = (b - a) / n;
    double integral = f(a, params) + f(b, params);

    for (int i = 1; i < n; i += 2) {
        integral += 4 * f(a + i * h, params);
    }

    for (int i = 2; i < n - 1; i += 2) {
        integral += 2 * f(a + i * h, params);
    }

    return integral * h / 3.0;
}

// Double integration using Simpson's rule
double integrate2D(std::function<double(double, double, const std::vector<double>&)> f, double ax, double bx, double ay, double by, int nx, int ny, const std::vector<double>& params) {
    auto inner_func = [&](double y, const std::vector<double>& params) {
        return integrate1D([&](double x, const std::vector<double>& params) { return f(x, y, params); }, ax, bx, nx, params);
    };

    return integrate1D(inner_func, ay, by, ny, params);
}

