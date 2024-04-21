#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

// Function to integrate
double function(double x, double y, const std::vector<double>& params) {
    // Example function: f(x, y) = a*x^2 + b*y^2
    double result = 0.0;
    if (params.size() >= 2) {
        result = params[0] * x * x + params[1] * y * y;
    }
    return result;
}

double function1(double r, double z, const std::vector<double>& params) {
    double fx = 0.0;
    double m = 0.0;
    if (params.size() >= 2) {
        m = pow(r * r + z * z / (params[0] * params[0]), 0.5);
        fx = 2*M_PI*r*pow(M_E, -1.0 * m * m / (params[1] * params[1]))*pow(m/params[2], -1.8);
    }
    return fx;
}

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

int main() {
    double ax = 0.01, bx = 1.0; // Integration limits for x
    double ay = -10.0, by = 10.0; // Integration limits for y
    int nx = 100000; // Number of subdivisions for x
    int ny = 100; // Number of subdivisions for y

    // Parameters for the function
    std::vector<double> params = {0.6, 1.9, 1.0}; // Example coefficients a and b

    // Perform double integration
    double result = integrate2D(function1, ax, bx, ay, by, nx, ny, params);
    std::cout << "Result of double integration: " << result << std::endl;

    return 0;
}
