#include <iostream>
#include <fstream>
#include <cmath>
//do sprawdzenia
const double m = 1.0;
const double l1 = 1.0;
const double l2 = 1.0 / std::sqrt(8);
const double dt = 0.01; //krok czasowy
const double t_max = 100.0; // maksymalny czas symulacji
const double dx = 0.001; // dla obliczenia pochodnej siły

// Funkcja potencjału
double potential(double x) {
    return -std::exp(-x*x / (l1*l1)) - 8 * std::exp(-std::pow(x - 2, 2) / (l2*l2));
}

// Obliczanie siły jako -dφ/dx
double force(double x) {
    return -(potential(x + dx) - potential(x - dx)) / (2 * dx);
}

// Metoda Eulera
void euler(std::ofstream &file) {
    double x = 2.8, v = 0.0;
    for (double t = 0; t <= t_max; t += dt) {
        double a = force(x) / m;
        x += v * dt;
        v += a * dt;
        file << t << " " << x << " " << v << std::endl;
    }
}

// Metoda Verleta
void verlet(std::ofstream &file) {
    double x = 2.8, v = 0.0;
    double a = force(x) / m;
    for (double t = 0; t <= t_max; t += dt) {
        double x_new = x + v * dt + 0.5 * a * dt * dt;
        double a_new = force(x_new) / m;
        double v_new = v + 0.5 * (a + a_new) * dt;
        x = x_new;
        v = v_new;
        a = a_new;
        file << t << " " << x << " " << v << std::endl;
    }
}

// Metoda RK4
void rk4(std::ofstream &file) {
    double x = 2.8, v = 0.0;
    for (double t = 0; t <= t_max; t += dt) {
        double k1_x = v;
        double k1_v = force(x) / m;

        double k2_x = v + 0.5 * dt * k1_v;
        double k2_v = force(x + 0.5 * dt * k1_x) / m;

        double k3_x = v + 0.5 * dt * k2_v;
        double k3_v = force(x + 0.5 * dt * k2_x) / m;

        double k4_x = v + dt * k3_v;
        double k4_v = force(x + dt * k3_x) / m;

        x += dt / 6.0 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
        v += dt / 6.0 * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
        file << t << " " << x << " " << v << std::endl;
    }
}

int main() {
    std::ofstream euler_file("euler.txt"), verlet_file("verlet.txt"), rk4_file("rk4.txt");

    euler(euler_file);
    verlet(verlet_file);
    rk4(rk4_file);

    euler_file.close();
    verlet_file.close();
    rk4_file.close();

    std::cout << "Symulacja zakończona. Wyniki zapisane do plików." << std::endl;
    return 0;
}
