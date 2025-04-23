#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
const double m = 1.0;
const double l1 = 1.0;
const double l2 = 1.0 / std::sqrt(8);
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


void verlet(std::ofstream &file) {
    double x = 2.8, v = 0.0, t = 0.0;
    double dt = 0.001;
    const double tol = 1e-6, c = 0.9;
    const int d = 2;

    while (t <= t_max) {
        double a = force(x) / m;

        // Duży krok (2dt)
        double x_big = x + v * 2*dt + 0.5 * a * 4*dt*dt;
        double a_big = force(x_big) / m;
        double v_big = v + 0.5 * (a + a_big) * 2*dt;

        // Dwa małe kroki
        double x_half = x + v * dt + 0.5 * a * dt * dt;
        double a_half = force(x_half) / m;
        double v_half = v + 0.5 * (a + a_half) * dt;

        double x_small = x_half + v_half * dt + 0.5 * a_half * dt * dt;
        double a_small = force(x_small) / m;
        double v_small = v_half + 0.5 * (a_half + a_small) * dt;

        double err_x = std::abs((x_small - x_big) / (std::pow(2, d) - 1));
        double err_v = std::abs((v_small - v_big) / (std::pow(2, d) - 1));
        double err = std::max(err_x, err_v);

        if (err <= tol) {
            x = x_small;
            v = v_small;
            t += 2 * dt;
            file << t << " " << x << " " << v << " " << dt << std::endl;
        }

        dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
        if (dt < 1e-8) dt = 1e-8;
        // if (t + 2 * dt > t_max) dt = (t_max - t) / 2.0;
    }
}

void rk4(std::ofstream &file) {
    double x = 2.8, v = 0.0, t = 0.0;
    double dt = 0.01;
    const double tol = 1e-6, c = 0.9;
    const int d = 4;

    while (t <= t_max) {
        auto rk4_step = [](double x, double v, double dt) {
            double k1_x = v;
            double k1_v = force(x) / m;

            double k2_x = v + 0.5 * dt * k1_v;
            double k2_v = force(x + 0.5 * dt * k1_x) / m;

            double k3_x = v + 0.5 * dt * k2_v;
            double k3_v = force(x + 0.5 * dt * k2_x) / m;

            double k4_x = v + dt * k3_v;
            double k4_v = force(x + dt * k3_x) / m;

            double x_next = x + dt / 6.0 * (k1_x + 2*k2_x + 2*k3_x + k4_x);
            double v_next = v + dt / 6.0 * (k1_v + 2*k2_v + 2*k3_v + k4_v);
            return std::make_pair(x_next, v_next);
        };

        // Duży krok
        auto [x_big, v_big] = rk4_step(x, v, 2 * dt);

        // Dwa małe kroki
        auto [x_half, v_half] = rk4_step(x, v, dt);
        auto [x_small, v_small] = rk4_step(x_half, v_half, dt);

        double err_x = std::abs((x_small - x_big) / (std::pow(2, d) - 1));
        double err_v = std::abs((v_small - v_big) / (std::pow(2, d) - 1));
        double err = std::max(err_x, err_v);

        if (err <= tol) {
            x = x_small;
            v = v_small;
            t += 2 * dt;
            file << t << " " << x << " " << v << " " << dt << std::endl;
        }

        dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
        if (dt < 1e-8) dt = 1e-8;
        // if (t + 2 * dt > t_max) dt = (t_max - t) / 2.0;
    }
}



void euler(std::ofstream &file) {
    double x = 2.8, v = 0.0;
    double t = 0.0;
    double dt = 0.0001;
    const double tol = 1e-6;
    const double c = 0.9;
    const int d = 1;

    while (t <= t_max) {
        // Jeden duży krok (2*dt)
        double a = force(x) / m;
        double x_big = x + v * 2 * dt;
        double v_big = v + a * 2 * dt;

        // Dwa małe kroki (dt + dt)
        double x_half = x + v * dt;
        double v_half = v + a * dt;
        double a_half = force(x_half) / m;

        double x_small = x_half + v_half * dt;
        double v_small = v_half + a_half * dt;

        // Szacowanie błędu
        double err_x = std::abs((x_small - x_big) / (std::pow(2, d) - 1));
        double err_v = std::abs((v_small - v_big) / (std::pow(2, d) - 1));
        double err = std::max(err_x, err_v);

        if (err <= tol) {
            // Krok zaakceptowany
            x = x_small;
            v = v_small;
            t += 2 * dt;
            file << t << " " << x << " " << v << " " << dt << std::endl;
        } else {
            // Za duży błąd – nie przesuwamy t
          // std::cout << "Za duży błąd (" << err << ") przy t = " << t << ", zmniejszamy dt\n";
        }

        // Aktualizacja kroku
        dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
        if (dt < 1e-8) dt = 1e-8; // ograniczenie minimalnego kroku
        //if (t + 2 * dt > t_max) dt = (t_max - t) / 2.0; // końcówka symulacji
    }
}


int main(int argc, char ** argv) {
    if(strcmp(argv[1], "euler") == 0) {
        std::ofstream euler_adaptive_file("euler_adaptive.txt");
        euler(euler_adaptive_file);
        euler_adaptive_file.close();
    }
    else if(strcmp(argv[1],"verlet")==0){
        std::ofstream verlet_file("verlet_adaptive.txt");
        verlet(verlet_file);
        verlet_file.close();
    }
    else if(strcmp(argv[1],"rk4")==0){
        std::ofstream rk4_file("rk4_adaptive.txt");
        rk4(rk4_file);
        rk4_file.close();
    }
    else {
        std::cerr << "Unknown method: " << argv[1] << std::endl;
        return 0;
    }

    std::cout << "Symulacja zakończona. Wyniki zapisane do plików." << std::endl;
    return 0;
}
