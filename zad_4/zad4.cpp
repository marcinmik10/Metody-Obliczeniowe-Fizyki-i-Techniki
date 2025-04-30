#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
const double m = 1.0;
const double l1 = 1.0;
const double l2 = 1.0 / std::sqrt(8);
const double t_max = 100.0; 
const double dx = 0.001; 

double potential(double x) {
    return -std::exp(-x*x / (l1*l1)) - 8 * std::exp(-std::pow(x - 2, 2) / (l2*l2));
}

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
            file << t << " " << x << " " << v <<" "<<m*v*v/2+potential(x)<< std::endl;

        }

        dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
        if (dt < 1e-8) dt = 1e-8;
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
            file << t << " " << x << " " << v <<" "<<m*v*v/2+potential(x)<< std::endl;

        }

        dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
        if (dt < 1e-8) dt = 1e-8;
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

        double err_x = std::abs((x_small - x_big) / (std::pow(2, d) - 1));
        double err_v = std::abs((v_small - v_big) / (std::pow(2, d) - 1));
        double err = std::max(err_x, err_v);

        if (err <= tol) {
            x = x_small;
            v = v_small;
            t += 2 * dt;
            file << t << " " << x << " " << v <<" "<<m*v*v/2+potential(x)<< std::endl;

        } else {
            // Za duży błąd – nie przesuwamy t
        }

        dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
        if (dt < 1e-8) dt = 1e-8; // ograniczenie minimalnego kroku
    }
}

void euler_damped(std::ofstream &file, double alpha) {
    double x = 2.8, v = 0.0, t = 0.0;
    double dt = 0.0001;
    const double tol = 1e-6, c = 0.9;
    const int d = 1;

    while (t <= t_max) {
        double a = force(x) / m - alpha * v;

        // Duży krok (2*dt)
        double x_big = x + v * 2 * dt;
        double v_big = v + a * 2 * dt;

        // Dwa małe kroki
        double x_half = x + v * dt;
        double v_half = v + a * dt;
        double a_half = force(x_half) / m - alpha * v_half;
        double x_small = x_half + v_half * dt;
        double v_small = v_half + a_half * dt;

        double err_x = std::abs((x_small - x_big) / (std::pow(2, d) - 1));
        double err_v = std::abs((v_small - v_big) / (std::pow(2, d) - 1));
        double err = std::max(err_x, err_v);

        if (err <= tol) {
            x = x_small;
            v = v_small;
            t += 2 * dt;
            file << t << " " << x << " " << v <<" "<<m*v*v/2+potential(x)<< std::endl;

        }

        dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
        // if (dt < 1e-8) dt = 1e-8; // ograniczenie minimalnego kroku
    }
}

void rk4_damped(std::ofstream &file, double alpha) {
    double x = 2.8, v = 0.0, t = 0.0;
    double dt = 0.01;
    const double tol = 1e-6, c = 0.9;
    const int d = 4;

    while (t <= t_max) {
        auto acceleration = [&](double x, double v) {
            return force(x) / m - alpha * v;
        };

        // Duży krok (2*dt)
        double k1x = v;
        double k1v = acceleration(x, v);
        double k2x = v + dt * k1v;
        double k2v = acceleration(x + dt * k1x, v + dt * k1v);
        double x_big = x + dt * (k1x + k2x) / 2.0;
        double v_big = v + dt * (k1v + k2v) / 2.0;

        // Dwa małe kroki
        double x1 = x, v1 = v;
        for (int i = 0; i < 2; ++i) {
            double k1x = v1;
            double k1v = acceleration(x1, v1);
            double k2x = v1 + 0.5 * dt * k1v;
            double k2v = acceleration(x1 + 0.5 * dt * k1x, v1 + 0.5 * dt * k1v);
            double k3x = v1 + 0.5 * dt * k2v;
            double k3v = acceleration(x1 + 0.5 * dt * k2x, v1 + 0.5 * dt * k2v);
            double k4x = v1 + dt * k3v;
            double k4v = acceleration(x1 + dt * k3x, v1 + dt * k3v);

            x1 += dt / 6.0 * (k1x + 2*k2x + 2*k3x + k4x);
            v1 += dt / 6.0 * (k1v + 2*k2v + 2*k3v + k4v);
        }

        double err_x = std::abs((x1 - x_big) / (std::pow(2, d) - 1));
        double err_v = std::abs((v1 - v_big) / (std::pow(2, d) - 1));
        double err = std::max(err_x, err_v);

        if (err <= tol) {
            x = x1;
            v = v1;
            t += 2 * dt;
            file << t << " " << x << " " << v <<" "<<m*v*v/2+potential(x)<< std::endl;

        }

        dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
        // if (dt < 1e-8) dt = 1e-8; // ograniczenie minimalnego kroku
    }
}
void trapezoidal(std::ofstream &file, double alpha) {
    double x = 2.8, v = 0.0, t = 0.0;
    double dt = 0.01;
    const double tol = 1e-6, c = 0.9;
    const int d = 2;
    const int max_iter = 20;

    while (t <= t_max) {
        // Jedno duże przejście (2dt)
        auto step = [&](double x0, double v0, double dt_step) {
            double x1 = x0;
            double v1 = v0;

            for (int iter = 0; iter < max_iter; ++iter) {
                double F1 = x1 - x0 - 0.5 * dt_step * (v1 + v0);
                double F2 = v1 - v0 - 0.5 * dt_step * ((force(x1) / m - alpha * v1) + (force(x0) / m - alpha * v0));

                // Przybliżenie drugiej pochodnej potencjału
                double d2phi = (potential(x1 + dx) - 2 * potential(x1) + potential(x1 - dx)) / (dx * dx);

                double J11 = 1;
                double J12 = -0.5 * dt_step;
                double J21 = -0.5 * dt_step * d2phi / m;
                double J22 = 1 + 0.5 * dt_step * alpha;

                double det = J11 * J22 - J12 * J21;
                if (std::abs(det) < 1e-12) break;

                double dx_corr = (-F1 * J22 + F2 * J12) / det;
                double dv_corr = (-J11 * F2 + J21 * F1) / det;

                x1 += dx_corr;
                v1 += dv_corr;

                if (std::abs(dx_corr) < 1e-10 && std::abs(dv_corr) < 1e-10) break;
            }

            return std::make_pair(x1, v1);
        };

        auto [x_big, v_big] = step(x, v, 2 * dt);

        auto [x_half, v_half] = step(x, v, dt);
        auto [x_small, v_small] = step(x_half, v_half, dt);

        double err_x = std::abs((x_small - x_big) / (std::pow(2, d) - 1));
        double err_v = std::abs((v_small - v_big) / (std::pow(2, d) - 1));
        double err = std::max(err_x, err_v);

        if (err <= tol) {
            x = x_small;
            v = v_small;
            t += 2 * dt;
            file << t << " " << x << " " << v <<" "<<m*v*v/2+potential(x)<< std::endl;

        }

        dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
        if (dt < 1e-8) dt = 1e-8;
    }
}



int main(int argc, char ** argv) {
    // if(strcmp(argv[1], "euler") == 0) {
    //     std::ofstream euler_adaptive_file("euler_adaptive.txt");
    //     euler(euler_adaptive_file);
    //     euler_adaptive_file.close();
    // }
    // else if(strcmp(argv[1],"verlet")==0){
    //     std::ofstream verlet_file("verlet_adaptive.txt");
    //     verlet(verlet_file);
    //     verlet_file.close();
    // }
    // else if(strcmp(argv[1],"rk4")==0){
    //     std::ofstream rk4_file("rk4_adaptive.txt");
    //     rk4(rk4_file);
    //     rk4_file.close();
    // }
    // else if(strcmp(argv[1], "euler_damped") == 0) {
    //     std::ofstream file("euler_damped.txt");
    //     euler_damped(file, 0.5); // lub 5.0
    //     file.close();
    // }
    // else if(strcmp(argv[1], "rk4_damped") == 0) {
    //     std::ofstream file("rk4_damped.txt");
    //     rk4_damped(file, 0.5); 
    //     file.close();
    // }
    // else 
    if(strcmp(argv[1], "trapez") == 0){
        std::ofstream file("trapezoidal_alpha05.txt");
        trapezoidal(file, 0.5);
    }    
    else{
        std::cerr << "Unknown method: " << argv[1] << std::endl;
        return 0;
    }

    std::cout << "Symulacja zakończona. Wyniki zapisane do plików." << std::endl;
    return 0;
}
