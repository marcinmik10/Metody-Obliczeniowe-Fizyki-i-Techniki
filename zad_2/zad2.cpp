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

// Metoda Eulera
void euler(std::ofstream &file) {
    double x = 2.8, v = 0.0;
    double dt=0.00001;
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
    double dt = 0.001; //krok czasowy
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
    double dt = 0.01; //krok czasowy
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
// void euler_adaptive(std::ofstream &file) {
//     double x = 2.8, v = 0.0;
//     double t = 0.0;
//     double dt = 0.001;
//     const double tol = 1e-6;
//     const double c = 0.9;
//     const int d = 1;
//     double err=1;
//     double x_small=1;
//     double v_small=1;
//     while (t <= t_max) {
//         while (err>tol){
//             // Jeden duży krok (2*dt)
//             double a = force(x) / m;
//             double x_big = x + v * 2 * dt;
//             double v_big = v + a * 2 * dt;

//             // Dwa małe kroki (dt + dt)
//             double x_half = x + v * dt;
//             double v_half = v + a * dt;
//             double a_half = force(x_half) / m;

//             double x_small = x_half + v_half * dt;//wartosci po 2 malych krokach
//             double v_small = v_half + a_half * dt;

//             // Szacowanie błędu
//             double err_x = std::abs((x_small - x_big) / (std::pow(2, d) - 1));
//             double err_v = std::abs((v_small - v_big) / (std::pow(2, d) - 1));
//             err = std::max(err_x, err_v);

//             // Aktualizacja kroku
//             dt = c * dt * std::pow(tol / err, 1.0 / (d + 1));
//             if (dt < 1e-8) dt = 1e-8; // ograniczenie minimalnego kroku
//             if (t + 2 * dt > t_max) dt = (t_max - t) / 2.0; // końcówka symulacji
//         }
//         // Krok zaakceptowany
//         x = x_small;
//         v = v_small;
//         t += 2 * dt;
//         file << t << " " << x << " " << v << std::endl;}}
void euler_adaptive(std::ofstream &file) {
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
        //if (dt < 1e-8) dt = 1e-8; // ograniczenie minimalnego kroku
        //if (t + 2 * dt > t_max) dt = (t_max - t) / 2.0; // końcówka symulacji
    }
}


int main(int argc, char ** argv) {
    if(strcmp(argv[1], "euler") == 0) {
        std::ofstream euler_file("euler.txt");
        euler(euler_file);
        euler_file.close();
    }
    else if(strcmp(argv[1],"verlet")==0){
        std::ofstream verlet_file("verlet.txt");
        verlet(verlet_file);
        verlet_file.close();
    }
    else if(strcmp(argv[1],"rk4")==0){
        std::ofstream rk4_file("rk4.txt");
        rk4(rk4_file);
        rk4_file.close();
    }
    else if(strcmp(argv[1], "euler_adaptive") == 0) {
        std::ofstream euler_adaptive_file("euler_adaptive.txt");
        euler_adaptive(euler_adaptive_file);
        euler_adaptive_file.close();
    }
    else {
        std::cerr << "Unknown method: " << argv[1] << std::endl;
        return 0;
    }

    std::cout << "Symulacja zakończona. Wyniki zapisane do plików." << std::endl;
    return 0;
}
