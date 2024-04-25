#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <stdexcept>

class Ensemble {
private:
    int dim;
    std::vector<std::vector<int>> arr;
    double coupling_const;
    double beta;
    double mag_field;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
    std::uniform_int_distribution<int> int_distribution;

    double calc_site_energy(int i, int j) {
        double site_energy = 0.0;
        int nrows = arr.size();
        int ncols = arr[0].size();

        int offsets[2] = {-1, 1};

        for (auto n : offsets) {
            int x = (i + n + nrows) % nrows;
            int y = (j + n + ncols) % ncols;
            site_energy += arr[x][j] + arr[i][y];
        }

        site_energy *= -coupling_const * arr[i][j];
        return site_energy;
    }

    void update_site(int i, int j) {
        if (arr[i][j] == -1)
            arr[i][j] = 1;
        else if (arr[i][j] == 1)
            arr[i][j] = -1;
        else
            throw std::invalid_argument("Invalid value, cannot update state");
    }

public:
    Ensemble(int dim, double coupling_const, double beta, double mag_field)
        : dim(dim), coupling_const(coupling_const), beta(beta), mag_field(mag_field),
          distribution(0.0, 1.0), int_distribution(0, dim - 1) {
        arr.resize(dim, std::vector<int>(dim));
        for (auto& row : arr)
            for (auto& cell : row)
                cell = (distribution(generator) < 0.5) ? -1 : 1;
    }

    void step() {
        int nrows = arr.size();
        int ncols = arr[0].size();

        int i = int_distribution(generator);
        int j = int_distribution(generator);
        double r = distribution(generator);

        double exponent = calc_site_energy(i, j) * -2.0;
        double accept_ratio = 0;

        if (exponent <= 0)
            accept_ratio = 1;
        else {
            accept_ratio = std::exp(-beta * exponent);
        }

        if (r < accept_ratio) {
            update_site(i, j);
        }
    }

    int size() const {
        return dim * dim;
    }
};

int main(int argc, char* argv[]) {
    int dim = 10; // Default dimension
    if (argc > 1) {
	try {
            dim = std::stoi(argv[1]);
        } catch (const std::invalid_argument& ia) {
            std::cerr << "Invalid argument:" << ia.what() << ". Using default dimension." << std::endl;
        } catch (const std::out_of_range& oor) {
            std::cerr << "Argument out of range:" << oor.what() << ". Using default dimension." << std::endl;
        }
    }

    double coupling_const = 1.0;
    double beta = 1.0;
    double mag_field = 1.0;

    Ensemble my_ensemble(dim, coupling_const, beta, mag_field);

    auto t_start = std::chrono::high_resolution_clock::now();

    for (int _ = 0; _ < my_ensemble.size(); ++_) {
        my_ensemble.step();
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count();

    std::cout << "Total Duration: " << duration_us << " us for " << my_ensemble.size() << " steps\n";

    return 0;
}
