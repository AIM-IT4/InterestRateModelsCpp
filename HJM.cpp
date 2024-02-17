#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

// Function to simulate the HJM model
void simulateHJM(const double& alpha, const double& sigma, const double& T, const double& dt, int num_paths) {
    // Set up random number generation
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> normalDistribution(0.0, 1.0);

    // Calculate the number of time steps
    int num_steps = static_cast<int>(T / dt);

    // Initialize vectors to store time and forward rate values
    std::vector<double> time(num_steps + 1);
    std::vector<std::vector<double>> forward_rates(num_paths, std::vector<double>(num_steps + 1));

    // Simulate the HJM model
    for (int i = 1; i <= num_steps; ++i) {
        // Update time
        time[i] = i * dt;

        // Generate a random increment
        double dW = normalDistribution(generator);

        // Update the forward rate using the HJM SDE
        for (int path = 0; path < num_paths; ++path) {
            forward_rates[path][i] = forward_rates[path][i - 1] + alpha * dt + sigma * std::sqrt(dt) * dW;
        }
    }

    // Output the results to a CSV file
    std::ofstream outputFile("hjm_simulation.csv");
    outputFile << "Time,ForwardRate1,ForwardRate2,...,ForwardRateN\n";
    for (int i = 0; i <= num_steps; ++i) {
        outputFile << time[i];
        for (int path = 0; path < num_paths; ++path) {
            outputFile << "," << forward_rates[path][i];
        }
        outputFile << "\n";
    }
    outputFile.close();
}

int main() {
    // Parameters for the HJM model
    double alpha = 0.1;
    double sigma = 0.02;
    double T = 1.0;
    double dt = 0.01;
    int num_paths = 5;

    // Simulate the HJM model
    simulateHJM(alpha, sigma, T, dt, num_paths);

    std::cout << "Simulation completed. Results saved to hjm_simulation.csv" << std::endl;

    return 0;
}
