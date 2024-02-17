#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

// Function to simulate the Hull and White model
void simulateHullAndWhite(const double& r0, const std::vector<double>& theta,
                           const std::vector<double>& alpha, const std::vector<double>& sigma,
                           const double& T, const double& dt, const std::string& outputPath) {
    // Set up random number generation
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> normalDistribution(0.0, 1.0);

    // Calculate the number of time steps
    int numSteps = static_cast<int>(T / dt);

    // Initialize vectors to store time and interest rate values
    std::vector<double> time(numSteps + 1);
    std::vector<double> interestRate(numSteps + 1);

    // Initialize the initial interest rate
    interestRate[0] = r0;

    // Simulate the Hull and White model
    for (int i = 1; i <= numSteps; ++i) {
        // Update time
        time[i] = i * dt;

        // Generate random increments
        double dW = normalDistribution(generator);

        // Calculate the integral terms for the explicit solution
        double integral1 = 0.0;
        double integral2 = 0.0;
        double integral3 = 0.0;

        for (int j = 0; j < i; ++j) {
            double u = j * dt;
            integral1 += alpha[j] * dt;
            integral2 += theta[j] * std::exp(-integral1);
            integral3 += sigma[j] * std::exp(-integral1) * dW;
        }

        // Update the interest rate using the Hull and White SDE
        interestRate[i] = r0 * std::exp(-integral1) + integral2 + integral3;
    }

    // Output the results to a CSV file
    std::ofstream outputFile(outputPath);
    outputFile << "Time,InterestRate\n";
    for (int i = 0; i <= numSteps; ++i) {
        outputFile << time[i] << "," << interestRate[i] << "\n";
    }
    outputFile.close();
}

int main() {
    // Parameters for the Hull and White model
    double r0 = 0.02;  // Replace with the desired initial interest rate
    std::vector<double> theta = {0.03, 0.02, 0.025};  // Replace with the desired deterministic function of time
    std::vector<double> alpha = {0.01, 0.015, 0.012};  // Replace with the desired deterministic function of time
    std::vector<double> sigma = {0.01, 0.015, 0.02};    // Replace with the desired deterministic function of time
    double T = 1.0;
    double dt = 0.01;
    std::string outputPath = "hull_and_white_simulation.csv";

    // Simulate the Hull and White model
    simulateHullAndWhite(r0, theta, alpha, sigma, T, dt, outputPath);

    std::cout << "Simulation completed. Results saved to " << outputPath << std::endl;

    return 0;
}
