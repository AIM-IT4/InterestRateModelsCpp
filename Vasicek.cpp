#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

// Function to simulate the Vasicek model
void simulateVasicek(const double& a, const double& b, const double& sigma,
                      const double& r0, const double& T, const double& dt, const std::string& outputPath) {
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

    // Simulate the Vasicek model
    for (int i = 1; i <= numSteps; ++i) {
        // Update time
        time[i] = i * dt;

        // Generate a random increment
        double dW = normalDistribution(generator);

        // Update the interest rate using the Vasicek SDE
        interestRate[i] = interestRate[i - 1] + (a - b * interestRate[i - 1]) * dt
                          + sigma * std::sqrt(dt) * dW;
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
    // Parameters for the Vasicek model
    double a = 0.1;
    double b = 0.2;
    double sigma = 0.02;
    double r0 = 0.05;
    double T = 1.0;
    double dt = 0.01;
    std::string outputPath = "vasicek_simulation.csv";

    // Simulate the Vasicek model
    simulateVasicek(a, b, sigma, r0, T, dt, outputPath);

    std::cout << "Simulation completed. Results saved to " << outputPath << std::endl;

    return 0;
}
