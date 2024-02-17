#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

// Function to simulate the CKLS model
void simulateCKLS(const double& alpha, const double& beta, const double& gamma, const double& sigma,
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

    // Simulate the CKLS model
    for (int i = 1; i <= numSteps; ++i) {
        // Update time
        time[i] = i * dt;

        // Generate a random increment
        double dW = normalDistribution(generator);

        // Update the interest rate using the CKLS SDE
        interestRate[i] = interestRate[i - 1] + (alpha - beta * interestRate[i - 1]) * dt
                          + sigma * std::pow(std::abs(interestRate[i - 1]), gamma) * std::sqrt(dt) * dW;
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
    // Parameters for the CKLS model
    double alpha = 0.1;
    double beta = 0.2;
    double gamma = 0.5;
    double sigma = 0.02;
    double r0 = 0.05;
    double T = 1.0;
    double dt = 0.01;
    std::string outputPath = "ckls_simulation.csv";

    // Simulate the CKLS model
    simulateCKLS(alpha, beta, gamma, sigma, r0, T, dt, outputPath);

    std::cout << "Simulation completed. Results saved to " << outputPath << std::endl;

    return 0;
}
