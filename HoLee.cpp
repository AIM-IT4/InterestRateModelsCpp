#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

// Function to simulate the Ho and Lee model
void simulateHoAndLee(const double& theta, const double& sigma, const double& T, const double& dt, const std::string& outputPath) {
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
    interestRate[0] = 0.0;

    // Simulate the Ho and Lee model
    for (int i = 1; i <= numSteps; ++i) {
        // Update time
        time[i] = i * dt;

        // Generate a random increment
        double dW = normalDistribution(generator);

        // Update the interest rate using the Ho and Lee SDE
        interestRate[i] = interestRate[i - 1] + theta * time[i] * dt + sigma * std::sqrt(dt) * dW;
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
    // Parameters for the Ho and Lee model
    double theta = 0.02;  // Replace with the desired deterministic function of time
    double sigma = 0.01;
    double T = 1.0;
    double dt = 0.01;
    std::string outputPath = "ho_and_lee_simulation.csv";

    // Simulate the Ho and Lee model
    simulateHoAndLee(theta, sigma, T, dt, outputPath);

    std::cout << "Simulation completed. Results saved to " << outputPath << std::endl;

    return 0;
}
