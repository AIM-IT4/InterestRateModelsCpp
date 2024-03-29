# Interest Rate Models in C++

This repository contains C++ code implementations for various interest rate models. The models covered include:

1. **Vasicek Model**
   - Description: The Vasicek model captures the mean reversion property of interest rates based on the Ornstein–Uhlenbeck process.
   - Code: [vasicek_model.cpp](Vasicek.cpp)

2. **Cox–Ingersoll–Ross (CIR) Model**
   - Description: The CIR model addresses the positivity problem encountered with the Vasicek model, introducing a nonconstant volatility.
   - Code: [cir_model.cpp](CIR.cpp)

3. **Constant Elasticity of Variance (CEV) Model**
   - Description: The CEV model accounts for nonconstant volatilities that can vary as a power of the underlying asset price.
   - Code: [cev_model.cpp](CEV.cpp)

4. **Hull-White (HW) Model**
   - Description: The Hull-White model is an extension of the Vasicek model, introducing time-dependent interest rate volatility.
   - Code: [hull_white_model.cpp](HW.cpp)

5. **Ho and Lee (HoLee) Model**
   - Description: The Ho and Lee model is based on a deterministic function of time, extending the Merton model.
   - Code: [ho_lee_model.cpp](HoLee.cpp)

6. **Chan–Karolyi–Longstaff–Sanders (CKLS) Model**
   - Description: The CKLS model is a parametrized interest rate model designed to account for nonconstant volatilities.
   - Code: [ckls_model.cpp](CKLS.cpp)
   - 
7. **Heath–Jarrow–Morton (HJM) Model**
   - Description: The HJM model represents the instantaneous forward rate with a stochastic differential equation, allowing for a multi-factor term structure model.
   - Code: [hjm_model.cpp](HJM.cpp)


## Usage

Each model has its own C++ file in the repository. To use a specific model, you can refer to the corresponding C++ file and adjust the parameters as needed for your application.

## Instructions

1. Clone the repository:

```bash
git clone https://github.com/AIM-IT4/interest-rate-models-cpp.git
