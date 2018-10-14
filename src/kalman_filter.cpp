#include <cmath>
#include <iostream>
#include "tools.h"
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

namespace {
    /**
     * Given the system's state in cartesian coordinates, calculates the state in polar coordinates
     * suitable for comparison with the RADAR measurements.
     * @param state The 2D position and velocity state.
     * @return A vector consisting of the range/Distance (rho), bearing (phi) and radial velocity (rho') of the object.
     */
    VectorXd cartesianToPolar(const VectorXd &state) {
        const auto px = state[0];
        const auto py = state[1];
        const auto vx = state[2];
        const auto vy = state[3];

        const auto rho = sqrt(px*px + py*py);
        const auto phi = atan2(py, px);
        const auto rho_dot = (px * vx + py * vy) / rho;

        VectorXd polar(3);
        polar << rho, phi, rho_dot;
        return polar;
    }

    /**
     * Ensures that a given angle is in the -pi..pi range.
     * @param phi The angle.
     * @return The angle in -pi..pi range.
     */
    float clampAngle(float phi) {
        const auto twoPi = 2 * M_PI;
        while (phi > M_PI) {
            phi -= twoPi;
        }
        while (phi < -M_PI) {
            phi += twoPi;
        }
        return phi;
    }
}

KalmanFilter::KalmanFilter() {
    I_ = MatrixXd::Identity(4, 4);
}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    // Obtain prediction error
    const auto z_pred = H_ * x_;
    const auto y = z - z_pred;

    // Determine Kalman gain
    const auto Ht = H_.transpose();
    const auto PHt = P_ * Ht;
    const auto S = H_ * PHt + R_;
    const auto Si = S.inverse();
    const auto K = PHt * Si;

    // Calculate new estimate
    x_ = x_ + (K * y);
    P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // Obtain the prediction error
    VectorXd y = z - cartesianToPolar(x_);
    y[1] = clampAngle(y[1]);

    // Determine Kalman gain (note that H_ has been externally set to the Jacobian!)
    const auto S = H_ * P_ * H_.transpose() + R_;
    const auto K = P_ * H_.transpose() * S.inverse();

    // Calculate new estimate
    x_ = x_ + (K * y);
    P_ = (I_ - K * H_) * P_;
}
