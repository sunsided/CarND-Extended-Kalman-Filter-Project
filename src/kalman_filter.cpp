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
}

KalmanFilter::KalmanFilter() {
    I_ = MatrixXd::Identity(4, 4);
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    // Obtain prediction error
    const VectorXd z_pred = H_ * x_;
    const VectorXd y = z - z_pred;

    // Determine Kalman gain
    const MatrixXd Ht = H_.transpose();
    const MatrixXd PHt = P_ * Ht;
    const MatrixXd S = H_ * PHt + R_;
    const MatrixXd K = PHt *  S.inverse();

    // Calculate new estimate
    x_ = x_ + (K * y);
    P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // Obtain the prediction error
    const VectorXd z_pred = cartesianToPolar(x_);
    VectorXd y = z - z_pred;
    y[1] = Tools::ClampAngle(y[1]);

    // Determine Kalman gain (note that H_ has been externally set to the Jacobian!)
    const MatrixXd Ht = H_.transpose();
    const MatrixXd PHt = P_ * Ht;
    const MatrixXd S = H_ * PHt + R_;
    const MatrixXd K = PHt * S.inverse();

    // Calculate new estimate
    x_ = x_ + (K * y);
    P_ = (I_ - K * H_) * P_;
}
