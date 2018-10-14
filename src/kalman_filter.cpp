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

    /**
    * A helper method to calculate the Jacobian matrix used for a first-order
     * Taylor series approximation of the nonlinear function that maps the system's
     * state to radar measurements.
    */
    MatrixXd calculateJacobian(const VectorXd &x_state) {
        MatrixXd Hj(3, 4);
        // recover state parameters
        const auto px = x_state(0);
        const auto py = x_state(1);
        const auto vx = x_state(2);
        const auto vy = x_state(3);

        // pre-compute a set of terms to avoid repeated calculation
        const auto c1 = px * px + py * py;
        const auto c2 = sqrt(c1);
        const auto c3 = (c1 * c2);

        // Only c1 needs to be sanity checked:
        // If c1 is nonzero, it is always positive, thus c2 is always positive.
        // Since c3 uses both c1 ad c2, it is also always positive.
        if (fabs(c1) < 0.0001) {
            // Recall that px and py are the distances of the object to our sensor.
            // If both values are zero, the object is exactly at our location, which is
            // something that should never happen to begin with (especially in the context
            // of pedestrian tracking).
            std::cerr << "CalculateJacobian () - Error - Division by Zero" << std::endl;
            Hj.setZero();
            return Hj;
        }

        // pre-compute inverses
        const auto ic1 = 1.0f / c1;
        const auto ic2 = 1.0f / c2;
        const auto ic3 = 1.0f / c3;

        // pre-compute terms used twice
        const auto px_c2 = px * ic2;
        const auto py_c2 = py * ic2;

        // compute the Jacobian
        Hj << px_c2, py_c2, 0, 0,
                -(py * ic1), (px * ic1), 0, 0,
                py * (vx * py - vy * px) * ic3, px * (px * vy - py * vx) * ic3, px_c2, py_c2;

        return Hj;
    }
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

void KalmanFilter::Predict(const double delta_T) {
    assert(delta_T == TIMESTEPS_NOT_USED);
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
    const auto x_size = x_.size();
    const auto I = MatrixXd::Identity(x_size, x_size); // TODO: State is of fixed size, so this can be precalculated
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // Obtain the prediction error
    VectorXd y = z - cartesianToPolar(x_);
    y[1] = clampAngle(y[1]);

    // Determine the Jacobian
    const auto Hj = calculateJacobian(x_);

    // Determine Kalman gain
    const auto S = Hj * P_ * Hj.transpose() + R_;
    const auto K = P_ * Hj.transpose() * S.inverse();

    // Calculate new estimate
    const auto I = MatrixXd::Identity(4, 4); // TODO: State is of fixed size, so this can be precalculated
    x_ = x_ + (K * y);
    P_ = (I - K * Hj) * P_;
}
