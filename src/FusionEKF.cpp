#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    // LIDAR measurement covariance matrix
    R_laser_ <<
             0.0225, 0,
             0,      0.0225;

    // RADAR measurement covariance matrix
    R_radar_ <<
             0.09, 0,      0,
             0,    0.0009, 0,
             0,    0,      0.09;

    // LIDAR measurement matrix
    H_laser_ <<
             1, 0, 0, 0,
             0, 1, 0, 0;

    // RADAR Jacobian matrix
    Hj_ <<
        1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 1, 1;

    noise_ax = 5;
    noise_ay = 5;
}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        // The initial transition matrix F_
        ekf_.F_ = MatrixXd(4, 4);
        ekf_.F_ <<
                1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0,
                0, 0, 0, 1;

        // State covariance matrix P
        ekf_.P_ = MatrixXd(4, 4);
        ekf_.P_ <<
                1, 0, 0,    0,
                0, 1, 0,    0,
                0, 0, 1000, 0,
                0, 0, 0,    1000;

        if (measurement_pack.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            const auto rho = measurement_pack.raw_measurements_(0);
            const auto phi = measurement_pack.raw_measurements_(1);
            const auto rho_dot = measurement_pack.raw_measurements_(2);
            const auto cos_phi = cos(phi);
            const auto sin_phi = sin(phi);
            ekf_.x_ <<
                    rho * cos_phi,
                    rho * sin_phi,
                    rho_dot * cos_phi,
                    rho_dot * sin_phi;
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::SensorType::LASER) {
            /**
            Initialize state.
            */
            ekf_.x_ <<
                    measurement_pack.raw_measurements_(0),
                    measurement_pack.raw_measurements_(1),
                    0,
                    0;
        }

        // done initializing, no need to predict or update
        previous_timestamp_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
     TODO:
       * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
       * Update the process noise covariance matrix.
       * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */

    // Determine elapsed time in seconds.
    const auto dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    const auto dt2 = dt*dt;
    const auto dt3 = dt2*dt;
    const auto dt4 = dt3*dt;

    // Update the state transition matrix; specifically, position increases
    // with respect to velocity and time elapsed.
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    // Update the process covariance matrix.
    const auto dt4_4 = dt4 * 0.25;
    const auto dt3_2 = dt3 * 0.5;
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<
            dt4_4*noise_ax,    0,              dt3_2*noise_ax, 0,
            0,                 dt4_4*noise_ay, 0,              dt3_2*noise_ay,
            dt3_2*noise_ax,    0,              dt2*noise_ax,   0,
            0,                 dt3_2*noise_ay, 0,              dt2*noise_ay;

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (measurement_pack.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
        // Radar updates
        cout << "Processing RADAR ...." << endl;
        ekf_.H_ = Tools::CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser updates
        cout << "Processing LIDAR ...." << endl;
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
