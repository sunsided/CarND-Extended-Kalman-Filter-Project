#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;


VectorXd squareComponents(const VectorXd &vector) {
    return vector.array().square().matrix();
}


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // sanity check
    if(estimations.size() != ground_truth.size() || estimations.empty()){
        cerr << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    // accumulate squared residuals
    for(auto i=0; i < estimations.size(); ++i) {
        rmse += squareComponents(estimations[i] - ground_truth[i]);
    }

    // calculate the square root of the mean
    rmse /= estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);
    // recover state parameters
    const auto px = x_state(0);
    const auto py = x_state(1);
    const auto vx = x_state(2);
    const auto vy = x_state(3);

    // pre-compute a set of terms to avoid repeated calculation
    const auto c1 = px*px + py*py;
    const auto c2 = sqrt(c1);
    const auto c3 = (c1*c2);

    // Only c1 needs to be sanity checked:
    // If c1 is nonzero, it is always positive, thus c2 is always positive.
    // Since c3 uses both c1 ad c2, it is also always positive.
    if(fabs(c1) < 0.0001){
        // Recall that px and py are the distances of the object to our sensor.
        // If both values are zero, the object is exactly at our location, which is
        // something that should never happen to begin with (especially in the context
        // of pedestrian tracking).
        cerr << "CalculateJacobian () - Error - Division by Zero" << endl;
        Hj.setZero();
        return Hj;
    }

    // pre-compute inverses
    const auto ic1 = 1.0f / c1;
    const auto ic2 = 1.0f / c2;
    const auto ic3 = 1.0f / c3;

    // pre-compute terms used twice
    const auto px_c2 = px*ic2;
    const auto py_c2 = py*ic2;

    // compute the Jacobian
    Hj << px_c2,                  py_c2,                  0,     0,
        -(py*ic1),               (px*ic1),                0,     0,
          py*(vx*py - vy*px)*ic3, px*(px*vy - py*vx)*ic3, px_c2, py_c2;

    return Hj;
}
