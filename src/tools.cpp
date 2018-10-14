#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;

namespace {

    VectorXd squareComponents(const VectorXd &vector) {
        return vector.array().square().matrix();
    }

}

namespace Tools {

    VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                           const vector<VectorXd> &ground_truth) {
        VectorXd rmse(4);
        rmse << 0, 0, 0, 0;

        // sanity check
        if (estimations.size() != ground_truth.size() || estimations.empty()) {
            cerr << "Invalid estimation or ground_truth data" << endl;
            return rmse;
        }

        // accumulate squared residuals
        for (auto i = 0U; i < estimations.size(); ++i) {
            rmse += squareComponents(estimations[i] - ground_truth[i]);
        }

        // calculate the square root of the mean
        rmse /= estimations.size();
        rmse = rmse.array().sqrt();
        return rmse;
    }

}
