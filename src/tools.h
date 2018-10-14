#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

namespace Tools {

    /**
     * Ensures that a given angle is in the -pi..pi range.
     * @param phi The angle.
     * @return The angle in -pi..pi range.
     */
    float ClampAngle(float phi);

    /**
    * A helper method to calculate RMSE.
    */
    Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                                  const std::vector<Eigen::VectorXd> &ground_truth);

    /**
   * A helper method to calculate the Jacobian matrix used for a first-order
    * Taylor series approximation of the nonlinear function that maps the system's
    * state to radar measurements.
   */
    Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd &x_state);

}

#endif /* TOOLS_H_ */
