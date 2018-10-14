#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"


class Tools {
public:
    /**
    * Constructor.
    */
    Tools() = default;

    /**
    * Destructor.
    */
    virtual ~Tools() = default;

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

};

#endif /* TOOLS_H_ */
