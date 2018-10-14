#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"

class MeasurementPackage {
public:
    /**
     * Measurement timestamp in microseconds.
     */
    long long timestamp_;

    /**
     * Type of the measurement, e.g. LIDAR or RADAR.
     */
    enum class SensorType{
        LASER,
        RADAR
    } sensor_type_;

    /**
     * The raw measurement values to be interpreted according to sensor type.
     */
    Eigen::VectorXd raw_measurements_;
};

#endif /* MEASUREMENT_PACKAGE_H_ */
