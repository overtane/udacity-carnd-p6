#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "MeasurementPackage.h"
#include <vector>
#include <string>
#include <fstream>
#include "KalmanFilter.h"
#include "Tools.h"

class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
   * Initialize state vector with first valid measurement
   */
  void InitializeMeasurement(const MeasurementPackage &mp);
  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
   * Calculate process covariance matrix Q
   */
  void CalculateProcessCovariance(double dt);

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;

private:
  // check whether the tracking toolbox was initiallized or not (first measurement)
  bool is_initialized_;

  // previous prediction timestamp (timestamp of the measurement that was used for prediction)
  long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  
  MatrixXd R_laser_; // laser measurement covariance matrix
  MatrixXd R_radar_; // radar measurement covariance matrix
  MatrixXd H_laser_; // laser measurement matrix
  MatrixXd Hj_;

  float noise_ax_;
  float noise_ay_;

};

#endif /* FusionEKF_H_ */
