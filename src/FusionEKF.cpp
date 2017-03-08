#include "FusionEKF.h"
#include "Tools.h"
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

  // measurement covariance
  R_laser_ << .15, 0,
              0, .15;
    
  R_radar_ << .3, 0,  0,
              0, .03, 0,
              0, 0,   .3;

  // measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;


  // initial state covariance matrix P
  MatrixXd P = MatrixXd(4, 4);
  P << 1, 0, 0, 0,
       0, 1, 0, 0,
       0, 0, 10000, 0,
       0, 0, 0, 10000;

  // the initial transition matrix F
  MatrixXd F = MatrixXd(4, 4);
  F << 1, 0, 1, 0,
       0, 1, 0, 1,
       0, 0, 1, 0,
       0, 0, 0, 1;

  VectorXd x  = VectorXd(4);
  MatrixXd Q  = MatrixXd(4,4);

  // create the filter
  ekf_ = KalmanFilter(x, P, F, H_laser_, R_laser_, Q);

  noise_ax_ = 50; 
  noise_ay_ = 50;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    previous_timestamp_ = measurement_pack.timestamp_;

    // initialize the state vector
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      ekf_.x_ << ro * cos(phi), ro * sin(phi), 0, 0;
      is_initialized_ = true;

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      is_initialized_ = true;

    } else {
        cout << "Invalid sensor type" << endl;
    }

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Compute the time elapsed between the current and previous measurements, in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	
    
  // Modify the state transition matrix F so that the time is integrated
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
    
  // skip prediction, if measurement has (about) the same timestamp as previous
  // timestamps are microseconds
  if (previous_timestamp_ < measurement_pack.timestamp_ - 100) {

    float t2 = dt * dt;
    float t3 = (t2 * dt)/2;
    float t4 = (t3 * dt)/4;
  
    float t4ax = t4*noise_ax_;
    float t4ay = t4*noise_ay_;
    float t3ax = t3*noise_ax_;
    float t3ay = t3*noise_ay_;
    float t2ax = t2*noise_ax_;
    float t2ay = t2*noise_ay_;
    
    // Compute the process covariance matrix Q
    ekf_.Q_ << t4ax,    0,    t3ax, 0,
               0,       t4ay, 0,    t3ay,
               t3ax,    0,    t2ax, 0,
               0,       t3ay, 0,    t2ay;
    
    ekf_.Predict();

    // update prediction timestamp only if we did prediction
    previous_timestamp_ = measurement_pack.timestamp_;
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    
    float px = ekf_.x_(0);
    float py = ekf_.x_(1);
    float vx = ekf_.x_(2);
    float vy = ekf_.x_(3);

    // compute radar h(x) from predicted state 
    VectorXd h(3);
    h(0) = sqrt(px*px + py*py);         // r
    if (fabs(h(0) > 0.001)) {
        h(1) = atan2(py,px);            // phi
        h(2) = (px*vx + py*vy ) / h(0); // r_dot
    } else {
        h(1) = 0.0;
        h(2) = 0.0;
    }

    Hj_ = tools.CalculateJacobian(ekf_.x_);

    VectorXd y = measurement_pack.raw_measurements_ - h; 
    
    ekf_.Update(y, Hj_, R_radar_);

  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

    VectorXd y = measurement_pack.raw_measurements_ - (H_laser_ * ekf_.x_); 
    ekf_.Update(y, H_laser_, R_laser_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
