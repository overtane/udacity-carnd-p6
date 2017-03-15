#include "FusionEKF.h"
#include "Tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// return the argument that has greater absolute value
static double absmax(double a, double b) {
    return fabs(a) > fabs(b) ? a : b;
}

/*
 * Constructor.
 */
FusionEKF::FusionEKF() :

  is_initialized_(false),
  previous_timestamp_(0),
  noise_ax_(9),
  noise_ay_(9),
  tools()
{
  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance
  R_laser_ << .0225,     0,
                  0, .0225;
    
  R_radar_ << .09,     0,    0,
                0, .0009,    0,
                0,     0,  .09;

  // measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // initial state covariance matrix P
  MatrixXd P = MatrixXd(4, 4);
  P << 1, 0,    0,    0,
       0, 1,    0,    0,
       0, 0, 1000,    0,
       0, 0,    0, 1000;

  // the initial transition matrix F
  MatrixXd F = MatrixXd(4, 4);
  F << 1, 0, 1, 0,
       0, 1, 0, 1,
       0, 0, 1, 0,
       0, 0, 0, 1;

  VectorXd x  = VectorXd(4);
  MatrixXd Q  = MatrixXd(4,4);

  // create filter object
  ekf_ = KalmanFilter(x, P, F, H_laser_, R_laser_, Q);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}


void FusionEKF::InitializeMeasurement(const MeasurementPackage &mp) {

    double px = 0.0;
    double py = 0.0;
    double vx = 0.0;
    double vy = 0.0;

    previous_timestamp_ = mp.timestamp_;

    // initialize the state vector
    if (mp.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = mp.raw_measurements_[0];  // distance
      double phi = mp.raw_measurements_[1];  // angle

      px = rho * cos(phi);
      py = rho * sin(phi);
      // TODO: we could calculate vx & vy from rho dot and phi.
      is_initialized_ = true;

    } else if (mp.sensor_type_ == MeasurementPackage::LASER) {
      px = mp.raw_measurements_[0];
      py = mp.raw_measurements_[1];
      is_initialized_ = true;

    } else {
      cout << "Invalid sensor type" << endl;
    }

    // initialize Kalman filter state vector
    ekf_.x_ << px, py, vx, vy;
}

void FusionEKF::CalculateProcessCovariance(double dt) {
	
    // Compute the process covariance matrix Q
    double t2 = dt*dt;
    double t3 = t2*dt/2;
    double t4 = t3*dt/2;
  
    double t4ax = t4*noise_ax_;
    double t4ay = t4*noise_ay_;
    double t3ax = t3*noise_ax_;
    double t3ay = t3*noise_ay_;
    double t2ax = t2*noise_ax_;
    double t2ay = t2*noise_ay_;
    
    ekf_.Q_ << t4ax,    0,    t3ax, 0,
               0,       t4ay, 0,    t3ay,
               t3ax,    0,    t2ax, 0,
               0,       t3ay, 0,    t2ay;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &mp) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    cout << "Initialize" << endl;
    InitializeMeasurement(mp);
    // no prediction for the first measurement
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Compute the time elapsed between the current and previous measurements, in seconds
  double dt = (mp.timestamp_ - previous_timestamp_) / 1000000.0L;	
    
  // skip prediction, if measurement has (about) the same timestamp as previous
  // timestamps are integers and in microseconds
  if (mp.timestamp_ - previous_timestamp_ > 100) {

      // Modify the state transition matrix F so that the time is integrated
      ekf_.F_(0,2) = dt;
      ekf_.F_(1,3) = dt;
    
      CalculateProcessCovariance(dt);
      ekf_.Predict();

      //cout << "Q:  " << ekf_.Q_ << endl;
      //cout << "P:  " << ekf_.P_ << endl;
      //cout << "x:  " << ekf_.x_.transpose() << endl;
      // update prediction timestamp only if we did prediction
  }

  previous_timestamp_ = mp.timestamp_;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (mp.sensor_type_ == MeasurementPackage::RADAR) {
    
    VectorXd h = tools.CalculateRadarMeasurementFunction(ekf_.x_);
    //cout << "h:  " << h.transpose() << endl;
    VectorXd y = mp.raw_measurements_ - h; 
    //cout << "y:  " << y.transpose() << endl;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    //cout << "Hj:  " << Hj_ << endl;
    ekf_.Update(y, Hj_, R_radar_);

  } else if (mp.sensor_type_ == MeasurementPackage::LASER) {

    VectorXd y = mp.raw_measurements_ - (H_laser_ * ekf_.x_); 
    ekf_.Update(y, H_laser_, R_laser_);
  }

#if 0
  // print the output
  cout << "x:  " << ekf_.x_.transpose() << endl;
  cout << "P:  " << ekf_.P_ << endl;
#endif
}
