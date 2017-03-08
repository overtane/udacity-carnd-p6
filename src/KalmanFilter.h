//
//  KalmanFilter.h
//  Extended Kalman Filters
//
//  Created by Olli Vertanen on 06/03/17.
//

#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class KalmanFilter {
public:
    
    ///* state vector
    VectorXd x_;
    
    ///* state covariance matrix
    MatrixXd P_;
    
    ///* state transistion matrix
    MatrixXd F_;
    
    ///* process covariance matrix
    MatrixXd Q_;
    
    ///* measurement matrix
    MatrixXd H_;
    
    ///* measurement covariance matrix
    MatrixXd R_;
    
    /**
     * Constructors
     */
    KalmanFilter(); // default constructor

    KalmanFilter(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, // initializing contructor
                 MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in);
    
    /**
     * Destructor
     */
    virtual ~KalmanFilter();
    
    /**
     * Predict Predicts the state and the state covariance
     * using the process model
     */
    void Predict();
    
    /**
     * Updates the state and the state covariance
     *
     * Method uses the latest measurement, and the instance variables for calculation. 
     *
     * @param z The measurement at k+1
     */
    void Update(const VectorXd &z);

    /*
     * Updates the state and the state covariance.
     *
     * This version passes in sensor dependent vectors and matrices, so this can be
     * used for all sensor types.
     *
     * Method uses pre-calculated delta vector y, which is the difference of measurement vector k+1
     * and state prediction k+1.
     *  
     * Also, a sensor type dependent measurement matrix is passed in. For linear measurement models, 
     * this is a fixed matrix, and for non-linear a calculated Jacobian matrix.
     *
     *
     * @param y delta of the measurement k+1 and predicted state k+1
     * @param H measurement matrix for state update
     * @param R noise covariance matrix of the sensor
     */
    void Update(const VectorXd &y, const MatrixXd &H, const MatrixXd &R);
    
private:

    MatrixXd I_;

};

#endif /* KALMAN_FILTER_H_ */

