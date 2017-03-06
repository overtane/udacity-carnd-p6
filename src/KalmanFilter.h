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
    KalmanFilter();
    
    KalmanFilter(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
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
     * @param z The measurement at k+1
     */
    void Update(const VectorXd &z);

    void Update(const VectorXd &z, const MatrixXd &H, const MatrixXd &R);
    
private:

    MatrixXd I_;

};

#endif /* KALMAN_FILTER_H_ */

