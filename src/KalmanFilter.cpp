//
//  KalmanFilter.cpp
//  Extended Kalman Filters
//
//

#include "KalmanFilter.h"
#include <iostream>

KalmanFilter::KalmanFilter() {
}

KalmanFilter::KalmanFilter(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                  MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) : 
    x_(x_in),
    P_(P_in),
    F_(F_in),
    H_(H_in),
    R_(R_in),
    Q_(Q_in)
{
    I_ = MatrixXd::Identity(x_.size(), x_.size());
}

KalmanFilter::~KalmanFilter() {
    // TODO handle identity vector initialisation in case of default constructor
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

    VectorXd y = z - (H_ * x_);
    Update(y, H_, R_);
}

void KalmanFilter::Update(const VectorXd &y, const MatrixXd &H, const MatrixXd &R) {
    
    // measurement update
    MatrixXd Ht = H.transpose();
    MatrixXd S = H * P_ * Ht + R;
    //std::cout << "S:  " <<  S << std::endl;
    //MatrixXd Si = S.inverse();
    //std::cout << "Si:  " << Si << std::endl;
    MatrixXd K = P_ * Ht * S.inverse();
    //std::cout << "K:  " << K << std::endl;

    // new estimate
    x_ = x_ + (K * y);
    P_ = (I_ - K * H) * P_;

    // update H and R for the object.
    H_ = H;
    R_ = R;
}



