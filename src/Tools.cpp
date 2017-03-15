#include <iostream>
#include "Tools.h"

// return the argument that has greater absolute value
static double absmax(double a, double b) {
    return fabs(a) > fabs(b) ? a : b;
}

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if (estimations.size() == 0 || (estimations.size() != ground_truth.size())) {
        return rmse;
    }
    
    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i) {
        
        VectorXd residual = estimations[i] - ground_truth[i];
        
        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }
    
    //calculate the mean
    rmse = rmse / estimations.size();
    //calculate the squared root
    rmse = rmse.array().sqrt();

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

    MatrixXd Hj(3,4);

    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    //check division by zero
    if (fabs(px)<0.0001 && fabs(py)<0.0001) {
	// handle possible div by zero situation
	// initialize to limits px->0, py->0
        Hj <<       0,    0, 0, 0,
                 1e+9, 1e+9, 0, 0,
                    0,    0, 0, 0;
    } else {
        //compute the Jacobian matrix
        float dvsr2 = px*px + py*py;
        float dvsr1 = sqrt(dvsr2);
        float dvsr3 = dvsr1 * dvsr2;
        
        Hj << px/dvsr1,               py/dvsr1,               0,        0,
              -py/dvsr2,              px/dvsr2,               0,        0,
              py*(vx*py-vy*px)/dvsr3, px*(vy*px-vx*py)/dvsr3, px/dvsr1, py/dvsr1;
    }
    
    return Hj;
}

VectorXd Tools::CalculateRadarMeasurementFunction(VectorXd &x) {    

    VectorXd h(3);

    double px = x(0);
    double py = x(1);
    double vx = x(2);
    double vy = x(3);

    // compute radar h(x) from predicted state 
    // in case of zero position, use an estimate close to zero in order to avoid div by zero
    h(0) = absmax(0.00001, sqrt(px*px + py*py)); // r
    h(1) = atan2(py,px);            // phi
    h(2) = (px*vx + py*vy ) / h(0); // r_dot
 
    return h;
}
