#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  
  x_ = F_ * x_;
  
  // Transpose F, next state function
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  
  // y = z - H * x'
  VectorXd y = z - (H_ * x_);
  
  // S = H * P' * H-transpose + R
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  
  // K = P' * H-transpose * S-inverse
  MatrixXd Si = S.inverse();
  MatrixXd K = (P_ * Ht) * Si;
  
  // x = x' + Ky
  x_ = x_ + (K * y);
  
  // P = (I - K*H) * P'
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - (K * H_)) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  // Define and Compute radar variables
  float rho = sqrt(px*px + py*py);
  float phi = atan2(py, px);
  float rho_dot = (px*vx + py*vy)/rho;
  
  // Compute y = z - h(x')
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;
  
  // Normalize angles to be between -pi and pi
  while (y(1) < -M_PI || y(1) > M_PI) {
    if (y(1) < -M_PI) {
      y(1) = y(1) + 2*M_PI;
    } else {
      y(1) = y(1) - 2*M_PI;
    }
  }
  
  // S = H * P' * H-transpose + R
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  
  // K = P' * H-transpose * S-inverse
  MatrixXd Si = S.inverse();
  MatrixXd K = (P_ * Ht) * Si;
  
  // x = x' + Ky
  x_ = x_ + (K * y);
  
  // P = (I - K*H) * P'
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - (K * H_)) * P_;
}
