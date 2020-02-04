#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  // Set H_laser_, Measurement function  matrix
  H_laser_ << 1, 0, 0, 0,
  			 0, 1, 0, 0;
  
  // Create x, State vector
  ekf_.x_ = VectorXd(4);
  
  // Create and Set P, State Covariance matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
  			 0, 1, 0, 0,
  			 0, 0, 1000, 0,
  			 0, 0, 0, 1000;
  
  // Create and Set F, Intitial Transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
  			 0, 1, 0, 1,
  			 0, 0, 1, 0,
  			 0, 0, 0, 1;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    cout << "First State Initialization: " << endl;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      
      // Pull ro, phi, and ro_dot variables from Measurement Package
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];
      
      // Initialize x, State vector
      ekf_.x_ << ro * cos(phi),
      			 ro * sin(phi),
      			 ro_dot * cos(phi),
      			 ro_dot * sin(phi);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      
      	// Initialize x, State vector
      ekf_.x_ << measurement_pack.raw_measurements_[0], 
           		 measurement_pack.raw_measurements_[1], 
              	 0, 
              	 0;
    }
    
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  cout << "Prediction Step" << endl;

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  //  Set noise_ax and noise_ay for Q process noise covariance matrix
  float noise_ax = 9.0;
  float noise_ay = 9.0;
  
  // Compute the change in time between the current and previous measurements in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  
  // Update previous timestamp
  previous_timestamp_ = measurement_pack.timestamp_;
  
  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;
  
  // Update state transition matrix F using new elapsed time
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // Create and Set Q, process noise covariance matrix
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt4 / 4 * noise_ax, 0, dt3 / 2 * noise_ax, 0,
  			 0, dt4 / 4 * noise_ay, 0, dt3 / 2 * noise_ay,
  			 dt3 / 2 * noise_ax, 0, dt2 * noise_ax, 0,
  			 0, dt3 / 2 * noise_ay, 0, dt2 * noise_ay;
  
  // Call Predict function in kalman_filter.cpp
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    cout << "Measurement Update - Radar" << endl;
    Tools tools;
  
  	// Compute Hj Matrix and set H to it for Radar measurement
  	Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // TODO: Laser updates
    cout << "Measurement Update - Laser" << endl;
    
    // Set H matrix to default H matrix for laser
    ekf_.H_ = H_laser_;
    
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
