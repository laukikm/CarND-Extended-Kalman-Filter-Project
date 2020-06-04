#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
//#include "cmath"

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

  /**
   * TODO: Finish initializing the FusionEKF

   * TODO: Set the process and measurement noises
   */

  H_laser_<< 1,0,0,0,
             0,1,0,0;

  //Q=MatrixXd::Identity(4,4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  Tools tools;
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    double cov_coeff=10000;
    ekf_.P_=cov_coeff*MatrixXd::Identity(4,4);

    ekf_.P_(0,2)=cov_coeff;
    ekf_.P_(1,3)=cov_coeff;
    ekf_.P_(2,0)=cov_coeff;
    ekf_.P_(3,1)=cov_coeff;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.

      double r=measurement_pack.raw_measurements_(0);
      double phi=measurement_pack.raw_measurements_(1);
      //double r_dot=measurement_pack.raw_measurements_(2);
      double dt=measurement_pack.timestamp_/1000000.0;

      double x=r*cos(phi);
      double y=r*sin(phi);

      double vx=x/dt;
      double vy=y/dt;

      ekf_.x_ << x, y, 0, 0;
      
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      double x=measurement_pack.raw_measurements_(0);
      double y=measurement_pack.raw_measurements_(1);
      ekf_.x_<<x,y,1,1;
      

    }

    previous_timestamp_=measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  double dt=(measurement_pack.timestamp_-previous_timestamp_)/1000000.0;
  double noise_ax=9;
  double noise_ay=9;

  F=MatrixXd::Identity(4,4);
  F << 1,0,dt,0,
       0,1,0,dt,
       0,0,1,0,
       0,0,0,1;
  //std::cout<<"F initialized"<<std::endl;
  

  Q=MatrixXd::Identity(4,4);

  /*
  Q<<noise_ax*pow(dt,4)/2,0,noise_ax*pow(dt,3)/2,0;
     0,noise_ay*pow(dt,4)/2,0,noise_ay*pow(dt,3)/2;
     noise_ax*pow(dt,3)/2,0,noise_ax*pow(dt,2),0,0;
     noise_ay*pow(dt,3)/2,0,0,noise_ay*pow(dt,2);
*/
  //std::cout<<"F and Q successfully initialized"<<endl;
  Q(0,0)=noise_ax*pow(dt,4)/2;
  
  
  Q(0,2)=noise_ax*pow(dt,3)/2;
  Q(1,1)=noise_ay*pow(dt,4)/2;
  Q(1,3)=noise_ay*pow(dt,3)/2;

  Q(2,0)=noise_ax*pow(dt,3)/2;
  Q(2,2)=noise_ax*pow(dt,2);
  
  
  
  Q(3,1)=noise_ay*pow(dt,3)/2;
  Q(3,3)=noise_ay*pow(dt,2);

  ekf_.F_=F;
  ekf_.Q_=Q;


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
    ekf_.R_=R_radar_;
    ekf_.H_=tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    cout<<"Sensor=RADAR"<<endl;

  } else {
    // TODO: Laser updates
    cout<<"Sensor=Laser"<<endl;
    ekf_.H_=H_laser_;
    ekf_.R_=R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  previous_timestamp_=measurement_pack.timestamp_;
  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
