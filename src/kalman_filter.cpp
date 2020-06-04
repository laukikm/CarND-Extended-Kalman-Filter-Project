#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  
  //std::cout<<"F="<<F_<<std::endl;  
  x_=F_*x_;  
  P_=F_*P_*F_.transpose()+Q_;

//  std::cout<<"Predicted Estimate="<<x_;
  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  //std::cout<<"Inside the laser function"<<std::endl;
  MatrixXd K=P_*H_.transpose()*(H_*P_*H_.transpose()+R_).inverse();

  VectorXd innovation=(z-H_*x_);

  double cutoff=25*(R_(0,0)+R_(1,1));

  if(innovation.norm()>cutoff){
    return;
  }

  x_=x_+K*(z-H_*x_);
  MatrixXd I=MatrixXd::Identity(4,4);
  P_=(I-K*H_)*P_;

  //std::cout<<"Laser Update="<<K;


}

VectorXd getRadarFromState(VectorXd& x_state){
  double px=x_state(0);
  double py=x_state(1);
  double vx=x_state(2);
  double vy=x_state(3);

  double r=sqrt(px*px+py*py);
  double phi=atan2(py,px); //Ranges from -180 to +180;
  double r_dot=(px*vx+py*vy)/r;

  VectorXd radar(3);
  radar<<r,phi,r_dot;
  return radar;
}
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  
  MatrixXd K(4,3);
  MatrixXd Temp(3,3);
  Temp=H_*P_*H_.transpose()+R_;
  
  
  K=P_*H_.transpose()*Temp.inverse();

  std::cout<<"K="<<K<<std::endl;
  VectorXd innovation=(z-getRadarFromState(x_));
  
  if(innovation(1)>M_PI){
    innovation(1)-=2*M_PI;
  }
  else if(innovation(1)<-M_PI){
    innovation(1)+=2*M_PI;
  }


  double cutoff=15*(R_(0,0)+R_(1,1)+R_(2,2));

  if(innovation.norm()>cutoff){
    return;
  }


  x_=x_+K*innovation; 
  std::cout<<"Error Term="<<innovation<<std::endl;

  MatrixXd I=MatrixXd::Identity(4,4);
  P_=(I-K*H_)*P_;

}
