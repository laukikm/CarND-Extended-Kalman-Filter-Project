#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using std::max;
using std::min;
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
  
  //std::cout<<"F="<<F_<<std::endl;  
  x_=F_*x_;  
  P_=F_*P_*F_.transpose()+Q_;
  
}

VectorXd KalmanFilter::saturateInnovationsLaser(VectorXd innovation){
  
  double cutoff=100*(R_(0,0)+R_(1,1));

  if(innovation.norm()>cutoff){
    innovation(0)=-1000;
  }  
  return innovation;
}

VectorXd KalmanFilter::saturateInnovationsRadar(VectorXd innovation){
  double cutoff=15*(R_(0,0)+R_(1,1)+R_(2,2));

  if(innovation.norm()>cutoff){
    innovation(0)=-1000;
  }
  return innovation;
}

bool KalmanFilter::isValid(VectorXd &innovation){
  //return false;
  return innovation(0)!=-1000;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  //std::cout<<"Inside the laser function"<<std::endl;
  MatrixXd K=P_*H_.transpose()*(H_*P_*H_.transpose()+R_).inverse();

  VectorXd innovation=(z-H_*x_);

  cout<<"Laser innovation="<<innovation<<std::endl;

  innovation=saturateInnovationsLaser(innovation);
  if(!isValid(innovation)){return;}

  x_=x_+K*innovation;
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

  
  VectorXd innovation=(z-getRadarFromState(x_));
  
  if(innovation(1)>M_PI){
    innovation(1)-=2*M_PI;
  }
  else if(innovation(1)<-M_PI){
    innovation(1)+=2*M_PI;
  }


  std::cout<<"Radar innovation="<<innovation<<std::endl;

  innovation=saturateInnovationsRadar(innovation);

  if(!isValid(innovation)){return;}

  x_=x_+K*innovation; 
  

  MatrixXd I=MatrixXd::Identity(4,4);
  P_=(I-K*H_)*P_;

}

KalmanFilterIS::KalmanFilterIS (){

  Lasersigma=VectorXd(2); 
  Lasersigma<<500,
              500; //Initialize to a high value so that it explores initially

  Radarsigma=VectorXd(3);
  Radarsigma<<100,5,200;

  gamma=0.5;

  mu=0.85;
}

bool KalmanFilterIS::isValid(VectorXd& innovation){
  //cout<<"All innovations are saturated, so we can include everything.";
  return true;
}

VectorXd KalmanFilterIS::saturateInnovationsRadar(VectorXd innovation){
  cout<<"Radarsigma="<<Radarsigma<<endl;
  for(int i=0;i<3;i++){
    innovation(i)=max(-sqrt(Radarsigma(i)),min(innovation(i),sqrt(Radarsigma(i))));
    Radarsigma(i)=mu*Radarsigma(i)+gamma*(innovation(i))*innovation(i);
  }
  return innovation;
}

VectorXd KalmanFilterIS::saturateInnovationsLaser(VectorXd innovation){
  cout<<"Lasersigma"<<Lasersigma<<endl;
  for(int i=0;i<2;i++){
    innovation(i)=max(-sqrt(Lasersigma(i)),min(innovation(i),sqrt(Lasersigma(i))));
    Lasersigma(i)=mu*Lasersigma(i)+gamma*(innovation(i))*innovation(i);
  }
  return innovation;
}