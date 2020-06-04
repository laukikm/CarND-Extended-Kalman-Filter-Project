#include "tools.h"
#include <iostream>
#include "cmath"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
	VectorXd RMSE=VectorXd(4);
	RMSE<<0,0,0,0;
	//std::cout<<"RMSE Initialized as:"<<RMSE<<std::endl;
	for(int i=0;i<estimations.size();i++){
		VectorXd error;
		error=estimations[i]-ground_truth[i];
		RMSE+=error.cwiseProduct(error);
	
	}

	RMSE=RMSE/estimations.size();
	RMSE(0)=sqrt(RMSE(0));
	RMSE(1)=sqrt(RMSE(1));
	RMSE(2)=sqrt(RMSE(2));
	RMSE(3)=sqrt(RMSE(3));
	std::cout<<"size="<<std::endl;
	return RMSE;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
	
  double px=x_state(0);
  double py=x_state(1);
  double vx=x_state(2);
  double vy=x_state(3);

  double r=sqrt(px*px+py*py);
  double phi=atan2(py,px); //Ranges from -180 to +180;
  
  double r_dot=(px*vx+py*vy)/r;
  

  MatrixXd H_=MatrixXd::Zero(3,4);

  double c2=r;
  double c1=r*r;
  double c3=r*r*r;

  H_ << (px/c2), (py/c2), 0, 0,
      -(py/c1), (px/c1), 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  return H_;
}
