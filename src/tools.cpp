#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE

    float pxy = px*px + py*py;
    float pxy_root = pow(pxy, 0.5);
    float pxy3 = pxy_root*pxy_root*pxy_root;

	//check division by zero
    if (pxy < 0.001)
    {
        cout << "Divide by zero";
        return Hj;
    }

	//compute the Jacobian matrix
    Hj(0, 0) = px/pxy_root;
    Hj(0, 1) = py/pxy_root;
    Hj(1, 0) = -py/pxy;
    Hj(1, 1) = px/pxy;
    Hj(2, 0) = py*(vx*py - vy*px)/pxy3;
    Hj(2, 1) = px*(vy*px - vx*py)/pxy3;
    Hj(2, 2) = px/pxy_root;
    Hj(2, 3) = py/pxy_root;

	return Hj;
}

VectorXd Tools::polarToCartesian(float rho, float phi)
{
	double x = rho*cos(phi);
	double y = rho*sin(phi);

	VectorXd xy = VectorXd(2);
	xy << x, y;

	return xy;
}

VectorXd Tools::cartesianToPolar(const VectorXd& c)
{
	double px = c(0);
	double py = c(1);
	double vx = c(2);
	double vy = c(3);

	double rho = sqrt(px*px + py*py);
	double phi = atan2(py, px);
	double rhodot = (px*vx + py*vy)/(sqrt(px*px + py*py));

	VectorXd polar = VectorXd(3);
	polar << rho, phi, rhodot;

	return polar;
}
