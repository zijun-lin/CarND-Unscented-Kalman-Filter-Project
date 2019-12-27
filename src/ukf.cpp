#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  is_initialized_ = false;
  use_laser_ = true;
  use_radar_ = true;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // set weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i = 1; i < 2*n_aug_+1; ++i) {
    weights_(i) = 0.5/(lambda_+n_aug_);
  }

  NIS_laser_ = MatrixXd(2, 2);
  NIS_radar_ = MatrixXd(3, 3);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
    if(!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER){
        return;
    }
    if(!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR){
        return;
    }

  if ( !is_initialized_ ) {
      // first measurement
      x_.fill(0.0);
      if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
          double px = meas_package.raw_measurements_(0);
          double py = meas_package.raw_measurements_(1);
          x_ << px, py, 0, 0, 0;
      } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
          double rho = meas_package.raw_measurements_(0);
          double phi = meas_package.raw_measurements_(1);
          double rho_dot = meas_package.raw_measurements_(2);
          x_ << 0, 0, rho, phi, rho_dot;
      }
      P_ = MatrixXd::Identity(n_x_, n_x_);

      time_us_ = meas_package.timestamp_;

      // done initializing, no need to predict or update
      is_initialized_ = true;

      return;
  }

  // Prediction
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

    // Update(Radar or Laser)
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    } else {
        std::cout << " Error Type! " << std::endl;
    }

//  std::cout << " x: " << x_ << std::endl;
//  std::cout << " P: " << P_ << std::endl;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

    /*** Augmentation Assignment ***/
    // create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    // create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

    // create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;// mean is zero
    // create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5, 5) = P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;
    // create square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    // create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug_; ++i) {
        Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }

    /*** Sigma Point Prediction ***/
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        double v    = Xsig_aug(2, i);
        double yaw  = Xsig_aug(3, i);
        double yawd = Xsig_aug(4, i);
        double nu_a = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);

        // delta in x state
        VectorXd delta_x = VectorXd(n_x_);
        if (std::abs(yawd) < 0.000001) {
            delta_x(0) = v*cos(yaw)*delta_t;
            delta_x(1) = v*sin(yaw)*delta_t;
        } else {
            delta_x(0) = v/yawd*( sin(yaw+yawd*delta_t) - sin(yaw));
            delta_x(1) = v/yawd*(-cos(yaw+yawd*delta_t) + cos(yaw));
        }
        delta_x(2) = 0;
        delta_x(3) = yawd*delta_t;
        delta_x(4) = 0;

        // noise
        VectorXd noise = VectorXd(n_x_);
        noise(0) = 0.5*delta_t*delta_t*cos(yaw)*nu_a;
        noise(1) = 0.5*delta_t*delta_t*sin(yaw)*nu_a;
        noise(2) = delta_t*nu_a;
        noise(3) = 0.5*delta_t*delta_t*nu_yawdd;
        noise(4) = delta_t*nu_yawdd;

        Xsig_pred_(0, i) = Xsig_aug(0, i) + delta_x(0) + noise(0);
        Xsig_pred_(1, i) = Xsig_aug(1, i) + delta_x(1) + noise(1);
        Xsig_pred_(2, i) = Xsig_aug(2, i) + delta_x(2) + noise(2);
        Xsig_pred_(3, i) = Xsig_aug(3, i) + delta_x(3) + noise(3);
        Xsig_pred_(4, i) = Xsig_aug(4, i) + delta_x(4) + noise(4);
    }

    /*** Predicted Mean and Covariance ***/
    // predict state mean
    x_.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }

    // predict state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        double loc_yaw = x_diff(3);
        if (loc_yaw > M_PI || loc_yaw < -M_PI) {
        double new_yaw = std::fmod(loc_yaw, 2.0*M_PI);
            loc_yaw = new_yaw < 0 ? new_yaw + 2.0*M_PI : new_yaw;
            x_diff(3) = loc_yaw;
        }
        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }
    /*** Predicted End ***/
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
    /*** Predict Radar Measurement ***/
    // set measurement dimension, laser can measure px and py
    long n_z = 2;
    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);

    // transform sigma points into measurement space
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        double px   = Xsig_pred_(0, i);
        double py   = Xsig_pred_(1, i);

        Zsig(0, i) = px;
        Zsig(1, i) = py;
    }

    // calculate mean predicted measurement
    z_pred.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        z_pred = z_pred + weights_(i)*Zsig.col(i);
    }

    // calculate innovation covariance matrix S
    S.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R.fill(0.0);
    R(0, 0) = std_laspx_ * std_laspx_;
    R(1, 1) = std_laspy_ * std_laspy_;
    S = S + R;

    /*** UKF Update ***/
    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    // calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        // residual
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        double loc_angle = x_diff(3);
        if (loc_angle > M_PI || loc_angle < -M_PI) {
            double new_angle = std::fmod(loc_angle, 2.0*M_PI);
            loc_angle = new_angle < 0 ? new_angle + 2.0*M_PI : new_angle;
            x_diff(3) = loc_angle;
        }

        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    // update state mean and covariance matrix
    // residual
    VectorXd z = meas_package.raw_measurements_;

    VectorXd z_diff = z - z_pred;

    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
    std::cout << "NIS_laser_: " << NIS_laser_ << std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

    /*** Predict Radar Measurement ***/
    // set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;
    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);

    // transform sigma points into measurement space
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        double px   = Xsig_pred_(0, i);
        double py   = Xsig_pred_(1, i);
        double v    = Xsig_pred_(2, i);
        double yaw  = Xsig_pred_(3, i);

        Zsig(0, i) = sqrt(px*px + py*py);
        Zsig(1, i) = atan2(py, px);
        Zsig(2, i) = (px*cos(yaw)*v + py*sin(yaw)*v) / sqrt(px*px + py*py);
    }

    // calculate mean predicted measurement
    z_pred.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        z_pred = z_pred + weights_(i)*Zsig.col(i);
    }

    // calculate innovation covariance matrix S
    S.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        double loc_angle = z_diff(1);
        if (loc_angle > M_PI || loc_angle < -M_PI) {
            double new_angle = std::fmod(loc_angle, 2.0*M_PI);
            loc_angle = new_angle < 0 ? new_angle + 2.0*M_PI : new_angle;
            z_diff(1) = loc_angle;
        }
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R.fill(0.0);
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radrd_ * std_radrd_;
    S = S + R;


    /*** UKF Update ***/
    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    // calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {
        // residual
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        double loc_angle = x_diff(3);
        if (loc_angle > M_PI || loc_angle < -M_PI) {
            double new_angle = std::fmod(loc_angle, 2.0*M_PI);
            loc_angle = new_angle < 0 ? new_angle + 2.0*M_PI : new_angle;
            x_diff(3) = loc_angle;
        }

        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        double loc_angle2 = z_diff(1);
        if (loc_angle2 > M_PI || loc_angle2 < -M_PI) {
            double new_angle2 = std::fmod(loc_angle2, 2.0*M_PI);
            loc_angle2 = new_angle2 < 0 ? new_angle2 + 2.0*M_PI : new_angle2;
            z_diff(1) = loc_angle2;
        }

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    // update state mean and covariance matrix
    // residual
    VectorXd z = meas_package.raw_measurements_;

    VectorXd z_diff = z - z_pred;
    // angle normalization
    double loc_angle = z_diff(1);
    if (loc_angle > M_PI || loc_angle < -M_PI) {
        double new_angle = std::fmod(loc_angle, 2.0*M_PI);
        loc_angle = new_angle < 0 ? new_angle+2.0*M_PI : new_angle;
        z_diff(1) = loc_angle;
    }

    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
    std::cout << "NIS_radar_: " << NIS_radar_ << std::endl;
}