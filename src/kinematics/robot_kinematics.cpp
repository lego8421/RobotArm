#include "robot_kinematics.h"

#include <math.h>


RobotKinematics::RobotKinematics(std::size_t dof)
    : dof_(dof) {

    q_.resize(dof_);
}

void RobotKinematics::setDH(std::vector<DhParameter> dh) {
    dh_ = dh;
}

void RobotKinematics::setQ(std::vector<double> q) {
    q_ = q;
}

std::vector<double> RobotKinematics::forwardKinematics(std::vector<double> q) {
    std::vector<double> x(6);

    setQ(q);

    Eigen::MatrixXd t = transformationMatrix();

    x[0] = t(0, 3);
    x[1] = t(1, 3);
    x[2] = t(2, 3);

    x[3] = atan2(t(2, 1), t(2, 2));
    x[4] = atan2(-t(2, 0), sqrt(t(2, 1) * t(2, 1) + t(2, 2) * t(2, 2)));
    x[5] = atan2(t(1, 0), t(0, 0));

    return x;
}

std::vector<double> RobotKinematics::inverseDifferentialKinematics(std::vector<double> xd, double sigma) {
    std::vector<double> qd(dof_);

    Eigen::MatrixXd geometry_jacobian = calculateGeometryJacobian();

    Eigen::VectorXd xd_vector(6);
    for (std::size_t i = 0; i < 6; i++) {
        xd_vector[i] = xd[i];
    }

    Eigen::VectorXd qd_vector = calculateInverseDLS(geometry_jacobian, sigma) * xd_vector;

    for (std::size_t i = 0; i < dof_; i++) {
        qd[i] = qd_vector[i];
    }

    return qd;
}

std::vector<double> RobotKinematics::getPoseError(const std::vector<double> desired, const std::vector<double> measured) {
    std::vector<double> error(6);

    // position error
    for(std::size_t i = 0; i < 3; i++) {
        error[i] = desired[i] - measured[i];
    }

    // orientation error
    Eigen::MatrixXd desired_rotation = rpyToRotation(std::vector<double>(desired.begin() + 3, desired.end()));
    Eigen::MatrixXd measured_rotation = rpyToRotation(std::vector<double>(measured.begin() + 3, measured.end()));
    std::vector<double> orientation_error = getOrientationError(desired_rotation, measured_rotation);

    for(std::size_t i = 0; i < 3; i++) {
        error[i + 3] = orientation_error[i];
    }

    return error;
}

std::vector<double> RobotKinematics::getOrientationError(const Eigen::MatrixXd desired, const Eigen::MatrixXd measured) {
    std::vector<double> error(3);
    double rotation[3][3];
    rotation[0][0] = desired(0, 0) * measured(0, 0) + desired(0, 1) * measured(0, 1) + desired(0, 2) * measured(0, 2);
    rotation[1][0] = desired(1, 0) * measured(0, 0) + desired(1, 1) * measured(0, 1) + desired(1, 2) * measured(0, 2);
    rotation[2][0] = desired(2, 0) * measured(0, 0) + desired(2, 1) * measured(0, 1) + desired(2, 2) * measured(0, 2);
    rotation[2][1] = desired(2, 0) * measured(1, 0) + desired(2, 1) * measured(1, 1) + desired(2, 2) * measured(1, 2);
    rotation[2][2] = desired(2, 0) * measured(2, 0) + desired(2, 1) * measured(2, 1) + desired(2, 2) * measured(2, 2);

    error[0] = atan2(rotation[2][1], rotation[2][2]);
    error[1] = atan2(-rotation[2][0], sqrt(rotation[2][1] * rotation[2][1] + rotation[2][2] * rotation[2][2]));
    error[2] = atan2(rotation[1][0], rotation[0][0]);

    return error;
}

Eigen::MatrixXd RobotKinematics::rpyToRotation(const std::vector<double> rpy) {
    Eigen::MatrixXd rotation(3, 3);

    rotation(0, 0) =  cos(rpy[2]) * cos(rpy[1]);
    rotation(0, 1) =  cos(rpy[2]) * sin(rpy[1]) * sin(rpy[0]) - sin(rpy[2]) * cos(rpy[0]);
    rotation(0, 2) =  cos(rpy[2]) * sin(rpy[1]) * cos(rpy[0]) + sin(rpy[2]) * sin(rpy[0]);

    rotation(1, 0) =  sin(rpy[2]) * cos(rpy[1]);
    rotation(1, 1) =  sin(rpy[2]) * sin(rpy[1]) * sin(rpy[0]) + cos(rpy[2]) * cos(rpy[0]);
    rotation(1, 2) =  sin(rpy[2]) * sin(rpy[1]) * cos(rpy[0]) - cos(rpy[2]) * sin(rpy[0]);

    rotation(2, 0) = -sin(rpy[1]);
    rotation(2, 1) =  cos(rpy[1]) * sin(rpy[0]);
    rotation(2, 2) =  0;

    return rotation;
}

Eigen::MatrixXd RobotKinematics::transformationMatrix(std::size_t jointIndex) {
    double sinQ = sin(q_[jointIndex]);
    double cosQ = cos(q_[jointIndex]);
    double sin_alpha = sin(dh_[jointIndex].alpha);
    double cos_alpha = cos(dh_[jointIndex].alpha);
    double a = dh_[jointIndex].a;
    double d = dh_[jointIndex].d;

    Eigen::MatrixXd t(4, 4);

    t(0, 0) = cosQ;
    t(0, 1) = -sinQ * cos_alpha;
    t(0, 2) = sinQ * sin_alpha;
    t(0, 3) = a * cosQ;

    t(1, 0) = sinQ;
    t(1, 1) = cosQ * cos_alpha;
    t(1, 2) = -cosQ * sin_alpha;
    t(1, 3) = a * sinQ;

    t(2, 0) = 0;
    t(2, 1) = sin_alpha;
    t(2, 2) = cos_alpha;
    t(2, 3) = d;

    t(3, 0) = 0;
    t(3, 1) = 0;
    t(3, 2) = 0;
    t(3, 3) = 1;

    return t;
}

Eigen::MatrixXd RobotKinematics::transformationMatrix() {
    Eigen::MatrixXd t = transformationMatrix(0);
    for (std::size_t i = 1; i < dof_; i++) {
        t = t * transformationMatrix(i);
    }
    return t;
}

Eigen::MatrixXd RobotKinematics::calculateInverse(Eigen::MatrixXd matrix) {
    Eigen::MatrixXd inverse(matrix.cols(), matrix.rows());
    if (matrix.rows() == matrix.cols()) {
        inverse = matrix.inverse();
    } else if (matrix.rows() > matrix.cols()) {
        // left
        inverse = (matrix.transpose() * matrix).inverse() * matrix.transpose();
    } else {
        // right
        inverse = matrix.transpose() * (matrix.transpose() * matrix).inverse();
    }
    return inverse;
}

Eigen::MatrixXd RobotKinematics::calculateInverseDLS(Eigen::MatrixXd matrix, double sigma) {
    if (matrix.cols() >= matrix.rows()) {
        Eigen::MatrixXd k = sigma * sigma * Eigen::MatrixXd::Identity(matrix.rows(), matrix.rows());
        Eigen::MatrixXd transpose = matrix.transpose();

        return transpose * calculateInverse(matrix * transpose + k);
    } else {
        Eigen::MatrixXd k = sigma * sigma * Eigen::MatrixXd::Identity(matrix.cols(), matrix.cols());
        return calculateInverse(matrix * matrix.transpose() + k) * matrix.transpose();
    }
}

Eigen::MatrixXd RobotKinematics::calculateGeometryJacobian() {
    Eigen::MatrixXd geometry_jacobian(6, dof_);
    Eigen::MatrixXd t[dof_];

    for (std::size_t i = 0; i < dof_; i++) {
        if (i == 0) {
            t[i] = transformationMatrix(i);
        } else {
            t[i] = t[i - 1] * transformationMatrix(i);
        }
    }

    Eigen::Vector3d p;
    for (std::size_t i = 0; i < 3; i++) {
        p[i] = t[dof_ - 1](i, 3);
    }

    for (std::size_t i = 0; i < dof_; i++) {
        Eigen::Vector3d mat_z;
        Eigen::Vector3d mat_p;

        if (i == 0) {
            mat_z = Eigen::Vector3d(0.0, 0.0, 1.0);
            mat_p = Eigen::Vector3d(0.0, 0.0, 0.0);
        } else {
            for (int j = 0; j < 3; j++) {
                mat_z[j] = t[i - 1](j, 2);
                mat_p[j] = t[i - 1](j, 3);
            }
        }

        Eigen::Vector3d mat_z_cross = mat_z.cross(p - mat_p);
        for (int j = 0; j < 3; j++) {
            geometry_jacobian(j, i) = mat_z_cross[j];
            geometry_jacobian(j + 3, i) = mat_z[j];
        }
    }

    return geometry_jacobian;
}
