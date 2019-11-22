#include "robot_kinematics.h"

#include <math.h>


RobotKinematics::RobotKinematics(std::size_t dof)
    : DOF(dof) {

    q_.resize(DOF);
}

void RobotKinematics::setDH(std::vector<std::array<double, 4>> dh) {
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
    std::vector<double> qd(DOF);

    Eigen::MatrixXd geometryJacobian = calculateGeometryJacobian();

    Eigen::VectorXd xdVector(6);
    for (std::size_t i = 0; i < 6; i++) {
        xdVector[i] = xd[i];
    }

    Eigen::VectorXd qdVector = calculateInverseDLS(geometryJacobian, sigma) * xdVector;

    for (std::size_t i = 0; i < DOF; i++) {
        qd[i] = qdVector[i];
    }

    return qd;
}

Eigen::MatrixXd RobotKinematics::transformationMatrix(std::size_t jointIndex) {
    double sinQ = sin(q_[jointIndex]);
    double cosQ = cos(q_[jointIndex]);
    double sinAlpha = sin(dh_[jointIndex][0]);
    double cosAlpha = cos(dh_[jointIndex][0]);
    double a = dh_[jointIndex][1];
    double d = dh_[jointIndex][2];

    Eigen::MatrixXd t(4, 4);

    t(0, 0) = cosQ;
    t(0, 1) = -sinQ * cosAlpha;
    t(0, 2) = sinQ * sinAlpha;
    t(0, 3) = a * cosQ;

    t(1, 0) = sinQ;
    t(1, 1) = cosQ * cosAlpha;
    t(1, 2) = -cosQ * sinAlpha;
    t(1, 3) = a * sinQ;

    t(2, 0) = 0;
    t(2, 1) = sinAlpha;
    t(2, 2) = cosAlpha;
    t(2, 3) = d;

    t(3, 0) = 0;
    t(3, 1) = 0;
    t(3, 2) = 0;
    t(3, 3) = 1;

    return t;
}

Eigen::MatrixXd RobotKinematics::transformationMatrix() {
    Eigen::MatrixXd t = transformationMatrix(0);
    for (std::size_t i = 1; i < DOF; i++) {
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
    Eigen::MatrixXd geometryJacobian(6, DOF);
    Eigen::MatrixXd t[DOF];

    for (std::size_t i = 0; i < DOF; i++) {
        if (i == 0) {
            t[i] = transformationMatrix(i);
        } else {
            t[i] = t[i - 1] * transformationMatrix(i);
        }
    }

    Eigen::Vector3d p;
    for (std::size_t i = 0; i < 3; i++) {
        p[i] = t[DOF - 1](i, 3);
    }

    for (std::size_t i = 0; i < DOF; i++) {
        Eigen::Vector3d matZ;
        Eigen::Vector3d matP;

        if (i == 0) {
            matZ = Eigen::Vector3d(0.0, 0.0, 1.0);
            matP = Eigen::Vector3d(0.0, 0.0, 0.0);
        } else {
            for (int j = 0; j < 3; j++) {
                matZ[j] = t[i - 1](j, 2);
                matP[j] = t[i - 1](j, 3);
            }
        }

        Eigen::Vector3d matZCross = matZ.cross(p - matP);
        for (int j = 0; j < 3; j++) {
            geometryJacobian(j, i) = matZCross[j];
            geometryJacobian(j + 3, i) = matZ[j];
        }
    }

    return geometryJacobian;
}
