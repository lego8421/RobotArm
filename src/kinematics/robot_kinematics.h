#ifndef ROBOTKINEMATICS_H
#define ROBOTKINEMATICS_H

#include <Eigen/Dense>
#include <vector>
#include <array>


class RobotKinematics {
public:
    RobotKinematics(std::size_t dof);

public:
    void setDH(std::vector<std::array<double, 4>> dh);

    void setQ(std::vector<double> q);

    std::vector<double> forwardKinematics(std::vector<double> q);
    std::vector<double> inverseDifferentialKinematics(std::vector<double> xd, double sigma);


private:
    Eigen::MatrixXd transformationMatrix(std::size_t jointIndex);
    Eigen::MatrixXd transformationMatrix();

    Eigen::MatrixXd calculateInverse(Eigen::MatrixXd matrix);
    Eigen::MatrixXd calculateInverseDLS(Eigen::MatrixXd matrix, double sigma);
    Eigen::MatrixXd calculateGeometryJacobian();

private:
    const std::size_t DOF;
    std::vector<std::array<double, 4>> dh_;

    std::vector<double> q_;
};

#endif // ROBOTKINEMATICS_H
