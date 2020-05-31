#ifndef ROBOTKINEMATICS_H
#define ROBOTKINEMATICS_H

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <cmath>


class RobotKinematics {
public:
    static constexpr double DEG2RAD = M_PI / 180.0;
    static constexpr double RAD2DEG = 180.0 / M_PI;

    static std::vector<double> toDegree(const std::vector<double> q) {
        std::vector<double> result;
        for (auto rad : q) {
            result.push_back(rad * RAD2DEG);
        }
        return result;
    }

    static std::vector<double> toRadian(const std::vector<double> q) {
        std::vector<double> result;
        for (auto deg : q) {
            result.push_back(deg * DEG2RAD);
        }
        return result;
    }

    union DhParameter {
        struct {
            double a;
            double alpha;
            double d;
            double theta;
        };
        struct {
            double p[4];
        };
    };

public:
    RobotKinematics(std::size_t dof);

public:
    void setDH(std::vector<DhParameter> dh);

    void setQ(std::vector<double> q);

    std::vector<double> forwardKinematics(std::vector<double> q);
    std::vector<double> inverseDifferentialKinematics(std::vector<double> xd, double sigma);
    static std::vector<double> getPoseError(const std::vector<double> desired, const std::vector<double> measured);
    static std::vector<double> getOrientationError(const Eigen::MatrixXd desired, const Eigen::MatrixXd measured);
    static Eigen::MatrixXd rpyToRotation(const std::vector<double> rpy);

private:
    Eigen::MatrixXd transformationMatrix(std::size_t joint_index);
    Eigen::MatrixXd transformationMatrix();

    Eigen::MatrixXd calculateInverse(Eigen::MatrixXd matrix);
    Eigen::MatrixXd calculateInverseDLS(Eigen::MatrixXd matrix, double sigma);
    Eigen::MatrixXd calculateGeometryJacobian();

private:
    const std::size_t dof_;
    std::vector<DhParameter> dh_;

    std::vector<double> q_;
};

#endif // ROBOTKINEMATICS_H
