#include "mainwindow.hpp"
#include "ui_mainwindow.h"

#include <QDebug>
#include <cmath>

#include "kinematics/robot_kinematics.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);

    RobotKinematics kinematics(7);

    std::vector<std::array<double, 4>> DH(7);

    DH[0][0] =  M_PI / 2.0;
    DH[0][1] =  0.0;
    DH[0][2] =  0.0;
    DH[0][3] = 0.0;
    DH[1][0] = -M_PI / 2.0;
    DH[1][1] =  0.0;
    DH[1][2] = -0.165;
    DH[1][3] = 0.0;
    DH[2][0] =  M_PI / 2.0;
    DH[2][1] =  0.0;
    DH[2][2] =  0.425;
    DH[2][3] = 0.0;
    DH[3][0] = -M_PI / 2.0;
    DH[3][1] =  0.0;
    DH[3][2] =  0.1705;
    DH[3][3] = 0.0;
    DH[4][0] =  M_PI / 2.0;
    DH[4][1] =  0.0;
    DH[4][2] =  0.376;
    DH[4][3] = 0.0;
    DH[5][0] = -M_PI / 2.0;
    DH[5][1] =  0.0;
    DH[5][2] = -0.135;
    DH[5][3] = 0.0;
    DH[6][0] =  0.0;
    DH[6][1] =  0.0;
    DH[6][2] = 0.482;
    DH[6][3] = 0.0;

    kinematics.setDH(DH);

    const double DEG2RAD = M_PI / 180.0;
    const double RAD2DEG = 180.0 / M_PI;

    std::vector<double> initQ = {
        0.0 * DEG2RAD,
        -20.0 * DEG2RAD,
        0.0 * DEG2RAD,
        -50.0 * DEG2RAD,
        0.0 * DEG2RAD,
        0.0 * DEG2RAD,
        0.0 * DEG2RAD
    };
    std::vector<double> targetX = {0.50855, 0.1295, 0.97869, 0, 0, 0};
    std::vector<double> currentQ = initQ;

    kinematics.setQ(initQ);
    for (int i = 0; i < 100; i++) {
        std::vector<double> error(6);
        std::vector<double> currentX = kinematics.forwardKinematics(currentQ);
        std::vector<double> inputXd(6);
        double kp = 100.0;

        for (std::size_t i = 0; i < 6; i++) {
            error[i] = targetX[i] - currentX[i];
            inputXd[i] = 0 + kp * error[i];
        }

        std::vector<double> qd = kinematics.inverseDifferentialKinematics(inputXd, 0.001);

        for (std::size_t i = 0; i < 7; i++) {
            currentQ[i] += qd[i] * 0.001;
        }
    }

    for (std::size_t i = 0; i < 7; i++) {
        qDebug() << i << ":" << currentQ[i] * RAD2DEG;
    }
}

MainWindow::~MainWindow() {
    delete ui;
}
