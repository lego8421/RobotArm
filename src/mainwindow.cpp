#include "mainwindow.hpp"
#include "ui_mainwindow.h"

#include <QDebug>
#include <cmath>

#include "kinematics/robot_kinematics.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);

    RobotKinematics kinematics(6);

    std::vector<std::array<double, 4>> DH(6);

    DH[0][0] =  M_PI / 2.0;
    DH[0][1] =  0.0;
    DH[0][2] =  0.3165;
    DH[0][3] =  0.0;
    DH[1][0] =  0.0;
    DH[1][1] =  0.661;
    DH[1][2] =  0.0;
    DH[1][3] =  0.0;
    DH[2][0] =  M_PI / 2.0;
    DH[2][1] =  0.0;
    DH[2][2] =  0.021;
    DH[2][3] =  0.0;
    DH[3][0] = -M_PI / 2.0;
    DH[3][1] =  0.0;
    DH[3][2] =  0.533;
    DH[3][3] =  0.0;
    DH[4][0] =  M_PI / 2.0;
    DH[4][1] =  0.0;
    DH[4][2] = -0.123;
    DH[4][3] =  0.0;
    DH[5][0] =  0.0;
    DH[5][1] =  0.0;
    DH[5][2] =  0.077;
    DH[5][3] =  0.0;

    kinematics.setDH(DH);

    std::vector<double> q = RobotKinematics::toRadian({
        30.2127,
        91.9934,
        -10.126,
        -179.995,
        -98.1331,
        -30.211
    });

    kinematics.setQ(q);

    std::vector<double> x = kinematics.forwardKinematics(q);

    qDebug() << "forward" << q << "to" << x;

    std::vector<double> initQ = RobotKinematics::toRadian({
        0.0,
        60.0,
        30.0,
        0.0,
        0.0,
        0.0
    });
    std::vector<double> targetX = {0.50855, 0.1295, 0.97869, 0, 0, 0};
    std::vector<double> currentQ = initQ;

    kinematics.setQ(initQ);
    while(true) {
        std::vector<double> currentX = kinematics.forwardKinematics(currentQ);
        std::vector<double> error = kinematics.getPoseError(targetX, currentX);
        std::vector<double> inputXd(6);
        double kp = 100.0;

        for (std::size_t i = 0; i < 6; i++) {
            error[i] = targetX[i] - currentX[i];
            inputXd[i] = 0 + kp * error[i];
        }

        std::vector<double> qd = kinematics.inverseDifferentialKinematics(inputXd, 0.001);

        for (std::size_t i = 0; i < 6; i++) {
            currentQ[i] += qd[i] * 0.001;
        }

        double rms = 0.0;
        for (std::size_t i = 0; i < error.size(); i++) {
            rms += pow(error[i], 2);
        }
        rms = sqrt(rms);
        if (rms < 0.02) {
            break;
        }
    }

    qDebug() << "inverse" << targetX << "to" << currentQ;
}

MainWindow::~MainWindow() {
    delete ui;
}
