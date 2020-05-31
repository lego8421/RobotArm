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

    std::vector<RobotKinematics::DhParameter> dh(6);

    dh[0].alpha = M_PI / 2.0;
    dh[0].a = 0.0;
    dh[0].d = 0.3165;
    dh[0].theta = 0.0;
    dh[1].alpha =  0.0;
    dh[1].a =  0.661;
    dh[1].d =  0.0;
    dh[1].theta =  0.0;
    dh[2].alpha =  M_PI / 2.0;
    dh[2].a =  0.0;
    dh[2].d =  0.021;
    dh[2].theta =  0.0;
    dh[3].alpha = -M_PI / 2.0;
    dh[3].a =  0.0;
    dh[3].d =  0.533;
    dh[3].theta =  0.0;
    dh[4].alpha =  M_PI / 2.0;
    dh[4].a =  0.0;
    dh[4].d = -0.123;
    dh[4].theta =  0.0;
    dh[5].alpha =  0.0;
    dh[5].a =  0.0;
    dh[5].d =  0.077;
    dh[5].theta =  0.0;

    kinematics.setDH(dh);

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

    std::vector<double> init_q = RobotKinematics::toRadian({
        0.0,
        60.0,
        30.0,
        0.0,
        0.0,
        0.0
    });
    std::vector<double> target_x = {0.50855, 0.1295, 0.97869, 0, 0, 0};
    std::vector<double> current_q = init_q;

    kinematics.setQ(init_q);
    while(true) {
        std::vector<double> current_x = kinematics.forwardKinematics(current_q);
        std::vector<double> error = kinematics.getPoseError(target_x, current_x);
        std::vector<double> input_xd(6);
        double kp = 100.0;

        for (std::size_t i = 0; i < 6; i++) {
            error[i] = target_x[i] - current_x[i];
            input_xd[i] = 0 + kp * error[i];
        }

        std::vector<double> qd = kinematics.inverseDifferentialKinematics(input_xd, 0.001);

        for (std::size_t i = 0; i < 6; i++) {
            current_q[i] += qd[i] * 0.001;
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

    qDebug() << "inverse" << target_x << "to" << RobotKinematics::toDegree(current_q);
}

MainWindow::~MainWindow() {
    delete ui;
}
