#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <utility>
#include "qcustomplot.h"
#include "algos.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);

private:
    void start();

    void set_color_map(QPointF const& left_bottom, QPointF const& right_top,
                                   QSize const& resolution,
                       auto const& f);
    void make_way(v<QCPCurveData>const& way, QString const& name);
    void set_points(v<pr<qreal, qreal>> const& points, QString const& name);
    void set_line(qreal const k, qreal const b, QString const& name);
    void set_func(auto const& f, QString const& name);

    QCustomPlot plot;
};

v<pr<qreal, qreal>> get_points_by_line(int n, qreal const k, qreal const b, qreal const x_min, qreal const x_max, qreal const delta);
v<pr<qreal, qreal>> get_points_by_func(int n, auto const& f, qreal const x_min, qreal const x_max, qreal const delta);

#endif // MAINWINDOW_H
