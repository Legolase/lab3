#include "mainwindow.h"
#include <set>
#include <tuple>
#include "rand.h"
#include <fstream>
#include <chrono>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    plot.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
    plot.axisRect()->setupFullAxesBox(true);
    plot.xAxis->setLabel("x");
    plot.yAxis->setLabel("y");

    start();

    plot.legend->setVisible(true);
    resize(700, 700);
    setCentralWidget(&plot);
}

void MainWindow::start()
{
    auto f = [](VectorXd const& x) {
        return 0.5 * (x[0] - 2) * (x[0] - 2) + (x[1] + 3) * (x[1] + 3);
    };

    auto grd = [](VectorXd const& x) -> VectorXd {
        VectorXd v(2);
        v << x[0] - 2, x[1] * 2 + 6;
        return v;
    };
    auto f_wrp = [&f](qreal x, qreal y) {
        VectorXd v(2);
        v << x, y;
        return f(v);
    };
    set_color_map({-5, -5}, {5, 5}, {500, 500}, f_wrp);
    VectorXd st(2);
    st << 0, 0;
    auto result = BFGS(f, grd, st, 100, 10.e-5);
//    auto result = LBFGS(f, grd, st, 100, 10.e-5, 10);
    qDebug() << "Steps: " << result.size();
    v<QCPCurveData> way = {
        {0, 0, 0},
        {0, 0.5, -0.9},
        {0, 1, -1.6},
        {0, 1.5, -2.30},
        {0, 2, -3}
    };
    make_way(way, "Way");
//    std::ifstream fin(".\\input.txt");
//    if (!fin) {
//        qDebug() << "File not found";
//    }
//    qreal x_min, x_max, delta, dlt, b, dff;
//    v<qreal> a;
//    int max_step, poly_deg, number_points, size;
//    fin >> x_min >> x_max >> delta >> max_step >> dlt >> poly_deg >> number_points >> dff >> size;
//    a.resize(size);
//    for (auto& elem : a) {
//        fin >> elem;
//    }
//    fin >> b;

//    auto f = [&a, &b](qreal const x) {
//        qreal res = 1;
//        for (auto const& elem : a) {
//            res *= x + elem;
//        }
//        return res / b;
////        return ((x + a[0]) * (x + a[1]) * (x + a[2]) * (x + a[3]) * (x + a[4]) * (x + a[5]) * (x + a[6]) * (x + a[7]) * (x + a[8]) * (x + a[9]))/b;
////        return ((x + 3) * (x + 2.8) * (x -3.3) * (x + 1.75) * (x - 0.1) * (x + 3.9) * (x - 2.76) * (x - 3.6) * (x - 1) * (x - 0.1))/2024;
////        return ((x + 4.7) * (x + 1.05) * (x + 1.4) * (x - 1.25) * (x - 2) * (x + 0.45) * (x + 4.16) * (x - 0.6) * (x + 3.7) * (x + 2.9))/(-554.);
//    };
//    v<pr<qreal, qreal>> points = get_points_by_func(number_points, f, x_min, x_max, delta);
//    set_points(points, "hehe");
//    set_func(f, "func");

//    auto result = gaussNewtonPolyRegression(points, max_step, dlt, poly_deg);
//    for (auto & elem: result.first) {
//        elem += random(-dff, dff, 2);
//    }
//    qDebug() << result.second;
//    auto f2 = [&result](qreal const x) {
//        qreal res = 0;
//        for (int i = 0; i < result.first.size(); ++i) {
//            res += std::pow(x, i) * result.first[i];
//        }
//        return res;
//    };
//    set_func(f2, "result");
//    auto db = qDebug();
//    for (auto const& elem : result.first) {
//        db << elem;
//    }
}

void MainWindow::set_color_map(QPointF const& left_bottom, QPointF const& right_top,
                               QSize const& resolution,
                               auto const& f)
{
    QCPColorMap* color_map = new QCPColorMap(plot.xAxis, plot.yAxis);
    plot.legend->removeItem(0);
    color_map->data()->setSize(resolution.width(), resolution.height());
    color_map->data()->setRange(QCPRange(left_bottom.x(), right_top.x()),
                                QCPRange(left_bottom.y(), right_top.y()));
    color_map->setTightBoundary(true);
    qreal dx, dy;

    // set heights
    QPoint cur{0, 0};
    for (cur.ry() = 0; cur.y() < resolution.height(); ++cur.ry()) {
        for (cur.rx() = 0; cur.x() < resolution.width(); ++cur.rx()) {
            color_map->data()->cellToCoord(cur.x(), cur.y(), &dx, &dy);
            color_map->data()->setCell(cur.x(), cur.y(), f(dx, dy));
        }
    }

    //set gradient
    QCPColorScale* color_scale = new QCPColorScale(&plot);
    color_map->setInterpolate(true);
    plot.plotLayout()->addElement(0, 1, color_scale);
    color_scale->setType(QCPAxis::atRight);
    color_scale->axis()->setLabel("MSE");



    color_map->setColorScale(color_scale);
    color_map->rescaleDataRange();

    QCPColorGradient grad(QCPColorGradient::gpCold);
    grad.setLevelCount(17);
    color_map->setGradient(grad);


    QCPMarginGroup *marginGroup = new QCPMarginGroup(&plot);
    plot.axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
    color_scale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
    color_map->rescaleAxes();
}

static QPen random_pen() {
    static int arr[] = {3, 7, 8, 9, 10, 11, 12, 5};
    static int cur = 1;
    cur %= 8;
    QPen pen = QPen(QBrush(static_cast<Qt::GlobalColor>(arr[cur++])), 1.5);
    return pen;
//    uchar main_colour = static_cast<uchar>(random(175, 255));
//    uchar sc1 = static_cast<uchar>(random(25, 150));
//    uchar sc2 = static_cast<uchar>(random(25, 150));
//    switch (random(0, 2)) {
//    case 0:
//        pen = QPen(QBrush(QColor(main_colour, sc1, sc2)), 1.5);
//        break;
//    case 1:
//        pen = QPen(QBrush(QColor(sc1, main_colour, sc2)), 1.5);
//        break;
//    case 2:
//        pen = QPen(QBrush(QColor(sc1, sc2, main_colour)), 1.5);
//        break;
//    }
//    return pen;
}

void MainWindow::make_way(v<QCPCurveData> const& way, QString const& name)
{
    QCPCurve *curve = new QCPCurve(plot.xAxis, plot.yAxis);
    curve->setPen(random_pen());
    curve->setName(name);
    curve->data()->set(way, true);
    plot.addGraph();
    plot.graph()->addData(way[0].key, way[0].value);
    plot.graph()->setName("Start");
    plot.graph()->setPen(QPen(QBrush(),0));

    QCPScatterStyle sc_begin(QCPScatterStyle::ssCross);
    sc_begin.setPen(QPen(QBrush(Qt::green),2));
    plot.graph()->setScatterStyle(sc_begin);

    plot.addGraph();
    plot.graph()->addData(way.back().key, way.back().value);
    plot.graph()->setName("End");
    plot.graph()->setPen(QPen(QBrush(),0));

    QCPScatterStyle sc_end(QCPScatterStyle::ssTriangle);
    sc_end.setPen(QPen(QBrush(Qt::red),2));
    plot.graph()->setScatterStyle(sc_end);
}

void MainWindow::set_points(const v<pr<qreal, qreal>> &points, QString const& name)
{
    plot.addGraph();
    plot.graph()->setScatterStyle(QCPScatterStyle::ssDisc);
    plot.graph()->setPen(Qt::NoPen);
    plot.graph()->setName(name);
    for (auto const& elem : points) {
        plot.graph()->addData(elem.first, elem.second);
    }
    plot.rescaleAxes();
    auto x_range = plot.xAxis->range(), y_range = plot.yAxis->range();
    qreal x_diff = x_range.upper - x_range.lower, y_diff = y_range.upper - y_range.lower;
    plot.xAxis->setRange(x_range.lower - 0.05 * x_diff, x_range.upper + 0.05 * x_diff);
    plot.yAxis->setRange(y_range.lower - 0.05 * y_diff, y_range.upper + 0.05 * y_diff);
}

void MainWindow::set_line(const qreal k, const qreal b, QString const& name)
{
    auto f = [&k, &b](qreal const x) { return k * x + b; };
    QPen pen = random_pen();
    qreal const x_min = plot.xAxis->range().lower;
    qreal const x_max = plot.xAxis->range().upper;
    plot.addGraph();
    plot.graph()->setPen(pen);
    plot.graph()->addData({x_min, x_max}, {f(x_min), f(x_max)});
    plot.graph()->setName(name);
}

void MainWindow::set_func(auto const& f, const QString &name)
{
    QPen pen = random_pen();
    qreal const x_min = plot.xAxis->range().lower;
    qreal const x_max = plot.xAxis->range().upper;
    qreal const step = (x_max - x_min)/500;
    plot.addGraph();
    plot.graph()->setPen(pen);
    for (qreal x = x_min; x < x_max; x += step) {
        plot.graph()->addData(x, f(x));
    }
    plot.graph()->setName(name);
}


v<pr<qreal, qreal>> get_points_by_line(int n, const qreal k, const qreal b, const qreal x_min, const qreal x_max, const qreal delta)
{
    v<pr<qreal, qreal>> points;
    qreal x_temp;
    while (n-- > 0) {
        x_temp = random(x_min, x_max, 5);
        points.emplace_back(x_temp, k * x_temp + b + random(-delta, delta, 5));
    }
    return points;
}

v<pr<qreal, qreal>> get_points_by_func(int n, auto const& f, const qreal x_min, const qreal x_max, const qreal delta)
{
    v<pr<qreal, qreal>> points;
    qreal x_temp;
    while (n-- > 0) {
        x_temp = random(x_min, x_max, 5);
        points.emplace_back(x_temp, f(x_temp) + random(-delta, delta, 5));
    }
    return points;
}
