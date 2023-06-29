#ifndef ALGOS_H
#define ALGOS_H

#include <QtGlobal>
#include <QVector>
#include <Eigen/Dense>
#include <cassert>
#include "qcustomplot.h"
#include <QQueue>

template<typename K, typename V>
using pr = std::pair<K, V>;

template<typename T>
using v = QVector<T>;

template<typename T>
using que = QQueue<T>;

using namespace Eigen;

void fillVector(v<qreal>& vec, int const x1, int const x2, auto const& f) {
    assert((x1 >= 0) && (x1 < x2) && (x2 <= vec.size()));
    for (int i = x1; i < x2; ++i) {
        vec[i] = f(i);
    }
}

void fillVector(VectorXd& vec, int const x1, int const x2, auto const& f) {
    assert((x1 >= 0) && (x1 < x2) && (x2 <= vec.size()));
    for (int i = x1; i < x2; ++i) {
        vec[i] = f(i);
    }
}

void fillMatrix(MatrixXd& mtr, int const x1, int const x2, int const y1, int const y2, auto const& f) {
    assert((x1 >= 0) && (y1 >= 0) && (x1 < x2) && (y1 < y2) && (x2 <= mtr.cols()) && (y2 <= mtr.rows()));
    for (int i = y1; i < y2; ++i) {
        for (int j = x1; j < x2; ++j) {
            mtr(i, j) = f(i, j);
        }
    }
}

pr<VectorXd, int> gaussNewtonPolyRegression(v<pr<qreal, qreal>> const& points, int const max_step, qreal const dlt, int poly_deg);

pr<VectorXd, int> powellDogLegRegression(v<pr<qreal, qreal>> const& points, int const max_step, qreal const dlt, int poly_deg);

static qreal BFGSStep(auto const& f, auto const& grd, VectorXd const& x, VectorXd const& step) {
    qreal lr = 1;
    while ((f(x + lr * step) > f(x) + 0.5 * lr * grd(x).dot(step))/*
           || (grd(x + lr * step).dot(step) < 0.9 * grd(x).dot(step))*/) {
        lr /= 2;
    }
    return lr;
}


v<QCPCurveData> BFGS(auto const& f, auto const& grd, VectorXd cur, const int max_step, const qreal dlt)
{
    auto norm = [](auto const& vec) {
        qreal result = 0;
        for (auto const& elem : vec) {
            result += elem * elem;
        }
        return std::sqrt(result);
    };
    int const n = cur.size();
    MatrixXd mtr = MatrixXd::Identity(n, n);
    VectorXd grad, step, next, diff, sy;
    v<QCPCurveData> way = {{-1, cur[0], cur[1]}};
    qreal lr;

    for (int i = 0; i < max_step; ++i) {
        grad = grd(cur);
        if (norm(grad) < dlt) {
            break;
        }
        step = (mtr * grad) * -1;
        lr = BFGSStep(f, grd, cur, step);
        next = cur + lr * step;
        diff = next - cur;
        sy = grd(next) - grad;
        if (sy.dot(diff) > 0) {
            VectorXd mtr_s = mtr * diff;
            qreal s_mtr_s = diff.dot(mtr_s);
            mtr += (sy * sy.transpose()) / sy.dot(diff) - (mtr_s * mtr_s.transpose()) / s_mtr_s;
        }
        cur = next;
        way.emplace_back(i, cur[0], cur[1]);
    }
    return way;
}

static VectorXd LBFGSStep(VectorXd const& grad, que<VectorXd> const& s_hist, que<VectorXd> const& y_hist, que<qreal> const& rho_hist, MatrixXd const& H0) {
    VectorXd q = grad;
    v<qreal> alpha(s_hist.size());
    for (int i = s_hist.size() - 1; i >= 0; --i) {
        alpha[i] = rho_hist[i] * s_hist[i].dot(q);
        fillVector(q, 0, q.size(), [&q, &alpha, &y_hist, &i](int j) { return q[i] - alpha[i] * y_hist[i][j]; });
    }
    VectorXd r = H0 * q;
    qreal beta;
    for (int i = 0; i < s_hist.size(); ++i) {
        beta = rho_hist[i] * y_hist[i].dot(r);
        fillVector(r, 0, r.size(), [&i, &alpha, &beta, &s_hist, &r](int j){ return r[j] + s_hist[i][j] * (alpha[i] - beta);});
    }
    return -r;
}

static MatrixXd updtH(MatrixXd const& H0, que<VectorXd> const& s_hist, que<VectorXd> const& y_hist) {
    qreal q = y_hist.back().dot(s_hist.back()) / (y_hist.back().dot(y_hist.back()) + 1.e-8);
    qreal r = 1.0 / (y_hist.back().dot(s_hist.back()) + 1.e-4);
    return H0 + (q * y_hist.back() * y_hist.back().transpose()) - (r * s_hist.back() * s_hist.back().transpose());
}


v<QCPCurveData> LBFGS(auto const& f, auto const& grd, VectorXd x0, const int max_step, const qreal dlt, int const m)
{
    auto norm = [](auto const& vec) {
        qreal result = 0;
        for (auto const& elem : vec) {
            result += elem * elem;
        }
        return std::sqrt(result);
    };
    int const n = x0.size();
    VectorXd x = x0;
    VectorXd grad = grd(x);
    MatrixXd H0 = MatrixXd::Identity(n, n);
    que<VectorXd> s_hist, y_hist;
    que<qreal> rho_hist;
    v<QCPCurveData> way = {{-1, x[0], x[1]}};
    for (int i = 0; i < max_step; ++i) {
        VectorXd p = LBFGSStep(grad, s_hist, y_hist, rho_hist, H0);
        VectorXd s = 0.1 * p;
        VectorXd x_new = x + s;
        VectorXd g_new = grd(x_new);
        if (norm(g_new) < dlt) {
            break;
        }
        VectorXd y = g_new - grad;
        qreal rho = 1. / (y.dot(s) + 1.e-4);
        if (rho > 0) {
            s_hist.push_back(s);
            y_hist.push_back(y);
            rho_hist.push_back(rho);
            if (s_hist.size() > m) {
                s_hist.pop_front();
                y_hist.pop_front();
                rho_hist.pop_front();
            }
        }
        H0 = updtH(H0, s_hist, y_hist);
        x = x_new;
        grad = g_new;
        way.emplace_back(i, x[0], x[1]);
    }
    return way;
}

#endif // ALGOS_H
