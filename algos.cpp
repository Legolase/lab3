#include "algos.h"
#include "rand.h"
#include <QDebug>
#include <cmath>

#define DEBUG_OUTPUT 0

static void printMatrixXd(MatrixXd const& mtr) noexcept {
    for (int i = 0; i < mtr.rows(); ++i) {
        auto qdebug = qDebug();
        for (int j = 0; j < mtr.rows(); ++j) {
            qdebug << mtr(i, j);
        }
    }
}

qreal norm(auto const& vec) {
    qreal result = 0;
    for (auto const& elem : vec) {
        result += elem * elem;
    }
    return std::sqrt(result);
}

pr<VectorXd, int> gaussNewtonPolyRegression(const v<pr<qreal, qreal> > &points, const int max_step, const qreal dlt, int poly_deg)
{

    MatrixXd mtr(points.size(), poly_deg+1);
    fillMatrix(mtr, 0, 1, 0, mtr.rows(), [](int i, int j) { return 1; });
    fillMatrix(mtr, 1, mtr.cols(), 0, mtr.rows(), [&points](int i, int j) { return std::pow(points[i].first, j); });
    MatrixXd trans_mtr = mtr.transpose();
    VectorXd result(poly_deg+1);
    VectorXd diff(points.size());
    qreal scalar;
    MatrixXd tmp(poly_deg + 1, poly_deg + 1);
    VectorXd delta;

    for (int i = 0; i < max_step; ++i) {
        for (int j = 0; j < points.size(); ++j) {
            scalar = 0;
            for (int k = 0; k <= poly_deg; ++k) scalar += result[k] * mtr(j, k);
            diff[j] = points[j].second - scalar;
        }
        tmp = trans_mtr * mtr + MatrixXd::Identity(poly_deg + 1, poly_deg + 1);
        delta = tmp.colPivHouseholderQr().solve(trans_mtr * diff);
        result += delta;
#if DEBUG_OUTPUT
        qDebug() << norm(delta);
#endif
        if (norm(delta) < dlt) {
            return {std::move(result), i};
        }
    }
    return {std::move(result), max_step};
}

static VectorXd powellDogLegGrad(const v<pr<qreal, qreal> > &points, MatrixXd const& mtr, MatrixXd const& trans_mtr, VectorXd const& cf) {
    VectorXd pred = mtr * cf;
    VectorXd diff(points.size());
    fillVector(diff, 0, points.size(), [&points, &pred](int i){ return pred[i] - points[i].second; });
    return 2 * (trans_mtr * diff);
}

static VectorXd powellDogLegStep(VectorXd const& grad, MatrixXd const& hmtr, qreal const r) {
    VectorXd ustep = -grad / norm(grad);
    VectorXd bstep = -(hmtr.colPivHouseholderQr().solve(grad));
    if (norm(bstep) <= r) {
        return bstep;
    }
    qreal ub = ustep.dot(bstep);
    qreal uu = ustep.dot(ustep);
    qreal bb = bstep.dot(bstep);
    qreal res, temp = ub - uu + std::sqrt((ub - uu) * (ub - uu) + (r * r - uu) * bb);
    res = ((std::abs(temp) == 0) ? 0 : ((r * r - uu) / (ub - uu + std::sqrt((ub - uu) * (ub - uu) + (r * r - uu) * std::abs(bb)))));
    return res * ustep + (1 - res) * bstep;
}

pr<VectorXd, int> powellDogLegRegression(const v<pr<qreal, qreal> > &points, const int max_step, const qreal dlt, int poly_deg)
{
    MatrixXd mtr(points.size(), poly_deg+1);
    fillMatrix(mtr, 0, 1, 0, mtr.rows(), [](int i, int j) { return 1; });
    fillMatrix(mtr, 1, mtr.cols(), 0, mtr.rows(), [&points](int i, int j) { return std::pow(points[i].first, j); });
    MatrixXd trans_mtr = mtr.transpose();
    VectorXd gs(poly_deg+1), cf(poly_deg+1);
    qreal r = 1, fig;
    VectorXd grad = powellDogLegGrad(points, mtr, trans_mtr, cf);
    MatrixXd hmtr = 2 * (trans_mtr * mtr);
    for (int i = 1; i <= max_step; ++i) {
        if (norm(grad) < dlt) {
            return {std::move(cf), i};
        }
        VectorXd temp = powellDogLegStep(grad, hmtr, r);
        fig = std::pow(norm(temp), 2) / temp.dot(hmtr * temp);
        cf += fig * temp;
        grad = powellDogLegGrad(points, mtr, trans_mtr, cf);
        hmtr = 2 * (trans_mtr * mtr);
    }
    return {std::move(cf), max_step};
}
