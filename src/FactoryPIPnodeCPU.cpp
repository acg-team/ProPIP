/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of miniJATI
 *
 * miniJATI is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * miniJATI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with miniJATI. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file FactoryPIPnodeCPU.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 19 02 2018
 * @version 1.0.7
 * @maintainer Lorenzo Gatti
 * @email lg@lorenzogatti.me
 * @maintainer Massimo Maiolo
 * @email massimo.maiolo@zhaw.ch
 * @status Development
 *
 * @brief
 * @details
 * @pre
 * @bug
 * @warning
 *
 * @see For more information visit:
 */
#include <random>

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <glog/logging.h>

#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"

using namespace bpp;

bool nodeCPU::_index_of_max(double m,
                            double x,
                            double y,
                            double epsilon,
                            std::default_random_engine &generator,
                            std::uniform_real_distribution<double> &distribution,
                            int &index,
                            double &val) {

    double random_number;

    if (std::isinf(m) & std::isinf(x) & std::isinf(y)){
        LOG(FATAL)
                << "\nSomething went wrong during the comparison of m,x,y variables in function pPIP::_index_of_max. Check call stack below. ";
        return false;
    }

    if (not(std::isinf(m)) & not(std::isinf(x)) & (fabs((m - x)) < epsilon)) {
        x = m;
    }

    if (not(std::isinf(m)) & not(std::isinf(y)) & (fabs((m - y)) < epsilon)) {
        y = m;
    }

    if (not(std::isinf(x)) & not(std::isinf(y)) & (fabs((x - y)) < epsilon)) {
        y = x;
    }

    if (m > x) {
        if (m > y) {
            index = int(MATCH_STATE);
            val = m;
            return true;
        } else if (y > m) {
            index = int(GAP_Y_STATE);
            val = y;
            return true;
        } else {
            if (abs(m - y) < epsilon) {
                //m or y
                random_number = distribution(generator);
                if (random_number < (1.0 / 2.0)) {
                    index = int(MATCH_STATE);
                    val = m;
                    return true;
                } else {
                    index = int(GAP_Y_STATE);
                    val = y;
                    return true;
                }
            } else {
                LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    } else if (x > m) {
        if (x > y) {
            index = int(GAP_X_STATE);
            val = x;
            return true;
        } else if (y > x) {
            index = int(GAP_Y_STATE);
            val = y;
            return true;
        } else {
            if (abs(x - y) < epsilon) {
                //x or y
                random_number = distribution(generator);
                if (random_number < (1.0 / 2.0)) {
                    index = int(GAP_X_STATE);
                    val = x;
                    return true;
                } else {
                    index = int(GAP_Y_STATE);
                    val = y;
                    return true;
                }
            } else {
                LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    } else {

        double mx = x;
        if (mx > y) {
            //m or x
            random_number = distribution(generator);
            if (random_number < (1.0 / 2.0)) {
                index = int(MATCH_STATE);
                val = m;
                return true;
            } else {
                index = int(GAP_X_STATE);
                val = x;
                return true;
            }
        } else if (y > mx) {
            index = int(GAP_Y_STATE);
            val = y;
            return true;
        } else {
            if (abs(mx - y) < epsilon) {
                //m or x or y
                random_number = distribution(generator);
                if (random_number < (1.0 / 3.0)) {
                    index = int(MATCH_STATE);
                    val = m;
                    return true;
                } else if (random_number < (2.0 / 3.0)) {
                    index = int(GAP_X_STATE);
                    val = x;
                    return true;
                } else {
                    index = int(GAP_Y_STATE);
                    val = y;
                    return true;
                }
            } else {
                LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    }

}

double nodeCPU::max_of_three(double a, double b, double c, double epsilon) {

    if (fabs(a) < epsilon) {
        a = -std::numeric_limits<double>::infinity();
    }
    if (fabs(b) < epsilon) {
        b = -std::numeric_limits<double>::infinity();
    }
    if (fabs(c) < epsilon) {
        c = -std::numeric_limits<double>::infinity();
    }

    if (std::isinf(a) && std::isinf(b) && std::isinf(c)) {
        LOG(FATAL) << "\nSomething went wrong during the comparison in function pPIP::max_of_three. Check call stack below.";
    }

    if (a > b) {
        if (a > c) {
            return a;
        }
        return c;
    } else {
        if (b > c) {
            return b;
        }
        return c;
    }

}

void nodeCPU::set_indices_M(int &up_corner_i,
                            int &up_corner_j,
                            int &bot_corner_i,
                            int &bot_corner_j,
                            int level,
                            int h,
                            int w) {

    if (level == 0) {
        up_corner_i = 0;
        up_corner_j = 0;
        bot_corner_i = 0;
        bot_corner_j = 0;
    } else {
        up_corner_i = 1 + level - std::min(w - 1, level);
        up_corner_j = std::min(w - 1, level);
        bot_corner_i = std::min(h - 1, level);
        bot_corner_j = 1 + level - std::min(h - 1, level);
    }

}

void nodeCPU::set_indices_X(int &up_corner_i,
                            int &up_corner_j,
                            int &bot_corner_i,
                            int &bot_corner_j,
                            int level,
                            int h,
                            int w) {

    if (level == 0) {
        up_corner_i = 0;
        up_corner_j = 0;
        bot_corner_i = 0;
        bot_corner_j = 0;
    } else {
        up_corner_i = 1 + level - 1 - std::min(w - 1, level - 1);
        up_corner_j = std::min(w - 1, level - 1);
        bot_corner_i = std::min(h - 1, level);
        bot_corner_j = level - std::min(h - 1, level);
    }

}

void nodeCPU::set_indices_Y(int &up_corner_i,
                            int &up_corner_j,
                            int &bot_corner_i,
                            int &bot_corner_j,
                            int level,
                            int h,
                            int w) {

    if (level == 0) {
        up_corner_i = 0;
        up_corner_j = 0;
        bot_corner_i = 0;
        bot_corner_j = 0;
    } else {
        up_corner_i = level - std::min(w - 1, level);
        up_corner_j = std::min(w - 1, level);
        bot_corner_i = std::min(h - 1, level - 1);
        bot_corner_j = 1 + level - 1 - std::min(h - 1, level - 1);
    }

}

signed long nodeCPU::get_indices_M(int nx,
                                   int ny,
                                   int up_corner_i,
                                   int up_corner_j,
                                   int bot_corner_i,
                                   int bot_corner_j,
                                   int m,
                                   int h,
                                   int w) {

    signed long idx;

    set_indices_M(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, m, h, w);

    if (is_inside(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, nx, ny)) {

        int dx, sx;

        dx = nx - up_corner_i + 1;

        sx = ((dx + 1) * dx / 2) - 1;

        idx = sx + (ny - up_corner_j);
    } else {
        idx = ERR_STATE;
    }

    return idx;

}

signed long nodeCPU::get_indices_X(int nx,
                                   int ny,
                                   int up_corner_i,
                                   int up_corner_j,
                                   int bot_corner_i,
                                   int bot_corner_j,
                                   int m,
                                   int h,
                                   int w) {

    signed long idx;

    set_indices_X(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, m, h, w);

    if (is_inside(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, nx, ny)) {

        int dx, sx;

        dx = nx - up_corner_i + 1;

        sx = ((dx + 1) * dx / 2) - 1;

        idx = sx + (ny - up_corner_j);
    } else {
        idx = ERR_STATE;
    }

    return idx;

}

signed long nodeCPU::get_indices_Y(int nx,
                                   int ny,
                                   int up_corner_i,
                                   int up_corner_j,
                                   int bot_corner_i,
                                   int bot_corner_j,
                                   int m,
                                   int h,
                                   int w) {

    signed long idx;

    set_indices_Y(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, m, h, w);

    if (is_inside(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, nx, ny)) {

        int dx, sx;

        dx = nx - up_corner_i + 1;

        sx = ((dx + 1) * dx / 2) - 1;

        idx = sx + (ny - up_corner_j);
    } else {
        idx = ERR_STATE;
    }

    return idx;

}

void nodeCPU::set_indices_T(int &up_corner_i,
                            int &up_corner_j,
                            int &bot_corner_i,
                            int &bot_corner_j,
                            int level,
                            int h,
                            int w) {

    int up_corner_i_x;
    int up_corner_i_y;

    int up_corner_j_x;
    int up_corner_j_y;

    int bot_corner_i_x;
    int bot_corner_i_y;

    int bot_corner_j_x;
    int bot_corner_j_y;

    set_indices_X(up_corner_i_x, up_corner_j_x, bot_corner_i_x, bot_corner_j_x, level, h, w);

    set_indices_Y(up_corner_i_y, up_corner_j_y, bot_corner_i_y, bot_corner_j_y, level, h, w);

    int delta_i, delta_j;

    delta_i = bot_corner_i_x - up_corner_i_y;
    delta_j = up_corner_j_y - bot_corner_j_x;

    if (delta_i > delta_j) {
        up_corner_i = up_corner_i_y;
        up_corner_j = up_corner_j_y;
        bot_corner_i = up_corner_i_y + delta_i;
        bot_corner_j = up_corner_j_y - delta_i;
    } else {
        up_corner_i = bot_corner_i_x - delta_j;
        up_corner_j = bot_corner_j_x + delta_j;
        bot_corner_i = bot_corner_i_x;
        bot_corner_j = bot_corner_j_x;
    }

}

void nodeCPU::reset_corner(int &up_corner_i,
                           int &up_corner_j,
                           int &bot_corner_i,
                           int &bot_corner_j,
                           int h,
                           int w) {

    int delta;

    if (up_corner_j >= w) {
        delta = up_corner_j - w + 1;
        up_corner_j -= delta;
        up_corner_i += delta;
    }
    if (bot_corner_i >= h) {
        delta = bot_corner_i - h + 1;
        bot_corner_i -= delta;
        bot_corner_j += delta;
    }

}

int nodeCPU::get_indices_T(int nx,
                           int ny,
                           int up_corner_i,
                           int up_corner_j,
                           int bot_corner_i,
                           int bot_corner_j,
                           int m,
                           int h,
                           int w) {

    set_indices_T(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, m, h, w);

    reset_corner(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, h, w);

    int idx;
    int dx, sx;

    dx = nx - up_corner_i + 1;

    sx = ((dx + 1) * dx / 2) - 1;

    idx = sx + (ny - up_corner_j);

    return idx;

}

bool nodeCPU::is_inside(int x0, int y0, int xf, int yf, int xt, int yt) {

    if ((xt < x0) || (yt > y0) || (xt > xf) || (yt < yf)) {
        return false;
    }

    if ((y0 - yt) > (xt - x0)) {
        return false;
    }

    return true;
}

bool nodeCPU::checkboundary(int up_corner_i,
                            int up_corner_j,
                            int bot_corner_i,
                            int bot_corner_j,
                            int h,
                            int w) {

    if ((up_corner_i >= 0) & (up_corner_i < h) & \
       (up_corner_j >= 0) & (up_corner_j < w) & \
       (bot_corner_i >= 0) & (bot_corner_i < h) & \
       (bot_corner_j >= 0) & (bot_corner_j < w)) {
        return true;
    }

    return false;
}

std::string nodeCPU::createGapCol(int len) {

    // create an MSA column full of gaps
    std::string colMSA(len, '-');

    return colMSA;
}

double nodeCPU::computeLK_M_local(double NU,
                                  double valM,
                                  double valX,
                                  double valY,
                                  bpp::Node *node,
                                  MSAcolumn_t &sL,
                                  MSAcolumn_t &sR,
                                  int m,
                                  std::map<MSAcolumn_t, double> &lkM,
                                  std::vector< std::vector<double> > &lkM_pattern) {

    double log_pr;

  //  bpp::Node *sonLeft = childL->_getBnode();
//    int sonLeftID = sonLeft->getId();

//    bpp::Node *sonRight = childR->_getBnode();
  //  int sonRightID = sonRight->getId();

    // create left + right column
    MSAcolumn_t s;
    s.append(sL);
    s.append(sR);

    // If the columns has been already seen in the alignment, then do not compute the value, but get it from the map
    bool is_found;
    std::map<MSAcolumn_t, double>::iterator it;

    it = lkM.find(s);
    is_found = (it != lkM.end());

    if (!is_found) {
        // is the first time that it computes the lk of this column

        int idx;

        // number of discrete gamma categories
        int num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

        double pr = 0.0;
        for (int catg = 0; catg < num_gamma_categories; catg++) {

            // computes the recursive Felsenstein's peeling weight on the left subtree
            idx = 0;
            bpp::ColMatrix<double> fvL = dynamic_cast<nodeCPU *>(childL)->computeFVrec(sL, idx, catg);

            // computes the recursive Felsenstein's peeling weight on the right subtree
            idx = 0;
            bpp::ColMatrix<double> fvR = dynamic_cast<nodeCPU *>(childR)->computeFVrec(sR, idx, catg);

            // PrfvL = Pr_L * fv_L
            bpp::ColMatrix<double> PrfvL;
            bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL, PrfvL);

            // PrfvR = Pr_R * fv_R
            bpp::ColMatrix<double> PrfvR;
            bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR, PrfvR);

            // fv = PrfvL * PrfvR
            bpp::ColMatrix<double> fv;
            bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

            // fv0 = pi * fv
            double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

            // match probability with gamma
            double p = progressivePIP_->rDist_->getProbability((size_t) catg) * \
                   iotasNode_[catg] * \
                   betasNode_[catg] * \
                   fv0;

            // marginal lk, marginalized over N gamma discrete classes
            pr += p;
        }

        log_pr = log((long double) pr);

        lkM[s] = log_pr;

    } else {
        // the lk of a column like this has been already computed,
        // the lk value can be retrieved from the map
        log_pr = it->second;
    }

    return log(NU) - log((double) m) + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON);
}

double nodeCPU::computeLK_X_local(double NU,
                                  double valM,
                                  double valX,
                                  double valY,
                                  bpp::Node *node,
                                  MSAcolumn_t &sL,
                                  MSAcolumn_t &col_gap_R,
                                  int m,
                                  std::map<MSAcolumn_t, double> &lkX,
                                  std::vector< std::vector<double> > &lkX_pattern) {

    double log_pr;


//    bpp::Node *sonLeft = childL->_getBnode();
//    int sonLeftID = sonLeft->getId();

  //  bpp::Node *sonRight = childR->_getBnode();
  //  int sonRightID = sonRight->getId();

    // create left + right column
    MSAcolumn_t s;
    bool is_found;
    std::map<MSAcolumn_t, double>::iterator it;

    s.append(sL);
    s.append(col_gap_R);

    it = lkX.find(s);
    is_found = (it != lkX.end());

    if (!is_found) {
        // is the first time that it computes the lk of this column

        int idx;

        // number of discrete gamma categories
        int numCatg = progressivePIP_->numCatg_;

        double pL;
        double pr = 0.0;
        for (int catg = 0; catg < numCatg; catg++) {

            // computes the recursive Felsenstein's peeling weight on the left subtree
            idx = 0;
            bpp::ColMatrix<double> fvL = dynamic_cast<nodeCPU *>(childL)->computeFVrec(sL, idx, catg);

            // computes the recursive Felsenstein's peeling weight on the right subtree
            idx = 0;
            bpp::ColMatrix<double> fvR = dynamic_cast<nodeCPU *>(childR)->computeFVrec(col_gap_R, idx, catg);

            // PrfvL = Pr_L * fv_L
            bpp::ColMatrix<double> PrfvL;
            bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL, PrfvL);

            // PrfvR = Pr_R * fv_R
            bpp::ColMatrix<double> PrfvR;
            bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR, PrfvR);

            // fv = PrfvL * PrfvR
            bpp::ColMatrix<double> fv;
            bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

            // fv0 = pi * fv
            double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

            // gapX probability with gamma
            double p0 = progressivePIP_->rDist_->getProbability((size_t) catg) * \
                    iotasNode_[catg] * \
                    betasNode_[catg] * \
                    fv0;

            pL = dynamic_cast<nodeCPU *>(childR)->_compute_lk_down(sL, catg);

            pr += p0 + pL;
        }

        log_pr = log((long double) pr);

        lkX[s] = log_pr;

    } else {
        // the lk of a column like this has been already computed,
        // the lk value can be retrieved from the map
        log_pr = it->second;
    }

    return log(NU) - log((double) m) + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON);
}

std::vector<double> nodeCPU::computeLK_GapColumn_local(MSAcolumn_t &sL,
                                                       MSAcolumn_t &sR) {

    // number of discrete gamma categories
    int num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    // array of lk (for each gamma rate) of a single column full of gaps
    std::vector<double> pc0;
    pc0.resize(num_gamma_categories);

    //bpp::Node *sonLeft = childL->_getBnode();
//    int sonLeftID = sonLeft->getId();

   // bpp::Node *sonRight = childR->_getBnode();
  //  int sonRightID = sonRight->getId();

    for (int catg = 0; catg < num_gamma_categories; catg++) {

        int idx;

        // computes the recursive Felsenstein's peeling weight on the left subtree
        idx = 0;
        bpp::ColMatrix<double> fvL = dynamic_cast<nodeCPU *>(childL)->computeFVrec(sL, idx, catg);

        // computes the recursive Felsenstein's peeling weight on the right subtree
        idx = 0;
        bpp::ColMatrix<double> fvR = dynamic_cast<nodeCPU *>(childR)->computeFVrec(sR, idx, catg);

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL, PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR, PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        // lk at the actual node (considered as root node => beta = 1.0)
        double p0 = iotasNode_[catg] * fv0;

        double pL,pR;
        pL = dynamic_cast<nodeCPU *>(childL)->compute_lk_gap_down(sL, catg);
        pR = dynamic_cast<nodeCPU *>(childR)->compute_lk_gap_down(sR, catg);

        pc0.at(catg) = p0 + pL + pR;
    }

    return pc0;
}

std::vector<double> nodeCPU::computeLK_GapColumn_local(int nodeID,
                                                       int sonLeftID,
                                                       int sonRightID,
                                                       std::vector< bpp::ColMatrix<double> > &fvL,
                                                       std::vector< bpp::ColMatrix<double> > &fvR,
                                                       std::vector< bpp::ColMatrix<double> > &Fv_gap) {

    // number of discrete gamma categories
    int num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    double fv0;
    double p0;
    double pL = 0.0;
    double pR = 0.0;

    // array of lk (for each gamma rate) of a single column full of gaps
    std::vector<double> pc0;
    pc0.resize(num_gamma_categories);

    for (int catg = 0; catg < num_gamma_categories; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_gap.at(catg) = fv;

        // fv0 = pi * fv
        fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        // lk at the actual node (considered as root node => beta = 1.0)
        p0 = iotasNode_[catg] * fv0;

//sis        pL=childL->log_lk_empty_down_[catg];

//sis        pR=childR->log_lk_empty_down_[catg];

        pc0.at(catg) = p0 + pL + pR;
    }

    return pc0;
}

double nodeCPU::computeLK_Y_local(double NU,
                                  double valM,
                                  double valX,
                                  double valY,
                                  bpp::Node *node,
                                  MSAcolumn_t &col_gap_L,
                                  MSAcolumn_t &sR,
                                  int m,
                                  std::map<MSAcolumn_t, double> &lkY,
                                  std::vector< std::vector<double> > &lkY_pattern) {

    double log_pr;

//    bpp::Node *sonLeft = childL->_getBnode();
//    int sonLeftID = sonLeft->getId();
//
//    bpp::Node *sonRight = childR->_getBnode();
//    int sonRightID = sonRight->getId();

    // create left + right column
    MSAcolumn_t s;
    bool is_found;
    std::map<MSAcolumn_t, double>::iterator it;

    s.append(col_gap_L);
    s.append(sR);

    it = lkY.find(s);
    is_found= (it != lkY.end());

    if (!is_found) {
        // is the first time that it computes the lk of this column

        // number of discrete gamma categories
        int num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

        int idx;

        double pR;
        double pr = 0.0;
        for (int catg = 0; catg < num_gamma_categories; catg++) {

            // computes the recursive Felsenstein's peeling weight on the left subtree
            idx = 0;
            bpp::ColMatrix<double> fvL = dynamic_cast<nodeCPU *>(childL)->computeFVrec(col_gap_L, idx, catg);

            // computes the recursive Felsenstein's peeling weight on the right subtree
            idx = 0;
            bpp::ColMatrix<double> fvR = dynamic_cast<nodeCPU *>(childR)->computeFVrec(sR, idx, catg);

            // PrfvL = Pr_L * fv_L
            bpp::ColMatrix<double> PrfvL;
            bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL, PrfvL);

            // PrfvR = Pr_R * fv_R
            bpp::ColMatrix<double> PrfvR;
            bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR, PrfvR);

            // fv = PrfvL * PrfvR
            bpp::ColMatrix<double> fv;
            bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

            // fv0 = pi * fv
            double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

            // gapY probability with gamma
            double p0 = progressivePIP_->rDist_->getProbability((size_t) catg) * \
                        iotasNode_[catg] * \
                        betasNode_[catg] * \
                        fv0;

            pR = dynamic_cast<nodeCPU *>(childR)->_compute_lk_down(sR, catg);

            pr += p0 + pR;

        }

        log_pr = log((long double) pr);

        lkY[s] = log_pr;

    } else {
        // the lk of a column like this has been already computed,
        // the lk value can be retrieved from the map

        log_pr = it->second;
    }

    return log(NU) - log((double) m) + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON);
}

bpp::ColMatrix<double> nodeCPU::computeFVrec(MSAcolumn_t &s, int &idx, int catg) {

    bpp::ColMatrix<double> fv;

    if (bnode_->isLeaf()) {

        //aggiustare
         //fv = fv_observed(s, idx);

    } else {

//        bpp::Node *sonLeft = childL->_getBnode();
//        int sonLeftID = sonLeft->getId();
//
//        bpp::Node *sonRight = childR->_getBnode();
//        int sonRightID = sonRight->getId();

        // computes the recursive Felsenstein's peeling weight on the left subtree
        bpp::ColMatrix<double> fvL = dynamic_cast<nodeCPU *>(childL)->computeFVrec(s, idx, catg);

        // computes the recursive Felsenstein's peeling weight on the right subtree
        bpp::ColMatrix<double> fvR = dynamic_cast<nodeCPU *>(childR)->computeFVrec(s, idx, catg);

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL, PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR, PrfvR);

        // fv = PrfvL * PrfvR
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

    }

    return fv;
}

void nodeCPU::allgaps(std::string &s, int &idx, bool &flag) {

    // flag is true if all the leaves of the subtree rooted in node contain a gap

    if (bnode_->isLeaf()) {
        char ch = s[idx];

        idx++;

        if (ch != '-') {
            flag = false;
        }

    } else {

//        bpp::Node *sonLeft = childL->_getBnode();
//        int sonLeftID = sonLeft->getId();
//
//        bpp::Node *sonRight = childR->_getBnode();
//        int sonRightID = sonRight->getId();

        dynamic_cast<nodeCPU *>(childL)->allgaps(s, idx, flag);
        dynamic_cast<nodeCPU *>(childR)->allgaps(s, idx, flag);
    }

}

double nodeCPU::compute_lk_gap_down(MSAcolumn_t &s, int catg) {

    int idx;

    if (bnode_->isLeaf()) {

        idx = 0;
        bpp::ColMatrix<double> fv = computeFVrec(s, idx, catg);

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        // lk of non survival till the leaves
        // lk = iota(v,r) - iota(v,r)*beta(v,r) + iota(v,r)*beta(v,r)*fv
        double pr = iotasNode_[catg] - \
             iotasNode_[catg] * betasNode_[catg] + \
             iotasNode_[catg] * betasNode_[catg] * fv0;

        return pr;

    }

   // bpp::Node *sonLeft = childL->_getBnode();
  //  int sonLeftID = sonLeft->getId();

 //   bpp::Node *sonRight = childR->_getBnode();
//    int sonRightID = sonRight->getId();

    idx = 0;
    bpp::ColMatrix<double> fv = computeFVrec(s, idx, catg);

    // fv0 = pi * fv
    double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

    // lk of non survival till the leaves
    // lk = iota(v,r) - iota(v,r)*beta(v,r) + iota(v,r)*beta(v,r)*fv
    double pr = iotasNode_[catg] - \
         iotasNode_[catg] * betasNode_[catg] + \
         iotasNode_[catg] * betasNode_[catg] * fv0;


    bool flagL = true;
    bool flagR = true;
    idx = 0;
    dynamic_cast<nodeCPU *>(childL)->allgaps(s, idx, flagL);

    int ixx = idx;
    dynamic_cast<nodeCPU *>(childR)->allgaps(s, idx, flagR);

    int len;

    MSAcolumn_t sL;
    len = ixx;
    sL = s.substr(0, len);
    double pL = dynamic_cast<nodeCPU *>(childL)->compute_lk_gap_down(sL, catg);

    MSAcolumn_t sR;
    sR = s.substr(ixx);
    double pR = dynamic_cast<nodeCPU *>(childR)->compute_lk_gap_down(sR, catg);

    return pr + pL + pR;
}

double nodeCPU::_compute_lk_down_rec(int idx,double lk){

    /*

    int num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    for (int catg=0; catg<num_gamma_categories; catg++) {
        lk = lk + progressivePIP_->rDist_->getProbability((size_t)catg) * \
             iotasNode_.at(catg) * \
             betasNode_.at(catg) * \
             fv_sigma_.at(idx).at(catg);
    }

    if (!bnode_->isLeaf()) {

        int idx_tr = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->rev_map_compressed_seqs_.at(idx);

        int tr = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.at(idx_tr);

        if(tr == (int)GAP_X_STATE){

            bpp::Node *sonLeft = childL->_getBnode();
            int sonLeftID = sonLeft->getId();

//            idx = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_mapL_.at(idx_tr);

            //int sub_i;// = subMSAidx_.at(LEFT).at(position);

            idx = static_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->map_compressed_seqs_.at(idx);

            idx = static_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->map_compressed_seqs_.at(idx);

            lk = childL->_compute_lk_down_rec(idx,lk);

        }else if(tr == (int)GAP_Y_STATE) {

            bpp::Node *sonRight = childR->_getBnode();
            int sonRightID = sonRight->getId();

//            idx = static_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_mapR_.at(idx_tr);

            //lk = childL->_compute_lk_down_rec(idx,lk);int sub_i;// = subMSAidx_.at(RIGHT).at(position);

            idx = static_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->map_compressed_seqs_.at(idx);

            idx = static_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->map_compressed_seqs_.at(idx);

            lk = childR->_compute_lk_down_rec(idx,lk);

        }

    }

    return lk;

     */
}

void nodeCPU::_computeAllFvEmptySigmaRec(){

}

double nodeCPU::_compute_lk_down(MSAcolumn_t &s, int catg) {

    int idx;

    if (bnode_->isLeaf()) {

        idx = 0;
        bpp::ColMatrix<double> fv = computeFVrec(s, idx, catg);

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        double pr = iotasNode_[catg] * betasNode_[catg] * fv0;

        return pr;

    }

//    bpp::Node *sonLeft = childL->_getBnode();
//    int sonLeftID = sonLeft->getId();

  //  bpp::Node *sonRight = childR->_getBnode();
  //  int sonRightID = sonRight->getId();

    idx = 0;
    bpp::ColMatrix<double> fv = computeFVrec(s, idx, catg);

    // fv0 = pi * fv
    double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

    double pr = iotasNode_[catg] * betasNode_[catg] * fv0;

    bool flagL = true;
    bool flagR = true;
    idx = 0;
    dynamic_cast<nodeCPU *>(childL)->allgaps(s, idx, flagL);
    int ixx = idx;
    dynamic_cast<nodeCPU *>(childR)->allgaps(s, idx, flagR);

    int len;
    if (flagR) {
        std::string sL;
        len = ixx;
        sL = s.substr(0, len);
        return pr + dynamic_cast<nodeCPU *>(childL)->_compute_lk_down(sL, catg);
    }

    if (flagL) {
        std::string sR;
        sR = s.substr(ixx);
        return pr + dynamic_cast<nodeCPU *>(childR)->_compute_lk_down(sR, catg);
    }

    return pr;
}

std::vector<double> nodeCPU::_compute_lk_down(){

    std::vector<double> lk_down;

    int MSAlen = dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->rev_map_compressed_seqs_.size();

    lk_down.resize(MSAlen);

    for(int idx=0;idx<MSAlen;idx++){
        double lk = 0.0;
        lk_down.at(idx) = _compute_lk_down_rec(idx,lk);
    }

    return lk_down;
}

void nodeCPU::DP3D_PIP_leaf() {

}

void nodeCPU::DP3D_PIP_node() {

    // number of discrete gamma categories
    size_t num_gamma_categories = progressivePIP_->rDist_->getNumberOfCategories();

    int up_corner_i;
    int up_corner_j;
    int bot_corner_i;
    int bot_corner_j;
    int lw;
    int h, w;

    int tr;

  //  bpp::Node *sonLeft = childL->_getBnode();
    //int sonLeftID = sonLeft->getId();

//    bpp::Node *sonRight = childR->_getBnode();
    //int sonRightID = sonRight->getId();

    // Compute dimensions of the 3D block at current internal node.
    h = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->msa_.size() + 1; // dimension of the alignment on the left side
    w = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->msa_.size() + 1; // dimension of the alignment on the riht side

    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix

    // lk of a single empty column (full of gaps) with rate variation (gamma distribution)
    std::vector<double> pc0;

    // MSA columns
    MSAcolumn_t sLs; // left column
    MSAcolumn_t sRs; // right column
    MSAcolumn_t col_gap_Ls; // left column full of gaps
    MSAcolumn_t col_gap_Rs; //right column full of gaps

    int numLeavesLeft = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->seqNames_.size(); // number of leaves in the left sub-tree
    int numLeavesRight = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->seqNames_.size(); // number of leaves in the right sub-tree

    col_gap_Ls = createGapCol(numLeavesLeft); // create column of gaps for the left sub-tree
    col_gap_Rs = createGapCol(numLeavesRight); // create column of gaps for the right sub-tree

    std::default_random_engine generator(progressivePIP_->getSeed()); // jatiapp seed
    std::uniform_real_distribution<double> distribution(0.0, 1.0); // Uniform distribution for the selection of lks with the same value

    auto epsilon = DBL_EPSILON;

    std::vector< std::vector<double> > lkM_pattern;
    std::vector< std::vector<double> > lkX_pattern;
    std::vector< std::vector<double> > lkY_pattern;

    //***************************************************************************************
    //***************************************************************************************
    // compute the lk of a column full of gaps
    pc0 = computeLK_GapColumn_local(col_gap_Ls, col_gap_Rs);
    //***************************************************************************************
    //***************************************************************************************

    //=================================================================================================
    std::vector< std::vector<double> > LogM; // DP sparse matrix for MATCH case (only 2 layer are needed)
    std::vector< std::vector<double> > LogX; // DP sparse matrix for GAPX case (only 2 layer are needed)
    std::vector< std::vector<double> > LogY; // DP sparse matrix for GAPY case (only 2 layer are needed)
    std::vector< std::vector<int> > TR; // 3D traceback matrix
    LogM.resize(2);
    LogX.resize(2);
    LogY.resize(2);
    TR.resize(d);
    //=================================================================================================

    // val num of cells occupied in a layer
    int numcells = int((w * (h + 1)) / 2);

    //=================================================================================================
    // allocate memory for the 2 layers
    LogM[0].resize(numcells);
    LogM[1].resize(numcells);
    LogX[0].resize(numcells);
    LogX[1].resize(numcells);
    LogY[0].resize(numcells);
    LogY[1].resize(numcells);
    //=================================================================================================

    //============================================================
    // marginal likelihood for all empty columns with rate variation (gamma distribution)
    // phi(m,pc0,r) depends on the MSA length m

    // marginal phi marginalized over gamma categories
    double log_phi_gamma;
    double prev_log_phi_gamma; // to store old value

    //=================================================================================================
    std::vector< std::vector<double> > PHI;
    PHI.resize(d);
    for (int i = 0; i < d; i++) {
        PHI[i].resize(num_gamma_categories);
    }
    //=================================================================================================

    double PC0 = 0.0;
    double NU = 0.0;
    for (int catg = 0; catg < num_gamma_categories; catg++) {
        // log( P_gamma(r) * phi(0,pc0(r),r) ): marginal lk for all empty columns of an alignment of size 0
        //PHI[0][catg] = log(rDist_->getProbability((size_t)catg)) + (nu_.at(catg) * (pc0.at(catg) - 1.0));
        PC0 += progressivePIP_->rDist_->getProbability((size_t) catg) * pc0.at(catg);
//        NU += progressivePIP_->rDist_->getProbability((size_t) catg) * nu_.at(catg);
    }

    // computes the marginal phi marginalized over all the gamma categories
    log_phi_gamma = NU * (PC0 - 1);
    //log_phi_gamma = PHI[0][0];
    //for (int catg = 1; catg < num_gamma_categories; catg++) {
    //    log_phi_gamma=pPIPUtils::add_lns(log_phi_gamma,PHI[0][catg]);
    //}
    //============================================================

    LogM[0][0] = log_phi_gamma;
    LogX[0][0] = log_phi_gamma;
    LogY[0][0] = log_phi_gamma;

    //=================================================================================================
    TR[0].resize(1);
    //=================================================================================================

    TR[0][0] = STOP_STATE;

    double max_of_3 = -std::numeric_limits<double>::infinity();

    signed long level_max_lk = INT_MIN;
    double val;
    int m_binary_this;
    int m_binary_prev;

    double valM;
    double valX;
    double valY;

    signed long idx;

    int coordSeq_1;
    int coordSeq_2;
    int coordTriangle_this_i;
    int coordTriangle_this_j;
    int coordTriangle_prev_i;
    int coordTriangle_prev_j;

    double score = -std::numeric_limits<double>::infinity();

    int depth;

    //int last_d = d - 1;
    int size_tr, tr_up_i, tr_up_j, tr_down_i, tr_down_j;
    std::map<MSAcolumn_t, double> lkM;
    std::map<MSAcolumn_t, double> lkX;
    std::map<MSAcolumn_t, double> lkY;

    //============================================================
    // early stop condition flag
    bool flag_exit = false;
    int counter_to_early_stop;
    int max_decrease_before_stop = 10;
    double prev_lk = -std::numeric_limits<double>::infinity();

    // ============================================================
    // For each slice of the 3D cube, compute the values of each cell

    for (int m = 1; m < d; m++) {

        if (flag_exit) {
            break;
        }

        // alternate the two layers
        m_binary_this = m % 2;
        m_binary_prev = (m + 1) % 2;

        //***************************************************************************************
        //***************************************************************************************
        for (int catg = 0; catg < num_gamma_categories; catg++) {
            // computes the marginal phi(m,pc0(r),r) with gamma by multiplying the starting value
            // phi(0,pco(r),r) = log( P_gamma(r) * exp( nu(r) * (pc0(r)-1) ) ) with
            // 1/m * nu(r) at each new layer
//            PHI[m][catg] = PHI[m - 1][catg] - log((long double) m) + log((long double) nu_.at(catg));
        }

        // store old value
        prev_log_phi_gamma = log_phi_gamma;

        // computes the marginal phi marginalized over all the gamma categories
        log_phi_gamma = PHI[m][0];
        for (int catg = 1; catg < num_gamma_categories; catg++) {
            log_phi_gamma = progressivePIPutils::add_lns(log_phi_gamma, PHI[m][catg]);
        }
        //***************************************************************************************
        //***************************************************************************************

        //***************************************************************************************
        //***************************************************************************************
        // COMPUTES MATCH LK
        set_indices_M(up_corner_i,
                      up_corner_j,
                      bot_corner_i,
                      bot_corner_j,
                      m, h, w);

        if (checkboundary(up_corner_i,
                          up_corner_j,
                          bot_corner_i,
                          bot_corner_j,
                          h, w)) {

            lw = 0;
            for (int i = up_corner_i; i <= bot_corner_i; i++) {

                coordTriangle_this_i = i;
                coordSeq_1 = coordTriangle_this_i - 1;
                coordTriangle_prev_i = coordTriangle_this_i - 1;

                // get left MSA column
                sLs = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->msa_.at(coordSeq_1);

                for (int j = 0; j <= lw; j++) {

                    coordTriangle_this_j = up_corner_j - j;
                    coordSeq_2 = coordTriangle_this_j - 1;
                    coordTriangle_prev_j = coordTriangle_this_j - 1;

                    // get right MSA column
                    sRs = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->msa_.at(coordSeq_2);

                    idx = get_indices_M(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valM = LogM[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valM = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_X(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valX = LogX[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valX = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_Y(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valY = LogY[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valY = -std::numeric_limits<double>::infinity();
                    }

                    if (std::isinf(valM) && std::isinf(valX) && std::isinf(valY)) {
                        LOG(FATAL) << "\nSomething went wrong during the comparison of valM, valX, valY in function pPIP::DP3D_PIP. Check call stack below.";
                    }

                    val = computeLK_M_local(NU,
                                            valM,
                                            valX,
                                            valY,
                                            bnode_,
                                            sLs,
                                            sRs,
                                            m,
                                            lkM,
                                            lkM_pattern);

                    if (std::isinf(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is infinite. Check call stack below.";
                    }

                    if (std::isnan(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                    }

                    idx = get_indices_M(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);

                    LogM[m_binary_this][idx] = val;
                }
                lw++;
            }
        }
        //***************************************************************************************
        //***************************************************************************************
        // COMPUTES GAPX LK
        set_indices_X(up_corner_i,
                      up_corner_j,
                      bot_corner_i,
                      bot_corner_j,
                      m, h, w);
        tr_down_i = bot_corner_i;
        tr_down_j = bot_corner_j;
        if (checkboundary(up_corner_i,
                          up_corner_j,
                          bot_corner_i,
                          bot_corner_j,
                          h, w)) {

            lw = 0;
            for (int i = up_corner_i; i <= bot_corner_i; i++) {

                coordTriangle_this_i = i;
                coordTriangle_prev_i = coordTriangle_this_i - 1;
                coordSeq_1 = coordTriangle_this_i - 1;

                // get left MSA column
                sLs = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->msa_.at(coordSeq_1);

                for (int j = 0; j <= lw; j++) {

                    coordTriangle_this_j = up_corner_j - j;
                    coordTriangle_prev_j = coordTriangle_this_j;

                    idx = get_indices_M(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valM = LogM[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valM = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_X(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valX = LogX[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valX = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_Y(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valY = LogY[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valY = -std::numeric_limits<double>::infinity();
                    }

                    if (std::isinf(valM) && std::isinf(valX) && std::isinf(valY)) {
                        LOG(FATAL) << "\nSomething went wrong during the comparison of valM, valX, valY in function pPIP::DP3D_PIP. Check call stack below.";
                    }

                    val = computeLK_X_local(NU,
                                            valM,
                                            valX,
                                            valY,
                                            bnode_,
                                            sLs,
                                            col_gap_Rs,
                                            m,
                                            lkX,
                                            lkX_pattern);

                    if (std::isinf(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is infinite. Check call stack below.";

                    }

                    if (std::isnan(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";
                    }
                    idx = get_indices_X(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);

                    LogX[m_binary_this][idx] = val;
                }
                lw++;
            }

        }
        //***************************************************************************************
        //***************************************************************************************
        // COMPUTES GAPY LK
        set_indices_Y(up_corner_i,
                      up_corner_j,
                      bot_corner_i,
                      bot_corner_j,
                      m, h, w);
        tr_up_i = up_corner_i;
        tr_up_j = up_corner_j;
        if (checkboundary(up_corner_i,
                          up_corner_j,
                          bot_corner_i,
                          bot_corner_j,
                          h, w)) {

            lw = 0;
            for (int i = up_corner_i; i <= bot_corner_i; i++) {
                coordTriangle_this_i = i;
                coordTriangle_prev_i = coordTriangle_this_i;
                for (int j = 0; j <= lw; j++) {

                    coordTriangle_this_j = up_corner_j - j;
                    coordTriangle_prev_j = coordTriangle_this_j - 1;
                    coordSeq_2 = coordTriangle_this_j - 1;

                    // get right MSA column
                    sRs = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->msa_.at(coordSeq_2);

                    idx = get_indices_M(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valM = LogM[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valM = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_X(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valX = LogX[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valX = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_Y(coordTriangle_prev_i,
                                        coordTriangle_prev_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m - 1, h, w);
                    if (idx >= 0) {
                        valY = LogY[m_binary_prev][idx];
                    } else {
                        // unreachable region of the 3D matrix
                        valY = -std::numeric_limits<double>::infinity();
                    }

                    if (std::isinf(valM) && std::isinf(valX) && std::isinf(valY)) {
                        LOG(FATAL) << "\nSomething went wrong during the comparison of valM, valX, valY in function pPIP::DP3D_PIP. Check call stack below.";
                    }

                    val = computeLK_Y_local(NU,
                                            valM,
                                            valX,
                                            valY,
                                            bnode_,
                                            col_gap_Ls,
                                            sRs,
                                            m,
                                            lkY,
                                            lkY_pattern);

                    if (std::isinf(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is infinite. Check call stack below.";
                    }

                    if (std::isnan(val)) {
                        LOG(FATAL) << "\nSomething went wrong function pPIP::DP3D_PIP. The value of 'val' is nan. Check call stack below.";

                    }

                    idx = get_indices_Y(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);

                    LogY[m_binary_this][idx] = val;
                }
                lw++;
            }

        }

        size_tr = (int) ceil((tr_down_i - tr_up_i + 1) * (tr_up_j - tr_down_j + 1 + 1) / 2);

        TR[m].resize(size_tr);

        set_indices_T(up_corner_i,
                      up_corner_j,
                      bot_corner_i,
                      bot_corner_j,
                      m, h, w);

        if (checkboundary(up_corner_i,
                          up_corner_j,
                          bot_corner_i,
                          bot_corner_j,
                          h, w)) {

            lw = 0;
            for (int i = up_corner_i; i <= bot_corner_i; i++) {
                coordTriangle_this_i = i;
                for (int j = 0; j <= lw; j++) {
                    coordTriangle_this_j = up_corner_j - j;

                    double mval;
                    double xval;
                    double yval;

                    idx = get_indices_M(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);
                    if (idx >= 0) {
                        mval = LogM[m_binary_this][idx];
                    } else {
                        mval = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_X(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);
                    if (idx >= 0) {
                        xval = LogX[m_binary_this][idx];
                    } else {
                        xval = -std::numeric_limits<double>::infinity();
                    }

                    idx = get_indices_Y(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);
                    if (idx >= 0) {
                        yval = LogY[m_binary_this][idx];
                    } else {
                        yval = -std::numeric_limits<double>::infinity();
                    }

                    // TODO:: remove these 3 lines
                    mval = fabs((long double) mval) < epsilon ? -std::numeric_limits<double>::infinity() : mval;
                    xval = fabs((long double) xval) < epsilon ? -std::numeric_limits<double>::infinity() : xval;
                    yval = fabs((long double) yval) < epsilon ? -std::numeric_limits<double>::infinity() : yval;

                    _index_of_max(mval, xval, yval, epsilon, generator, distribution, tr, max_of_3);

                    idx = get_indices_T(coordTriangle_this_i,
                                        coordTriangle_this_j,
                                        up_corner_i,
                                        up_corner_j,
                                        bot_corner_i,
                                        bot_corner_j,
                                        m, h, w);

                    if (TR[m][idx] != 0) {
                        LOG(FATAL) << "\nSomething went wrong in accessing TR at indices:[" << m << "][" << idx << "] in function pPIP::DP3D_PIP. Check call stack below.";
                    }

                    TR[m][idx] = tr;//max_val_index.index;

                    if ((coordTriangle_this_i == (h - 1)) & (coordTriangle_this_j == (w - 1))) {
                        // the algorithm is filling the last column of 3D DP matrix where
                        // all the characters are in the MSA

                        max_of_3 = max_of_three(mval, xval, yval, epsilon);

                        if (max_of_3 > score) {
                            score = max_of_3;
                            level_max_lk = m;
                        }

                        //=====================================================================
                        // early stop condition
                        if (score < prev_lk) {
                            prev_lk = score;
                            counter_to_early_stop++;
                            if (counter_to_early_stop > max_decrease_before_stop) {
                                // if for max_decrease_before_stop consecutive times
                                // the lk decrease then exit, the maximum lk has been reached
                                flag_exit = true;
                            }
                        } else {
                            counter_to_early_stop = 0;
                        }
                        //=====================================================================

                    }

                }
                lw++;
            }
        }
    }

    // level (k position) in the DP matrix that contains the highest lk value
    depth = level_max_lk;

    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->score_ = score;
    //==========================================================================================
    // start backtracing the 3 matrices (MATCH, GAPX, GAPY)
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.resize(depth);
    int id1 = h - 1;
    int id2 = w - 1;
    for (int lev = depth; lev > 0; lev--) {
        set_indices_T(up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, lev, h, w);
        idx = get_indices_T(id1, id2, up_corner_i, up_corner_j, bot_corner_i, bot_corner_j, lev, h, w);
//        int state = TR[lev][idx];
        switch (TR[lev][idx]) {
            case MATCH_STATE:
                id1 = id1 - 1;
                id2 = id2 - 1;
                dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.at(lev - 1) = (int)MATCH_STATE;
                break;
            case GAP_X_STATE:
                id1 = id1 - 1;
                dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.at(lev - 1) = (int)GAP_X_STATE;
                break;
            case GAP_Y_STATE:
                id2 = id2 - 1;
                dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->traceback_path_.at(lev - 1) = (int)GAP_Y_STATE;
                break;
            default:
                LOG(FATAL) << "\nSomething went wrong during the alignment reconstruction in function pPIP::DP3D_PIP. Check call stack below.";
        }
    }

    // converts traceback path into an MSA
    MSA_t *msaL = dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->_getMSA();
    MSA_t *msaR = dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->_getMSA();
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_build_MSA(*msaL,*msaR);

    // assigns the sequence names of the new alligned sequences to the current MSA
    std::vector<string> *seqNameL = &dynamic_cast<PIPmsaSingle *>(childL->MSA_)->pipmsa->seqNames_;
    std::vector<string> *seqNameR = &dynamic_cast<PIPmsaSingle *>(childR->MSA_)->pipmsa->seqNames_;
    dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa->_setSeqNameNode(*seqNameL,*seqNameR);
    //==========================================================================================

}

void nodeCPU::DP3D_PIP() {

    if (_isTerminalNode()) {
        DP3D_PIP_leaf();
    }else{
        DP3D_PIP_node();
    }

}


