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
 * @file PIPnode.cpp
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

#include <chrono>
#include <random>


#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <glog/logging.h>

#include "progressivePIP.hpp"
#include "PIPnode.hpp"
#include "CompositePIPnode.hpp"
#include "PIPlkData.hpp"

#define ERR_STATE (-999)
#define DBL_EPSILON std::numeric_limits<double>::min()
#define MATCH_STATE 1
#define GAP_X_STATE 2
#define GAP_Y_STATE 3
#define STOP_STATE 4
#define LEFT 0
#define RIGHT 1

using namespace bpp;

PIPnode::PIPnode(const progressivePIP *pPIP,
                 tshlib::VirtualNode *vnode,
                 bpp::Node *bnode){

    progressivePIP_ = pPIP;

    vnode_ = vnode; // tshlib node
    bnode_ = bnode; // bpp node

    nodeID_ = bnode->getId();

}

bool PIPnode::_isRootNode(){

    // true if PIPnode is the root node
    // (the last node to be aligned) otherwise false
    if(parent == nullptr){
        return true;
    }else{
        return false;
    }

}

bool PIPnode::_isTerminalNode(){

    // true if PIPnode is a leaf node
    // (a leaf node doesn't need to be aligned) otherwise false
    if(childL == nullptr || childR == nullptr){
        return true;
    }else{
        return false;
    }

}

void PIPnode::_getPrFromSubstitutionModel() {

    if (!bnode_->hasFather()) {
        // root PIPnode doesn't have Pr
    } else {

        prNode_.resize(progressivePIP_->numCatg_);

        for (int i = 0; i < progressivePIP_->rDist_->getNumberOfCategories(); i++) {
            // substitution/deletion probabilities with rate variation (gamma)
            // Pr = exp( branchLength * rateVariation * Q )
            double brlen = bnode_->getDistanceToFather();

            if(brlen<SMALL_DOUBLE){
                LOG(FATAL) << "\nERROR branch length too small";
            }

            prNode_.at(i) = progressivePIP_->substModel_->getPij_t(brlen * \
                            progressivePIP_->rDist_->getCategory(i));
        }
    }

}

void PIPnode::_computeLK_M(std::vector<bpp::ColMatrix<double> > &fvL,
                           std::vector<bpp::ColMatrix<double> > &fvR,
                           std::vector<bpp::ColMatrix<double> > &Fv_M_ij,
                           std::vector<double> &Fv_sigma_M_ij,
                           double &pr_match,
                           double &pr_match_full_path) {


    // COMPUTE THE MATCH LK

    // fvL: fv array of the left child
    // fvR: fv array of the right child
    // Fv_M_ij: result of Fv_M_ij = Pr_R * fv_R hadamard_prod Pr_L * fv_L
    // Fv_sigma_M_ij: result of Fv_M_ij dot pi
    // pr_match: match probability (stored for the next layer)
    // pr_match_full_path: full match probability (used at this layer)
    // it encompasses the probability of an insertion along
    // the whole path between the root and this node

    // number of discrete gamma categories
    int num_gamma_categories = progressivePIP_->numCatg_;

    pr_match = 0.0;
    pr_match_full_path = 0.0;
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

        Fv_M_ij[catg] = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        Fv_sigma_M_ij.at(catg) = fv0;

        // match probability with gamma
        double p = progressivePIP_->rDist_->getProbability((size_t) catg) * \
               iotasNode_[catg] * betasNode_[catg] * fv0;

        pr_match_full_path += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                              alphaNode_.at(catg) * fv0;

        pr_match += p;
    }

}

void PIPnode::_computeLK_X(std::vector<bpp::ColMatrix<double> > &fvL,
                           std::vector<bpp::ColMatrix<double> > &fvR,
                           std::vector<bpp::ColMatrix<double> > &Fv_X_ij,
                           std::vector<double> &Fv_sigma_X_ij,
                           double &pr_gapx,
                           double &pr_gapx_full_path) {

    // COMPUTE THE GAPX LK

    // fvL: fv array of the left child
    // fvR: fv array of the right child
    // Fv_X_ij: result of Fv_X_ij = Pr_R * fv_R hadamard_prod Pr_L * fv_L
    // Fv_sigma_X_ij: result of Fv_sigma_X_ij dot pi
    // pr_gapx: gapx probability (stored for the next layer)
    // pr_gapx_full_path: full gapx probability (used at this layer)
    // it encompasses the probability of an insertion along
    // the whole path between the root and this node

    // number of discrete gamma categories
    int numCatg = progressivePIP_->numCatg_;

    pr_gapx = 0.0;
    pr_gapx_full_path = 0.0;
    for (int catg = 0; catg < numCatg; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_X_ij[catg] = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        Fv_sigma_X_ij.at(catg) = fv0;

        // gapX probability with gamma
        double p0 = progressivePIP_->rDist_->getProbability((size_t) catg) * \
                    iotasNode_[catg] * betasNode_[catg] * fv0;

        pr_gapx_full_path += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                              alphaNode_.at(catg) * fv0;

        pr_gapx += p0;
    }

}

void PIPnode::_computeLK_Y(std::vector<bpp::ColMatrix<double> > &fvL,
                           std::vector<bpp::ColMatrix<double> > &fvR,
                           std::vector<bpp::ColMatrix<double> > &Fv_Y_ij,
                           std::vector<double> &Fv_sigma_Y_ij,
                           double &pr_gapy,
                           double &pr_gapy_full_path) {

    // COMPUTE THE GAPY LK

    // fvL: fv array of the left child
    // fvR: fv array of the right child
    // Fv_Y_ij: result of Fv_Y_ij = Pr_R * fv_R hadamard_prod Pr_L * fv_L
    // Fv_sigma_Y_ij: result of Fv_sigma_Y_ij dot pi
    // pr_gapy: gapx probability (stored for the next layer)
    // pr_gapy_full_path: full gapx probability (used at this layer)
    // it encompasses the probability of an insertion along
    // the whole path between the root and this node

    // number of discrete gamma categories
    int numCatg = progressivePIP_->numCatg_;

    pr_gapy = 0.0;
    pr_gapy_full_path = 0.0;
    for (int catg = 0; catg < numCatg; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        Fv_Y_ij[catg] = fv;

        // fv0 = pi * fv
        double fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        Fv_sigma_Y_ij.at(catg) = fv0;

        // gapY probability with gamma
        double p0 = progressivePIP_->rDist_->getProbability((size_t) catg) * \
                    iotasNode_[catg] * betasNode_[catg] * fv0;

        pr_gapy_full_path += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                              alphaNode_.at(catg) * fv0;

        pr_gapy += p0;

    }

}

std::vector<double> PIPnode::_computeLkEmptyNode(std::vector<bpp::ColMatrix<double> > &fvL,
                                                 std::vector<bpp::ColMatrix<double> > &fvR,
                                                 std::vector<bpp::ColMatrix<double> > &fv_empty_data,
                                                 std::vector<double> &fv_empty_sigma,
                                                 std::vector<double> &lk_emptyL,
                                                 std::vector<double> &lk_emptyR,
                                                 std::vector<double> &lk_empty){

    // number of discrete gamma categories
    int numCatg = progressivePIP_->numCatg_;

    double p0;
    double pL,pR;
    double fv0;

    // array of lk (for each gamma rate) of a single column full of gaps
    std::vector<double> pc0;
    pc0.resize(numCatg);

    for (int catg = 0; catg < numCatg; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(childL->prNode_.at(catg), fvL.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(childR->prNode_.at(catg), fvR.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        fv_empty_data.at(catg) = fv;

        // fv0 = pi * fv
        fv0 = MatrixBppUtils::dotProd(fv, progressivePIP_->pi_);

        fv_empty_sigma.at(catg) = fv0;

        if(_isRootNode()){ // root
            // lk at root node (beta = 1.0)
            p0 = iotasNode_.at(catg) * fv0;
        }else{ // internal node
            p0 = ( iotasNode_.at(catg) - \
                   iotasNode_.at(catg) * betasNode_.at(catg) + \
                   iotasNode_.at(catg) * betasNode_.at(catg) * fv0 );
        }

        pL = lk_emptyL.at(catg); // lk_empty DOWN left
        pR = lk_emptyR.at(catg); // lk_empty DOWN right

        pc0.at(catg) = etaNode_.at(catg) + \
                       alphaNode_.at(catg) * fv0 + \
                       pL + pR; // this lk_empty is used at this layer

        lk_empty.at(catg) = p0 + pL + pR; // here store the lk for the next layer (probability UP is not added here)
    }

    return pc0;

}

//void PIPnode::_DP2D(PIPmsa *msaL,
//                    PIPmsa *msaR,
//                    int h_compr,
//                    int w_compr,
//                    std::vector<vector<double> > &Log2DM,
//                    std::vector<double> &Log2DX,
//                    std::vector<double> &Log2DY,
//                    std::vector<vector<double> > &Log2DM_fp,
//                    std::vector<double> &Log2DX_fp,
//                    std::vector<double> &Log2DY_fp,
//                    std::vector<vector<vector<bpp::ColMatrix<double> > > > &Fv_M,
//                    std::vector<vector<bpp::ColMatrix<double> > > &Fv_X,
//                    std::vector<vector<bpp::ColMatrix<double> > > &Fv_Y,
//                    std::vector<vector<vector<double> > > &Fv_sigma_M,
//                    std::vector<vector<double> > &Fv_sigma_X,
//                    std::vector<vector<double> > &Fv_sigma_Y){
//
//    //***************************************************************************************
//    // 2D LK COMPUTATION
//    //***************************************************************************************
//
//    int i,j;
//
//    // computes the lk in the two subtrees
//    std::vector<double> &lk_down_L = msaL->log_lk_down_;
//    std::vector<double> &lk_down_R = msaR->log_lk_down_;
//
//    // MATCH2D
//    double pr_m;
//    double pr_m_fp;
//    for (i = 0; i < h_compr; i++) {
//        for (j = 0; j < w_compr; j++) {
//
//            _computeLK_M(msaL->fv_data_.at(i),
//                         msaR->fv_data_.at(j),
//                         Fv_M[i][j],
//                         Fv_sigma_M[i][j],
//                         pr_m,
//                         pr_m_fp);
//
//            Log2DM[i][j] = log(pr_m); // stored for the next layer
//            Log2DM_fp[i][j] = log(pr_m_fp); // used at this node
//        }
//    }
//    //***************************************************************************************
//    // GAPX2D
//    double pr_x;
//    double pr_x_fp;
//    for (i = 0; i < h_compr; i++) {
//
//        _computeLK_X(msaL->fv_data_.at(i),
//                     msaR->fv_empty_data_,
//                     Fv_X[i],
//                     Fv_sigma_X[i],
//                     pr_x,
//                     pr_x_fp);
//
//        Log2DX[i] = progressivePIPutils::add_lns(log(pr_x),lk_down_L.at(i)); // stored for the next layer
//        Log2DX_fp[i] = progressivePIPutils::add_lns(log(pr_x_fp),lk_down_L.at(i)); // used at this node
//    }
//    //***************************************************************************************
//    // GAPY2D
//    double pr_y;
//    double pr_y_fp;
//    for (j = 0; j < w_compr; j++) {
//
//        _computeLK_Y(msaL->fv_empty_data_,
//                     msaR->fv_data_.at(j),
//                     Fv_Y[j],
//                     Fv_sigma_Y[j],
//                     pr_y,
//                     pr_y_fp);
//
//        Log2DY[j] = progressivePIPutils::add_lns(log(pr_y),lk_down_R.at(j)); // stored for the next layer
//        Log2DY_fp[j] = progressivePIPutils::add_lns(log(pr_y_fp),lk_down_R.at(j)); // used at this node
//    }
//
//}

void PIPnode::_alignStateMatrices2D(PIPmsa *msaL,
                                    PIPmsa *msaR,
                                    LKdata &lkdata){

    //***************************************************************************************
    // LK COMPUTATION
    //***************************************************************************************

    int i,j;

    // computes the lk in the two subtrees
    std::vector<double> &lk_down_L = msaL->log_lk_down_;
    std::vector<double> &lk_down_R = msaR->log_lk_down_;

    // MATCH2D
    double pr_m;
    double pr_m_fp;
    for (i = 0; i < lkdata.h_compr_; i++) {
        for (j = 0; j < lkdata.w_compr_; j++) {

            _computeLK_M(msaL->fv_data_.at(i),
                         msaR->fv_data_.at(j),
                         lkdata.Fv_M[i][j],
                         lkdata.Fv_sigma_M[i][j],
                         pr_m,
                         pr_m_fp);

            lkdata.Log2DM[i][j] = log(pr_m); // stored for the next layer
            lkdata.Log2DM_fp[i][j] = log(pr_m_fp); // used at this node
        }
    }
    //***************************************************************************************
    // GAPX2D
    double pr_x;
    double pr_x_fp;
    for (i = 0; i < lkdata.h_compr_; i++) {

        _computeLK_X(msaL->fv_data_.at(i),
                     msaR->fv_empty_data_,
                     lkdata.Fv_X[i],
                     lkdata.Fv_sigma_X[i],
                     pr_x,
                     pr_x_fp);

        lkdata.Log2DX[i] = progressivePIPutils::add_lns(log(pr_x),lk_down_L.at(i)); // stored for the next layer
        lkdata.Log2DX_fp[i] = progressivePIPutils::add_lns(log(pr_x_fp),lk_down_L.at(i)); // used at this node
    }
    //***************************************************************************************
    // GAPY2D
    double pr_y;
    double pr_y_fp;
    for (j = 0; j < lkdata.w_compr_; j++) {

        _computeLK_Y(msaL->fv_empty_data_,
                     msaR->fv_data_.at(j),
                     lkdata.Fv_Y[j],
                     lkdata.Fv_sigma_Y[j],
                     pr_y,
                     pr_y_fp);

        lkdata.Log2DY[j] = progressivePIPutils::add_lns(log(pr_y),lk_down_R.at(j)); // stored for the next layer
        lkdata.Log2DY_fp[j] = progressivePIPutils::add_lns(log(pr_y_fp),lk_down_R.at(j)); // used at this node
    }





    //==== DEBUG ===============
//    std::cout<<"\n";
//
//    std::cout<<"M2D\n";
//    for(int ii=0;ii<lkdata.h_compr_;ii++){
//        for(int jj=0;jj<lkdata.w_compr_;jj++){
//            double lk;
//            if(std::isinf(lkdata.Log2DM[ii][jj])){
//                lk=-0.0;
//            }else{
//                lk=lkdata.Log2DM[ii][jj];
//            }
//            printf("%8.6lf ",lk);
//        }
//        std::cout<<"\n";
//    }
//
//    std::cout<<"\n";
//
//    std::cout<<"X2D\n";
//    for(int ii=0;ii<lkdata.h_compr_;ii++){
//        double lk;
//        if(std::isinf(lkdata.Log2DX[ii])){
//            lk=-0.0;
//        }else{
//            lk=lkdata.Log2DX[ii];
//        }
//        printf("%8.6lf ",lk);
//    }
//    std::cout<<"\n\n";
//
//    std::cout<<"Y2D\n";
//    for(int jj=0;jj<lkdata.w_compr_;jj++){
//        double lk;
//        if(std::isinf(lkdata.Log2DY[jj])){
//            lk=-0.0;
//        }else{
//            lk=lkdata.Log2DY[jj];
//        }
//        printf("%8.6lf ",lk);
//    }
//    std::cout<<"\n";
    //==== DEBUG ===============





}

// TODO:re-implement this method for nodeCPU
bool PIPnode::_index_of_max(double m,
                            double x,
                            double y,
                            double epsilon,
                            std::default_random_engine &generator,
                            std::uniform_real_distribution<double> &distribution,
                            int &index,
                            double &val) {

    // get the index and the max value among the three input values (m,x,y)

    // m:  match value
    // x:  gapx value
    // y:   gapy value
    // epsilon: small number for the comparison between to numbers
    // generator: random number generator (when two or three numbers have the same value)
    // distribution: uniform distribution
    // index: index of max (1: MATCH, 2: GAPX, 3: GAPY)
    // val: max value between the three (m,x,y)

    double random_number;

    if (std::isinf(m) & std::isinf(x) & std::isinf(y)){
        // if the three values are -inf than this cell is marked as
        // non-valid (STOP_STATE) and the max val is -inf
        index = int(STOP_STATE);
        val = -std::numeric_limits<double>::infinity();
        return true;
    }

    if (not(std::isinf(m)) & not(std::isinf(x)) & (fabs((m - x)) < epsilon)) {
        // m and x are both not -inf
        // they are identical (their difference is smaller than epsilon)
        x = m; // x is exactly equal to m
    }

    if (not(std::isinf(m)) & not(std::isinf(y)) & (fabs((m - y)) < epsilon)) {
        // y and m are both not -inf
        // they are identical (their difference is smaller than epsilon)
        y = m; // y is exactly equal to m
    }

    if (not(std::isinf(x)) & not(std::isinf(y)) & (fabs((x - y)) < epsilon)) {
        // y and x are both not -inf
        // they are identical (their difference is smaller than epsilon)
        y = x; // y is exactly equal to x
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
                // m and y are equal and have the same value,
                // the state is selected with a uniform random
                // distribution with 50% probability each
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
                LOG(FATAL) << "\nSomething went wrong during the comparison in function "
                              "pPIP::_index_of_max. Check call stack below.";
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
                // x and y are equal and have the same value,
                // the state is selected with a uniform random
                // distribution with 50% probability each
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
                LOG(FATAL) << "\nSomething went wrong during the comparison in function "
                              "pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    } else {

        double mx = x;
        if (mx > y) {
            //m or x
            random_number = distribution(generator);
            // m and x are equal and have the same value,
            // the state is selected with a uniform random
            // distribution with 50% probability each
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
                // m,x and y are equal and have the same value,
                // the state is selected with a uniform random
                // distribution with 1/3 probability each
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
                LOG(FATAL) << "\nSomething went wrong during the comparison in function "
                              "pPIP::_index_of_max. Check call stack below.";
                return false;
            }
        }
    }

}