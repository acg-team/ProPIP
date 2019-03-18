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
 * @file FactoryPIPnodeRAM.cpp
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
#include "CompositePIPmsa.hpp"

using namespace bpp;

#define EARLY_STOP_THR 10

double nodeRAM::max_of_three(double m,
                             double x,
                             double y,
                             double epsilon) {

    // get max value among the three input values (m,x,y)

    if (fabs(m) < epsilon) { // if a cell has been initialized to 0 then its lk is -inf
        m = -std::numeric_limits<double>::infinity();
    }
    if (fabs(x) < epsilon) { // if a cell has been initialized to 0 then its lk is -inf
        x = -std::numeric_limits<double>::infinity();
    }
    if (fabs(y) < epsilon) { // if a cell has been initialized to 0 then its lk is -inf
        y = -std::numeric_limits<double>::infinity();
    }

    if (std::isinf(m) && std::isinf(x) && std::isinf(y)) {
        // if all the three values are -inf then the max value is also -inf
        return -std::numeric_limits<double>::infinity();
    }

    if (m > x) {
        if (m > y) {
            // m greater than x and y
            return m;
        }
        // m greater than x but y greater than m
        return y;
    } else {
        if (x > y) {
            // x greater than m and y
            return x;
        }
        // x greater than m but y greater than x
        return y;
    }

}

double nodeRAM::_computeLK_MXY(double log_phi_gamma,
                               double valM,
                               double valX,
                               double valY,
                               double log_pr) {

    // compute the lk at a given matrix entry extending the previous best
    // lk (valM,valX,valY) together with the actual lk value (log_pr) and
    // the marginal lk of an empty column

    return log_phi_gamma + log_pr + max_of_three(valM, valX, valY, DBL_EPSILON);
}

void nodeRAM::DP3D_PIP_leaf() {

    //*******************************************************************************
    // ALIGNS LEAVES
    //*******************************************************************************

    // get vnode Id
    int vnodeId = (int) vnode_->vnode_seqid;

    // get sequence name from vnodeId
    std::string seqname = progressivePIP_->sequences_->getSequencesNames().at(vnodeId);

    // associates the sequence name to the leaf node
    MSA_->getMSA()->_setSeqNameLeaf(seqname);

    // get sequence from sequence name
    const bpp::Sequence *sequence = &progressivePIP_->sequences_->getSequence(seqname);

    // creates a vector containing the sequence associated to the leaf node
    MSA_->getMSA()->_setMSAleaf(sequence);

    // compresses sequence at the leaves
    MSA_->getMSA()->_compressMSA(progressivePIP_->alphabet_);

    // set fv_empty
    MSA_->getMSA()->_setFVemptyLeaf(progressivePIP_->numCatg_,
                                    progressivePIP_->alphabet_);

    // set fv_sigma_empty = fv_empty dot pi
    MSA_->getMSA()->_setFVsigmaEmptyLeaf(progressivePIP_->numCatg_);

    // computes the indicator values (fv values) at the leaves
    MSA_->getMSA()->_setFVleaf(progressivePIP_->numCatg_,
                               progressivePIP_->alphabet_);

    MSA_->getMSA()->_setFVsigmaLeaf(progressivePIP_->numCatg_,
                                    progressivePIP_->pi_);

    // compute the lk of an empty column
    MSA_->getMSA()->_computeLkEmptyLeaf(progressivePIP_,
                                        iotasNode_,
                                        betasNode_);

    // computes the lk for all the characters at the leaf
    MSA_->getMSA()->_computeLkLeaf(progressivePIP_,
                                   iotasNode_,
                                   betasNode_);

    // sets the traceback path at the leaf
    MSA_->getMSA()->_setTracebackPathLeaf();

}

void nodeRAM::DP3D(LKdata &lkdata,
                   double log_phi_gamma,
                   double log_nu_gamma,
                   double &curr_best_score, // best likelihood value at this node
                   int &level_max_lk) {      // depth in M,X,Y with the highest lk value


    //***************************************************************************************
    // 3D DYNAMIC PROGRAMMING
    //***************************************************************************************

//    printf("\n");
//    for(int i=0;i<21;i++){
//        for(int j=0;j<21;j++){
//            printf("%18.16lf ",progressivePIP_->substModel_->getGenerator().operator()(i,j));
//        }
//        printf("\n");
//    }
//    printf("\n");
//    printf("\n");
//    printf("\n");
//    for(int i=0;i<21;i++){
//        for(int j=0;j<21;j++){
//            printf("%18.16lf ",prNode_.at(0).operator()(i,j));
//        }
//        printf("\n");
//    }
//    printf("\n");

    //***************************************************************************************
    // DP VARIABLES
    //***************************************************************************************
    int m_binary_this; // Level Index during computation / current
    int m_binary_prev; // Level Index during computation / old
    auto epsilon = DBL_EPSILON; // very small number
    double min_inf = -std::numeric_limits<double>::infinity(); // -inf
    //***************************************************************************************
    // TRACEBACK VARIABLES
    //***************************************************************************************
    bool flag_exit = false; // early stop condition flag
    double prev_best_score = min_inf; // previuous best value at this node
    int tr_index = (int) STOP_STATE; // traceback index: 1=MATCH, 2=GAPX, 3=GAPY
    double max_lk_val = min_inf; // best lk value
    //***************************************************************************************
    // RANDOM NUMBERS GENERATOR
    //***************************************************************************************
    std::default_random_engine generator(progressivePIP_->getSeed()); // jatiapp seed
    std::uniform_real_distribution<double> distribution(0.0, 1.0); // uniform distribution for the selection
    // of lks with the same value
    //***************************************************************************************
    // EARLY STOP VARIABLES
    //***************************************************************************************
    int counter_to_early_stop = 0; // current number of consecutive steps where the lk decreases
    int max_decrease_before_stop = EARLY_STOP_THR; // hardcoded to prevent early-stops
    //***************************************************************************************
    // WORKING VARIABLES
    //***************************************************************************************
    int i, j;
    int id1m, id2m;
    int id1x, id2y;
    //***************************************************************************************
    // GET SONS
    //***************************************************************************************
    std::vector<int> *map_compr_L = &(childL->MSA_->getMSA()->map_compressed_seqs_);
    std::vector<int> *map_compr_R = &(childR->MSA_->getMSA()->map_compressed_seqs_);
    //***************************************************************************************
    // For each slice of the 3D cube, compute the values of each cell
    lkdata.Log3DM[0][0][0] = log_phi_gamma;
    lkdata.Log3DX[0][0][0] = log_phi_gamma;
    lkdata.Log3DY[0][0][0] = log_phi_gamma;
    lkdata.TR[0][0][0] = STOP_STATE;
    for (int m = 1; m < lkdata.d_; m++) {

        // if lk doesn't increase anymore for K steps (EARLY_STOP_THR) exit
        if (flag_exit) {
            break;
        }

        //***********************************************************************************
        // alternate the two layers
        m_binary_this = m % 2;
        m_binary_prev = (m + 1) % 2;
        //***********************************************************************************

        //***********************************************************************************
        // delta phi to add
        log_phi_gamma = -log((long double) m) + log_nu_gamma;
        //***********************************************************************************

        //***************************************************************************
        // GAPX[i][0]
        j = 0;
        for (i = 1; i < lkdata.h_; i++) {
            id1x = map_compr_L->at(i - 1);

            lkdata.Log3DM[m_binary_this][i - 1][j] = min_inf;
            lkdata.Log3DY[m_binary_this][i - 1][j] = min_inf;

            lkdata.Log3DX[m_binary_this][i][j] = _computeLK_MXY(log_phi_gamma,
                                                                min_inf,
                                                                lkdata.Log3DX[m_binary_prev][i - 1][j],
                                                                min_inf,
                                                                lkdata.Log2DX_fp[id1x]);;

            lkdata.TR[m][i][j] = (int) GAP_X_STATE;
        }
        //***********************************************************************************
        // GAPY[0][j]
        i = 0;
        for (j = 1; j < lkdata.w_; j++) {
            id2y = map_compr_R->at(j - 1);

            lkdata.Log3DM[m_binary_this][i][j - 1] = min_inf;
            lkdata.Log3DX[m_binary_this][i][j - 1] = min_inf;
            lkdata.Log3DY[m_binary_this][i][j] = _computeLK_MXY(log_phi_gamma,
                                                                min_inf,
                                                                min_inf,
                                                                lkdata.Log3DY[m_binary_prev][i][j - 1],
                                                                lkdata.Log2DY_fp[id2y]);

            lkdata.TR[m][i][j] = (int) GAP_Y_STATE;
        }
        //***********************************************************************************
        for (i = 1; i < lkdata.h_; i++) {
            for (j = 1; j < lkdata.w_; j++) {
                //***************************************************************************
                // MATCH[i][j]
                id1m = map_compr_L->at(i - 1);
                id2m = map_compr_R->at(j - 1);

                lkdata.Log3DM[m_binary_this][i][j] = _computeLK_MXY(log_phi_gamma,
                                                                    lkdata.Log3DM[m_binary_prev][i - 1][j - 1],
                                                                    lkdata.Log3DX[m_binary_prev][i - 1][j - 1],
                                                                    lkdata.Log3DY[m_binary_prev][i - 1][j - 1],
                                                                    lkdata.Log2DM_fp[id1m][id2m]);
                //***************************************************************************
                // GAPX[i][j]
                id1x = map_compr_L->at(i - 1);

                lkdata.Log3DX[m_binary_this][i][j] = _computeLK_MXY(log_phi_gamma,
                                                                    lkdata.Log3DM[m_binary_prev][i - 1][j],
                                                                    lkdata.Log3DX[m_binary_prev][i - 1][j],
                                                                    lkdata.Log3DY[m_binary_prev][i - 1][j],
                                                                    lkdata.Log2DX_fp[id1x]);
                //***************************************************************************
                // GAPY[i][j]
                id2y = map_compr_R->at(j - 1);

                lkdata.Log3DY[m_binary_this][i][j] = _computeLK_MXY(log_phi_gamma,
                                                                    lkdata.Log3DM[m_binary_prev][i][j - 1],
                                                                    lkdata.Log3DX[m_binary_prev][i][j - 1],
                                                                    lkdata.Log3DY[m_binary_prev][i][j - 1],
                                                                    lkdata.Log2DY_fp[id2y]);
                //***************************************************************************
                // TR[i][j]
                // Find which matrix contains the best value of LK found until this point.
                _index_of_max(lkdata.Log3DM[m_binary_this][i][j],
                              lkdata.Log3DX[m_binary_this][i][j],
                              lkdata.Log3DY[m_binary_this][i][j],
                              epsilon,
                              generator,
                              distribution,
                              tr_index,
                              max_lk_val);

                // Store the index for the traceback
                lkdata.TR[m][i][j] = tr_index;

                // If we reached the corner of the 3D cube, then:
                if ((m >= (lkdata.h_ - 1)) && (m >= (lkdata.w_ - 1)) && (i == (lkdata.h_ - 1)) &&
                    (j == (lkdata.w_ - 1))) {
                    // the algorithm is filling the last column of 3D DP matrix where
                    // all the characters are in the MSA

                    if (tr_index == (int) STOP_STATE) {
                        LOG(FATAL) << "\nSomething went wrong in reading the TR value. "
                                      "TR is neither MATCH, nor GAPX, nor GAPY. ";
                    }

                    if (max_lk_val > curr_best_score) {
                        curr_best_score = max_lk_val;
                        level_max_lk = m;
                    }

                    //***********************************************************************
                    // early stop condition
                    if (curr_best_score < prev_best_score) {
                        prev_best_score = curr_best_score;
                        counter_to_early_stop++;
                        if (counter_to_early_stop > max_decrease_before_stop) {
                            // if for max_decrease_before_stop consecutive times
                            // the lk decrease then exit, the maximum lk has been reached
                            flag_exit = true;
                        }
                    } else {
                        counter_to_early_stop = 0;
                    }
                    //***************************************************************************

                }
            }
        }
    }

}

void nodeRAM::DP3D_PIP_node() {

    //*******************************************************************************
    // ALIGNS INTERNAL NODE
    //*******************************************************************************

    DVLOG(1) << "DP3D_PIP at node: " << bnode_->getName();

    //***************************************************************************************
    // DP VARIABLES
    //***************************************************************************************
    double min_inf = -std::numeric_limits<double>::infinity(); // -inf
    //***************************************************************************************
    // TRACEBACK VARIABLES
    //***************************************************************************************
    double curr_best_score = min_inf; // best likelihood value at this node
    int level_max_lk = INT_MIN; // depth in M,X,Y with the highest lk value
    //***************************************************************************************
    // GAMMA VARIABLES
    //***************************************************************************************
    // number of discrete gamma categories
    size_t numCatg = progressivePIP_->numCatg_;
    //***************************************************************************************
    // GET SONS
    //***************************************************************************************
    std::vector<int> *map_compr_L = &(childL->MSA_->getMSA()->map_compressed_seqs_);
    std::vector<int> *map_compr_R = &(childR->MSA_->getMSA()->map_compressed_seqs_);
    //***************************************************************************************
    // DP SIZES
    //***************************************************************************************
    // Compute dimensions of the 3D block at current internal node.
    int h = childL->MSA_->getMSA()->_getMSAlength() + 1; // dimension of the alignment on the left side
    int w = childR->MSA_->getMSA()->_getMSAlength() + 1; // dimension of the alignment on the right side
    int d = (h - 1) + (w - 1) + 1; // third dimension of the DP matrix
    int h_compr = childL->MSA_->getMSA()->_getCompressedMSAlength(); // dimension of the compressed alignment on the left side
    int w_compr = childR->MSA_->getMSA()->_getCompressedMSAlength(); // dimension of the compressed alignment on the right side
    //***************************************************************************************
    // WORKING VARIABLES
    //***************************************************************************************
    int i, j;
    //***************************************************************************************
    // MEMORY ALLOCATION
    //***************************************************************************************
    // Initialisation of the data structure
    LKdata lkdata(d, h, h_compr, w, w_compr, numCatg, true);
    //***************************************************************************************
    // LK COMPUTATION OF AN EMPTY COLUMNS (FULL OF GAPS)
    //***************************************************************************************
    // computes the lk of an empty column in the two subtrees
    MSA_->getMSA()->lk_empty_.resize(numCatg);
    MSA_->getMSA()->fv_empty_data_.resize(numCatg);
    MSA_->getMSA()->fv_empty_sigma_.resize(numCatg);

    std::vector<bpp::ColMatrix<double> > &fvL = childL->MSA_->getMSA()->fv_empty_data_;
    std::vector<bpp::ColMatrix<double> > &fvR = childR->MSA_->getMSA()->fv_empty_data_;

    std::vector<bpp::ColMatrix<double> > &fv_empty_data = MSA_->getMSA()->fv_empty_data_;
    std::vector<double> &fv_empty_sigma = MSA_->getMSA()->fv_empty_sigma_;

    std::vector<double> &lk_emptyL = childL->MSA_->getMSA()->lk_empty_;
    std::vector<double> &lk_emptyR = childR->MSA_->getMSA()->lk_empty_;

    std::vector<double> &lk_empty = MSA_->getMSA()->lk_empty_;

    std::vector<double> pc0 = _computeLkEmptyNode(fvL,
                                                  fvR,
                                                  fv_empty_data,
                                                  fv_empty_sigma,
                                                  lk_emptyL,
                                                  lk_emptyR,
                                                  lk_empty);
    //***************************************************************************************
    // COMPUTES LOG(PHI(0))
    //***************************************************************************************
    // marginal likelihood for all empty columns with rate variation (gamma distribution)
    // phi(m,pc0,r) depends on the MSA length m
    // marginal phi marginalized over gamma categories
    double nu_gamma = 0.0;
    double log_phi_gamma = 0.0;
    for (int catg = 0; catg < numCatg; catg++) {
        // log( P_gamma(r) * phi(0,pc0(r),r) ): marginal lk for all empty columns of an alignment of size 0

        nu_gamma += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                    progressivePIP_->nu_.at(catg);

        log_phi_gamma += progressivePIP_->rDist_->getProbability((size_t) catg) * \
                         (progressivePIP_->nu_.at(catg) * \
                         (pc0.at(catg) - 1));
    }

    double log_nu_gamma = log(nu_gamma);
    //***************************************************************************************
    // 2D LK COMPUTATION
    //***************************************************************************************
    PIPmsa *pipmsaL = childL->MSA_->getMSA();
    PIPmsa *pipmsaR = childR->MSA_->getMSA();

    _alignStateMatrices2D(pipmsaL,
                          pipmsaR,
                          lkdata);
    //***************************************************************************************
    // 3D DYNAMIC PROGRAMMING
    //***************************************************************************************
    DP3D(lkdata,
         log_phi_gamma,
         log_nu_gamma,
         curr_best_score, // best likelihood value at this node
         level_max_lk);   // depth in M,X,Y with the highest lk value
    //***************************************************************************************
    // STORE THE SCORE
    //***************************************************************************************
    // level (k position) in the DP matrix that contains the highest lk value
    MSA_->getMSA()->score_ = curr_best_score;
    //***************************************************************************************
    // TRACEBACK ALGORITHM
    //***************************************************************************************
    std::vector<vector<bpp::ColMatrix<double> > > fv_data_not_compressed;
    fv_data_not_compressed.resize(level_max_lk);

    std::vector<std::vector<double>> fv_sigma_not_compressed;
    fv_sigma_not_compressed.resize(level_max_lk);

    std::vector<double> lk_down_not_compressed;
    lk_down_not_compressed.resize(level_max_lk);

    // start backtracing the 3 matrices (MATCH, GAPX, GAPY)
    MSA_->getMSA()->traceback_path_.resize(level_max_lk);

    i = h - 1;
    j = w - 1;
    int idmL, idmR;
    int state;
    for (int lev = level_max_lk; lev > 0; lev--) {
        state = lkdata.TR[lev][i][j];
        switch (state) {
            case MATCH_STATE:

                idmL = map_compr_L->at(i - 1);
                idmR = map_compr_R->at(j - 1);

                fv_data_not_compressed.at(lev - 1) = lkdata.Fv_M[idmL][idmR];
                fv_sigma_not_compressed.at(lev - 1) = lkdata.Fv_sigma_M[idmL][idmR];
                lk_down_not_compressed.at(lev - 1) = lkdata.Log2DM[idmL][idmR];

                i = i - 1;
                j = j - 1;

                MSA_->getMSA()->traceback_path_.at(lev - 1) = (int) MATCH_STATE;

                break;
            case GAP_X_STATE:

                idmL = map_compr_L->at(i - 1);

                fv_data_not_compressed.at(lev - 1) = lkdata.Fv_X[idmL];
                fv_sigma_not_compressed.at(lev - 1) = lkdata.Fv_sigma_X[idmL];
                lk_down_not_compressed.at(lev - 1) = lkdata.Log2DX[idmL];

                i = i - 1;

                MSA_->getMSA()->traceback_path_.at(lev - 1) = (int) GAP_X_STATE;

                break;
            case GAP_Y_STATE:

                idmR = map_compr_R->at(j - 1);

                fv_data_not_compressed.at(lev - 1) = lkdata.Fv_Y[idmR];
                fv_sigma_not_compressed.at(lev - 1) = lkdata.Fv_sigma_Y[idmR];
                lk_down_not_compressed.at(lev - 1) = lkdata.Log2DY[idmR];

                j = j - 1;

                MSA_->getMSA()->traceback_path_.at(lev - 1) = (int) GAP_Y_STATE;

                break;
            default:
                LOG(FATAL) << "\nSomething went wrong during the alignment reconstruction in function "
                              "pPIP::DP3D_PIP. Check call stack below.";
        }
    }
    //***************************************************************************************
    // BUILD NEW MSA
    //***************************************************************************************
    // converts traceback path into an MSA
    MSA_t *msaL = childL->MSA_->getMSA()->_getMSA();
    MSA_t *msaR = childR->MSA_->getMSA()->_getMSA();
    MSA_->getMSA()->_build_MSA(*msaL, *msaR);

    // assigns the sequence names of the new alligned sequences to the current MSA
    std::vector<string> *seqNameL = &(childL->MSA_->getMSA()->seqNames_);
    std::vector<string> *seqNameR = &(childR->MSA_->getMSA()->seqNames_);
    MSA_->getMSA()->_setSeqNameNode(*seqNameL, *seqNameR);
    //***************************************************************************************
    // COMPRESS INFO
    //***************************************************************************************
    // compress the MSA
    MSA_->getMSA()->_compressMSA(progressivePIP_->alphabet_);

    // compress fv values and lk_down
    MSA_->getMSA()->_compressLK(lk_down_not_compressed);
    MSA_->getMSA()->_compressFv(fv_data_not_compressed);
    MSA_->getMSA()->_compressFvSigma(fv_sigma_not_compressed);
    //***************************************************************************************
    // FREE MEMORY
    //***************************************************************************************
    // free memory
    lkdata.freeMemory(true);


    delete childL;
    delete childR;


}

void nodeRAM::DP3D_PIP() {

    if (_isTerminalNode()) {
        // align leaf (prepare data)
        DP3D_PIP_leaf();
    } else {
        // align internal node
        DP3D_PIP_node();
    }

}