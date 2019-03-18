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
 * @file FactoryPIPnodeCPU.hpp
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

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Node.h>
#include <random>
#include <Utree.hpp>
#include <glog/logging.h>

#include "Utilities.hpp"

#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"
#include "CompositePIPmsa.hpp"

#ifndef MINIJATI_FACTORYNODECPU_HPP
#define MINIJATI_FACTORYNODECPU_HPP

namespace bpp {

    class nodeCPU : public PIPnode {

    private:

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

        void DP3D_PIP_leaf();

        void DP3D_PIP_node();

        bool _index_of_max(double m,
                           double x,
                           double y,
                           double epsilon,
                           std::default_random_engine &generator,
                           std::uniform_real_distribution<double> &distribution,
                           int &index,
                           double &val);

        double max_of_three(double a, double b, double c, double epsilon);

        std::vector<double> computeLK_GapColumn_local(MSAcolumn_t &sL,
                                                      MSAcolumn_t &sR);

        std::vector<double> computeLK_GapColumn_local(int nodeID,
                                                      int sonLeftID,
                                                      int sonRightID,
                                                      std::vector<bpp::ColMatrix<double> > &fvL,
                                                      std::vector<bpp::ColMatrix<double> > &fvR,
                                                      std::vector<bpp::ColMatrix<double> > &Fv_gap);

//        double compute_lk_gap_down(MSAcolumn_t &s,
//                                   int catg);


        double computeLK_M_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn_t &sL,
                                 MSAcolumn_t &sR,
                                 int m,
                                 std::map<MSAcolumn_t, double> &lkM,
                                 std::vector<std::vector<double> > &lkM_pattern);

        double computeLK_X_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn_t &sL,
                                 MSAcolumn_t &col_gap_R,
                                 int m,
                                 std::map<MSAcolumn_t, double> &lkX,
                                 std::vector<std::vector<double> > &lkX_pattern);

        double computeLK_Y_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn_t &col_gap_L,
                                 MSAcolumn_t &sR,
                                 int m,
                                 std::map<MSAcolumn_t, double> &lkY,
                                 std::vector<std::vector<double> > &lkY_pattern);

//        bpp::ColMatrix<double> computeFVrec(MSAcolumn_t &s,
//                                            int &idx,
//                                            int catg);

        std::string createGapCol(int len);

        void set_indices_M(int &up_corner_i,
                           int &up_corner_j,
                           int &bot_corner_i,
                           int &bot_corner_j,
                           int level,
                           int h,
                           int w);

        void set_indices_X(int &up_corner_i,
                           int &up_corner_j,
                           int &bot_corner_i,
                           int &bot_corner_j,
                           int level,
                           int h,
                           int w);

        void set_indices_Y(int &up_corner_i,
                           int &up_corner_j,
                           int &bot_corner_i,
                           int &bot_corner_j,
                           int level,
                           int h,
                           int w);

        signed long get_indices_M(int nx,
                                  int ny,
                                  int up_corner_i,
                                  int up_corner_j,
                                  int bot_corner_i,
                                  int bot_corner_j,
                                  int m,
                                  int h,
                                  int w);

        signed long get_indices_X(int nx,
                                  int ny,
                                  int up_corner_i,
                                  int up_corner_j,
                                  int bot_corner_i,
                                  int bot_corner_j,
                                  int m,
                                  int h,
                                  int w);

        signed long get_indices_Y(int nx,
                                  int ny,
                                  int up_corner_i,
                                  int up_corner_j,
                                  int bot_corner_i,
                                  int bot_corner_j,
                                  int m,
                                  int h,
                                  int w);

        void set_indices_T(int &up_corner_i,
                           int &up_corner_j,
                           int &bot_corner_i,
                           int &bot_corner_j,
                           int level,
                           int h,
                           int w);

        void reset_corner(int &up_corner_i,
                          int &up_corner_j,
                          int &bot_corner_i,
                          int &bot_corner_j,
                          int h,
                          int w);

        int get_indices_T(int nx,
                          int ny,
                          int up_corner_i,
                          int up_corner_j,
                          int bot_corner_i,
                          int bot_corner_j,
                          int m,
                          int h,
                          int w);

        bool is_inside(int x0,
                       int y0,
                       int xf,
                       int yf,
                       int xt,
                       int yt);

        bool checkboundary(int up_corner_i,
                           int up_corner_j,
                           int bot_corner_i,
                           int bot_corner_j,
                           int h,
                           int w);




        //-----------------------------------------------------------
        bpp::ColMatrix<double> computeFVrec(MSAcolumn_t &s,
                                            int &idx,
                                            int catg);
        void allgaps(std::string &s,
                     int &idx,
                     bool &flag);

        double compute_lk_gap_down(MSAcolumn_t &s,
                                   int catg);

        double _compute_lk_down_rec(int idx,
                                    double lk);

        double _compute_lk_down(MSAcolumn_t &s,
                                int catg);

        std::vector<double> _compute_lk_down();

        std::vector<double> _compute_lk_down(int idx_sb);
        //-----------------------------------------------------------

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

//        nodeCPU *parent;
//        nodeCPU *childL;
//        nodeCPU *childR;

        //PIPmsaSingle *MSA_; //contains the MSA

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        nodeCPU(const progressivePIP *pPIP, tshlib::VirtualNode *vnode, bpp::Node *bnode) : PIPnode(pPIP, vnode,
                                                                                                    bnode) {
        }

        //virtual ~nodeCPU() = default;

        ~nodeCPU(){ delete MSA_; };

        void DP3D_PIP();

        void _computeAllFvEmptySigmaRec();
    };

}

#endif //MINIJATI_FACTORYNODECPU_HPP
