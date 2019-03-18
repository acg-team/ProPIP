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
 * @file FactoryPIPnodeRAM.hpp
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

#ifndef MINIJATI_FACTORYNODERAM_HPP
#define MINIJATI_FACTORYNODERAM_HPP

namespace bpp {

    class nodeRAM : public PIPnode {

    private:

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

        // get max value among the three input values (m,x,y)
        double max_of_three(double m,           // match value
                            double x,           // gapx value
                            double y,           // gapy value
                            double epsilon);    // small number for the comparison between to numbers

        // compute the lk at a given matrix entry extending the previous best lk (valM,valX,valY) together with the
        // actual lk value (log_pr) and the marginal lk of an empty column
        double _computeLK_MXY(double log_phi_gamma,
                              double valM,
                              double valX,
                              double valY,
                              double log_pr);

        void DP3D(LKdata &lkdata,
                  double log_phi_gamma,
                  double log_nu_gamma,
                  double &curr_best_score, // best likelihood value at this node
                  int &level_max_lk); // depth in M,X,Y with the highest lk value

        void DP3D_PIP_leaf(); // DP method to align a sequence at a leaf PIPnode
                              // (which reduces to data preparation)

        void DP3D_PIP_node(); // DP method to align 2 MSAs at an internal node

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        // constructor
        nodeRAM(const progressivePIP *pPIP, tshlib::VirtualNode *vnode, bpp::Node *bnode) : PIPnode(pPIP,
                                                                                                    vnode,
                                                                                                    bnode) {
                // create a PIPmsaSingle object
                MSA_  = new PIPmsaSingle();

                // create a new PIPmsa
                dynamic_cast<PIPmsaSingle *>(MSA_)->pipmsa = new PIPmsa();

        }

        //virtual ~nodeRAM() = default;

        ~nodeRAM(){ delete MSA_; };

        void DP3D_PIP(); // DP algorithm to align (leaf/internal node) under the PIP model

    };

}

#endif //MINIJATI_FACTORYNODERAM_HPP
