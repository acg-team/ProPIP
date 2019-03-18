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
 * @file pPIP.cpp
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

#include "CompositePIPnode.hpp"

using namespace bpp;

void CompositePIPnode::addPIPnode(PIPnode *pip_node){

    pip_nodes_.at(pip_node->_getId()) = pip_node; // put the PIPnode in the array at the position indicated
                                                  // by the PIPnode Id (same as bpp Node Id). The array
                                                  // of PIPnodes is pre-allocated in the constructor

};

CompositePIPnode::CompositePIPnode(int numNodes){

    pip_nodes_.resize(numNodes); // resize the array of PIPnodes

}

void CompositePIPnode::PIPnodeAlign() {

    // ALIGN UNDER THE PIP MODEL

    // get the number of PIPnodes
    size_t num_nodes = pip_nodes_.size();

    for (int k = 0; k < num_nodes; k++) {

        // traverses the list of nodes and aligns the MSAs on the left and right side
        // if nodeInterface is a leaf the resulting MSA is the sequence itself

        ApplicationTools::displayGauge(k, num_nodes);

        pip_nodes_[k]->DP3D_PIP(); // align the input sequences under the PIP model

    }

}