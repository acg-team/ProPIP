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
 * @file pPIP.hpp
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

#include <string>

#include "PIPnode.hpp"

#ifndef MINIJATI_COMPOSITEPIPNODE_HPP
#define MINIJATI_COMPOSITEPIPNODE_HPP

namespace bpp {

    class CompositePIPnode {

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        std::vector<PIPnode *> pip_nodes_; // array of PIPnodes

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        CompositePIPnode(int numNodes); // constructor

        void addPIPnode(PIPnode *pip_node); // put the PIPnode in the array at the position indicated
                                            // by the PIPnode Id (same as bpp Node Id). The array
                                            // of PIPnodes is pre-allocated in the constructor

        void PIPnodeAlign(); // align all the input sequences under the PIP model

    };

}

#endif //MINIJATI_COMPOSITEPIPNODE_HPP
