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
 * @file FactoryPIPnode.hpp
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

#ifndef MINIJATI_FACTORYPIPNODE_HPP
#define MINIJATI_FACTORYPIPNODE_HPP

#include "PIPnode.hpp"
#include "FactoryPIPnodeCPU.hpp"
#include "FactoryPIPnodeRAM.hpp"
#include "FactoryPIPnodeSB.hpp"

namespace bpp {

    class nodeFactory { // PIPnode Factory class

    public:

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        // factory methd
        PIPnode *getPIPnode(enumDP3Dversion DPversion,
                            const progressivePIP *pPIP,
                            tshlib::VirtualNode *vnode,
                            bpp::Node *bnode) {

            switch (DPversion) {
                case CPU:
                    return new nodeCPU(pPIP, vnode, bnode);
                    break;
                case RAM:
                    return new nodeRAM(pPIP, vnode, bnode);
                    break;
                case SB:
                    return new nodeSB(pPIP, vnode, bnode);
                    break;
                default:
                    LOG(FATAL) << "\nERROR DP3D version not recognized";
                    return NULL;
            }

        }

    };

}

#endif //MINIJATI_FACTORYPIPNODE_HPP
