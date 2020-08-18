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


#include "tbb/task.h"
#include "tbb/task_group.h"
#include "tbb/tick_count.h"

#ifndef MINIJATI_COMPOSITEPIPNODE_TBB_HPP
#define MINIJATI_COMPOSITEPIPNODE_TBB_HPP

namespace bpp {

    class CompositePIPnodeTBB : public CompositePIPnode {

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


        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        CompositePIPnodeTBB(int numNodes, bool doTasks = false, bool doParallelFor = false)
        :CompositePIPnode(numNodes)
        {
            nodeTBB::setDoTasks(doTasks);
            nodeTBB::setDoParallelFor(doParallelFor);
            if (doTasks) {
                LOG(INFO) << "Using CompositePIPnodeTBB Tasks doParallelFor=" << doParallelFor;
            }else{
                LOG(INFO) << "Using CompositePIPnodeTBB NoTasks doParallelFor=" << doParallelFor;
            }
        }

        //void PIPnodeAlign(); // align all the input sequences under the PIP model
        void PIPnodeAlign() {

            // ALIGN UNDER THE PIP MODEL
            if (nodeTBB::isDoTasks()) {
                // TASKS start execution from root node
                // get the number of PIPnodes
                size_t num_nodes = pip_nodes_.size();

#if 0
                // REMARK DF:  bottomUP approach all leafs in //
                std::vector< PIPnode * > leafs;
                int cnt = num_nodes;
                for (int i = 0; i < num_nodes; i++) {
                    PIPnode *cur = pip_nodes_[i];
                    if (cur->childL_ == nullptr && cur->childR_ == nullptr) {
                        leafs.push_back( cur );
                        cnt--;
                    }
                }
                LOG(INFO) << "CompositePIPnodeTBB::PIPnodeAlign() BottomUp approach for "
                    << num_nodes << " nodes starting from "
                    << leafs.size() << " leafs";


                tbb::task_group childrenGroup;
                for (PIPnode *elem : leafs) {
                    VLOG(2) << "enqueuing Node-> " << elem->bnode_->getName() << std::endl;

                    childrenGroup.run(  [=] {
                        nodeTBB* node = static_cast<nodeTBB*>(elem);
                        node->DP3D_PIP_bottomUp();
                    } );
                }

                childrenGroup.wait();

#else
                LOG(INFO) << "CompositePIPnodeTBB::PIPnodeAlign() DivideEtImpera approach for " << num_nodes << " nodes";
                // Get Root node and start execution: each internalNode will first spawn children tasks and then run it's own work
                nodeTBB *root = static_cast<nodeTBB*>(pip_nodes_[num_nodes - 1]);
                root->DP3D_PIP();
#endif

            } else {
                // NonTasks execute common order
                CompositePIPnode::PIPnodeAlign();

            }
        }
    };
}

#endif //MINIJATI_COMPOSITEPIPNODE_TBB_HPP
