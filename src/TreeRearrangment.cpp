/***************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 * Copyright (C) 2015-2019 by Lorenzo Gatti & Massimo Maiolo
 ***************************************************************************
 * This file is part of Tree Search Heuristics Library (TSH-LIB)
 *
 * TSH-LIB is a C++ Library whose purpose is generate alternative tree
 * toplogies applying NNI, SPR, and F-VFNI rearrangements on an initial tree
 *
 * TSH-LIB is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * TSH-LIB is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public
 * License along with TSH-LIB. If not, see <http://www.gnu.org/licenses/>.
 ***************************************************************************/

/**
 * @file TreeRearrangment.cpp
 * @author Lorenzo Gatti
 * @date 18 10 2017
 * @version 3.0.2
 * @maintainer Lorenzo Gatti
 * @email lg@lorenzogatti.me
 * @status Development
 *
 * @brief
 * @details
 * @pre
 * @bug
 * @warning
 *
 * @see For more information visit: https://bitbucket.org/lorenzogatti89/tshlib
 */
#include <numeric>
#include <limits>
#include <iomanip>
#include <iterator>
#include <chrono>
#include <algorithm>
#include <glog/logging.h>
#include <map>
#include "TreeRearrangment.hpp"


namespace tshlib {
    TreeRearrangment::TreeRearrangment() = default;

    TreeRearrangment::~TreeRearrangment() {

        for (std::vector<Move *>::reverse_iterator i = trMoveSet.rbegin(); i < trMoveSet.rend(); i++) {
            Move *move = *i;
            delete move;
        }
        std::vector<Move *>().swap(trMoveSet);


    };

    void TreeRearrangment::getNodesInRadiusDown(VirtualNode *node_source, int radius_min, int radius_curr, int radius_max, bool includeSelf,
                                                MoveDirections direction, bool allowDuplicatedMoves) {

        VirtualNode *node = node_source;

        if (!includeSelf) {
            includeSelf = true;

        } else {

            if (radius_curr <= (radius_max - radius_min) && radius_curr >= 0) {

                addMove(radius_max - radius_curr, node, direction, allowDuplicatedMoves);
            }
        }

        if (radius_curr < 0) {
            return;
        }

        if (!node->isTerminalNode()) {
            radius_curr--;
            // Direction is required to know where to look at in the tree when the move is performed
            getNodesInRadiusDown(node->getNodeLeft(), radius_min, radius_curr, radius_max, includeSelf, direction, allowDuplicatedMoves);
            getNodesInRadiusDown(node->getNodeRight(), radius_min, radius_curr, radius_max, includeSelf, direction, allowDuplicatedMoves);
        }


    }

    void
    TreeRearrangment::getNodesInRadiusUp(VirtualNode *node_source, int radius_min, int radius_curr, int radius_max, NodePosition traverse_direction,
                                         bool allowDuplicatedMoves) {

        VirtualNode *vnode;

        NodePosition idx;
        MoveDirections moving_direction;


        vnode = node_source;

        //TODO: check binary tree condition!
        if (radius_curr <= (radius_max - radius_min) && radius_curr >= 0) {

            switch (traverse_direction) {
                case NodePosition::left:
                    moving_direction = MoveDirections::up_left;
                    break;
                case NodePosition::right:
                    moving_direction = MoveDirections::up_right;
                    break;
                case NodePosition::up:
                    moving_direction = MoveDirections::up;
                    break;
                default:
                    moving_direction = MoveDirections::undef;

            }

            addMove(radius_max - radius_curr, vnode, moving_direction, allowDuplicatedMoves);
        }

        if (radius_curr > 0) {
            if (traverse_direction == NodePosition::left) {
                radius_curr--;
                if (vnode->getNodeUp() != nullptr) {

                    idx = vnode->indexOf();
                    getNodesInRadiusUp(vnode->getNodeUp(), radius_min, radius_curr, radius_max, idx, allowDuplicatedMoves);
                }

                getNodesInRadiusDown(vnode->getNodeRight(), radius_min, radius_curr, radius_max, true, MoveDirections::up, allowDuplicatedMoves);

            } else if (traverse_direction == NodePosition::right) {
                radius_curr--;
                if (vnode->getNodeUp() != nullptr) {

                    idx = vnode->indexOf();
                    getNodesInRadiusUp(vnode->getNodeUp(), radius_min, radius_curr, radius_max, idx, allowDuplicatedMoves);
                }

                getNodesInRadiusDown(vnode->getNodeLeft(), radius_min, radius_curr, radius_max, true, MoveDirections::up, allowDuplicatedMoves);

            } else if (traverse_direction == NodePosition::up) {

                // If the traverse_direction is 2, we are moving across the pseudoroot, therefore no need to decrease the radious
                if (vnode->getNodeUp() != nullptr) {

                    // Get the nodes in the radius from the node we reached after crossing the pseudoroot
                    getNodesInRadiusDown(vnode, radius_min, radius_curr - 1, radius_max - 1, false, MoveDirections::up, allowDuplicatedMoves);

                }

            }
        }
    }

    void TreeRearrangment::defineMoves(bool includeSelf, bool allowDuplicatedMoves) {
        // Flag the nodes according to their position on the tree (left or right or above the source node -- p node).
        // For each node within the radius extremities, define a move and add it to TreeRearrangment.
        // Start from the children of the current starting node (if any)

        trCandidateMovesFound_ = 0;

        if (!utree_->getNode(trSourceNode_)->isTerminalNode()) {

            getNodesInRadiusDown(utree_->getNode(trSourceNode_)->getNodeLeft(),
                                 trSearchRadius_min,
                                 trSearchRadius_max - 1,
                                 trSearchRadius_max,
                                 includeSelf,
                                 MoveDirections::down_left,
                                 allowDuplicatedMoves);

            getNodesInRadiusDown(utree_->getNode(trSourceNode_)->getNodeRight(),
                                 trSearchRadius_min,
                                 trSearchRadius_max - 1,
                                 trSearchRadius_max,
                                 includeSelf,
                                 MoveDirections::down_right,
                                 allowDuplicatedMoves);
        }
        // If the node is a leaf, then go up
        if (nullptr != utree_->getNode(trSourceNode_)->getNodeUp()) {

            getNodesInRadiusUp(utree_->getNode(trSourceNode_)->getNodeUp(),
                               trSearchRadius_min,
                               trSearchRadius_max - 1,
                               trSearchRadius_max,
                               utree_->getNode(trSourceNode_)->indexOf(),
                               allowDuplicatedMoves);
        }


        VLOG(1) << "[TSH Cycle]  Found " << trCandidateMovesFound_ << " candidate moves for node " << utree_->getNode(trSourceNode_)->getNodeName();
    }

    void TreeRearrangment::addMove(int radius, VirtualNode *targetNode, MoveDirections moveDirection, bool allowDuplicatedMoves) {

        std::vector<TreeSearchHeuristics> strategies;

        if (trStrategy == TreeSearchHeuristics::mixed) {
            strategies.push_back(TreeSearchHeuristics::swap);
            strategies.push_back(TreeSearchHeuristics::phyml);
        } else {
            strategies.push_back(trStrategy);
        }

        for (auto &ts_strategy:strategies) {

            // Skip SPR-like moves if the radius is insufficient
            if (ts_strategy == TreeSearchHeuristics::phyml && radius < 4) continue;

            auto moveInstance = new Move;
            moveInstance->initMove();
            moveInstance->setSourceNode(trSourceNode_);
            moveInstance->setDirection(moveDirection);
            moveInstance->setRadius(radius);
            moveInstance->setTargetNode(targetNode->getVnode_id());
            moveInstance->setClass(ts_strategy, (utree_->getNode(trSourceNode_)->isPseudoRootNode() || targetNode->isPseudoRootNode()));


            bool storeMove = true;

            if (!allowDuplicatedMoves) {
                for (auto &query:trMoveSet) {

                    if (query->getTargetNode() == targetNode->getVnode_id()
                        && query->getSourceNode() == trSourceNode_
                        && query->getMoveStrategy() == ts_strategy) {

                        storeMove = false;
                    }

                    if (targetNode->getVnode_id() == query->getSourceNode()
                        && trSourceNode_ == query->getTargetNode()
                        && query->getMoveStrategy() == ts_strategy) {
                        storeMove = false;
                    }

                }
            }

            if (storeMove) {
                moveInstance->moveUID_ = (int) trMoveSet.size();
                trMoveSet.push_back(moveInstance);
                trCandidateMovesFound_++;
            } else {
                delete moveInstance;
            }
        }
    }

    void TreeRearrangment::printMoves() {

        VLOG(2) << "[set " << trUID_ << "] " << getStrategy() << " strategy" << std::endl;
        VLOG(2) << "[class]\t(P\t; Q)" << std::endl;
        for (auto &nmove: trMoveSet) {
            VLOG(2) << "[" << nmove->moveClassDescription_ << "]\t" << nmove->moveRadius_
                    << "\t(" << utree_->getNode(trSourceNode_)->vnode_name << "\t; "
                    << utree_->getNode(nmove->getTargetNode())->vnode_name << ")\t"
                    << static_cast<int>(nmove->moveDirection_) << std::endl;
        }

    }

    unsigned long TreeRearrangment::getNumberOfMoves() {

        return trMoveSet.size();
    }

    bool TreeRearrangment::applyMove(unsigned long moveID, Utree &_utree__topology) {

        bool outcomeExecutionMove = false;

        //VirtualNode *pnode = trMoveSet.at(moveID)->getSourceNode();
        //VirtualNode *qnode = trMoveSet.at(moveID)->getTargetNode();

        int pnode = trMoveSet.at(moveID)->getSourceNode();
        int qnode = trMoveSet.at(moveID)->getTargetNode();

        VLOG(2) << "[tshlib::applyMove] Debug: [MOVE#" << trMoveSet.at(moveID)->getUID() << "] direction: " << trMoveSet.at(moveID)->getDirection();
        VLOG(2) << "[tshlib::applyMove] Debug: [MOVE#"
                << trMoveSet.at(moveID)->getUID()
                << "] S: " << _utree__topology.getNode(pnode)->getNodeName()
                << " T: " << _utree__topology.getNode(qnode)->getNodeName();


        bool revertRotations;

        switch (trMoveSet.at(moveID)->getType()) {

            case MoveType::VFNNI:

            case MoveType::FNNI:

            case MoveType::NNI:

                revertRotations = false;
                // Swap pnode with qnode according to the direction found during the move configuration
                // If the swap is performed correctly then the function returns true otherwise false
                outcomeExecutionMove = _utree__topology.getNode(pnode)->swapNode(_utree__topology.getNode(qnode),
                                                                                 trMoveSet.at(moveID)->moveDirection_,
                                                                                 revertRotations);

                break;

            case MoveType::SPR:
                // Swap pnode with qnode according to the direction found during the move configuration
                // If the swap is performed correctly then the function returns true otherwise false
                outcomeExecutionMove = _applySPR(trMoveSet.at(moveID), _utree__topology);
                _utree__topology.printAllNodesNeighbors();
                break;

            case MoveType::TBR:
                outcomeExecutionMove = false;
                LOG(WARNING) << "[tshlib::TreeRearrangment::applyMove] Move [" << moveID << "] is TBR. Not implemented method to apply it.";

                break;


            case MoveType::undef:
                outcomeExecutionMove = false;
                LOG(FATAL) << "[tshlib::TreeRearrangment::applyMove] Something went wrong during the application of the move [" << moveID
                           << "]. It looks like its type is undefined!";

                break;
        }

        return outcomeExecutionMove;

    }

    bool TreeRearrangment::revertMove(unsigned long moveID, Utree &_utree__topology) {

        bool outcomeExecutionMove = false;

        //VirtualNode *pnode = trMoveSet.at(moveID)->getTargetNode();
        //VirtualNode *qnode = trMoveSet.at(moveID)->getSourceNode();

        int pnode = trMoveSet.at(moveID)->getTargetNode();
        int qnode = trMoveSet.at(moveID)->getSourceNode();

        bool revertRotations;

        switch (trMoveSet.at(moveID)->getType()) {

            case MoveType::VFNNI:

            case MoveType::FNNI:

            case MoveType::NNI:
                revertRotations = true;
                // Swap pnode with qnode according to the direction found during the move configuration
                // If the swap is performed correctly then the function returns true otherwise false
                outcomeExecutionMove = _utree__topology.getNode(pnode)->swapNode(_utree__topology.getNode(qnode),
                                                                                 MoveDirections::up,
                                                                                 revertRotations);

                break;

            case MoveType::SPR:
                // Swap pnode with qnode according to the direction found during the move configuration
                // If the swap is performed correctly then the function returns true otherwise false
                outcomeExecutionMove = _revertSPR(trMoveSet.at(moveID), _utree__topology);
                _utree__topology.printAllNodesNeighbors();
                break;

            case MoveType::TBR:
                outcomeExecutionMove = false;
                LOG(WARNING) << "[tshlib::TreeRearrangment::applyMove] Move [" << moveID << "] is TBR. Not implemented method to revert it.";

                break;


            case MoveType::undef:
                outcomeExecutionMove = false;
                LOG(FATAL) << "[tshlib::TreeRearrangment::applyMove] Something went wrong during the application of the move [" << moveID
                           << "]. It looks like its type is undefined!";

                break;
        }

        return outcomeExecutionMove;


    }

    Move *TreeRearrangment::getMove(unsigned long moveID) {

        return trMoveSet.at(moveID) ?: nullptr;
    }

    Move *TreeRearrangment::selectBestMove(double value = -std::numeric_limits<double>::infinity()) {

        Move *selectedMove = nullptr;
        for (auto &move:trMoveSet) {
            if (move->moveScore_ > value) {
                selectedMove = move;
                value = move->moveScore_;
            }
        }

        return selectedMove;

    }

    void TreeRearrangment::commitMove(int moveID, Utree &_utree__topology) {

        // apply move
        applyMove(moveID, _utree__topology);

        // reset node rotations
        for (auto &nodeID:_utree__topology.getNodeIdsMap()) {
            _utree__topology.getNode(nodeID.first)->vnode_rotated = NodeRotation::undef;
        }

        // reset node rotations
        //for (auto &node:utree_->listVNodes) {
        //    utree_->getNode(node)->vnode_rotated = NodeRotation::undef;
        //}

    }

    void TreeRearrangment::storeMove(Move *inMove) {

        inMove->moveUID_ = (int) trMoveSet.size();
        trMoveSet.push_back(inMove);

    }

    void TreeRearrangment::setTree(Utree *inTree) {
        utree_ = inTree;
    }

    void TreeRearrangment::displayRearrangmentStatus(int idMove, Utree & _utree__topology, bool printTree) {

        std::string start_col_line, end_col_line;

        // ------------------------------------
        // Some abbellishments for the console output
        if (getMove(idMove)->moveScore_ > 0) {
            start_col_line = "\033[1;34m";
            end_col_line = "\033[0m";
        } else {
            start_col_line = "";
            end_col_line = "";

        }

        // ------------------------------------
        // Move exection details
        if (printTree) {
            VLOG(2) << "[test  move]\t" << getMove(idMove)->moveClassDescription_ << "." << std::setfill('0') << std::setw(3) << idMove
                    << " | " << start_col_line << getMove(idMove)->moveScore_ << end_col_line << "\t"
                    << " | (" << _utree__topology.getNode(getMove(idMove)->getSourceNode())->vnode_name << "->"
                    << _utree__topology.getNode(getMove(idMove)->getTargetNode())->vnode_name << ")"
                    << "\t[" << getMove(idMove)->moveRadius_ << "] | " << getTree()->printTreeNewick(true);
        } else {
            VLOG(2) << "[test  move]\t" << getMove(idMove)->moveClassDescription_ << "." << std::setfill('0') << std::setw(3) << idMove
                    << " | " << start_col_line << getMove(idMove)->moveScore_ << end_col_line << "\t"
                    << " | (" << _utree__topology.getNode(getMove(idMove)->getSourceNode())->vnode_name
                    << "->" << _utree__topology.getNode(getMove(idMove)->getTargetNode())->vnode_name << ")"
                    << "\t[" << getMove(idMove)->moveRadius_ << "]";
        }

    }

    const std::vector<int> TreeRearrangment::updatePathBetweenNodes(unsigned long moveID, std::vector<int> inPath) {

        // Retrieve move object
        Move *move = getMove(moveID);

        // Declaring temporary objects
        std::vector<int> tmpVector_B, tmpVector_C, updatedNodesInPath, outNodePath;
        std::ptrdiff_t pos;

        // According to the different kind of move, apply the specific method to compute the list of nodes
        switch (move->getType()) {
            case MoveType::SPR:

                tmpVector_B = inPath;

                if (move->moveDirection_ == MoveDirections::up_left || move->moveDirection_ == MoveDirections::up_right) {
                    pos = std::distance(tmpVector_B.begin(), std::find(tmpVector_B.begin(), tmpVector_B.end(), move->getTargetNode()));
                } else {
                    pos = std::distance(tmpVector_B.begin(), std::find(tmpVector_B.begin(), tmpVector_B.end(), move->getSourceNode()));
                }

                if (pos <= tmpVector_B.size()) {

                    switch (move->getMoveDirection()) {

                        case MoveDirections::down_right:
                        case MoveDirections::down_left:

                            // Copy all the VirtualNode-pointers from the point where Source is found to the end
                            for (std::ptrdiff_t i = pos; i < tmpVector_B.size(); i++) {
                                tmpVector_C.push_back(tmpVector_B.at(i));
                            }

                            // Revert the order of the elements in the vector
                            std::reverse(std::begin(tmpVector_C), std::end(tmpVector_C));

                            // Add
                            tmpVector_C.push_back(tmpVector_B.at(pos - 1));

                            // Add the beginning of the vector B to vector C
                            for (std::ptrdiff_t i = 1; i < pos - 1; i++) {
                                tmpVector_C.push_back(tmpVector_B.at(i));
                            }

                            break;

                        case MoveDirections::up_right:
                        case MoveDirections::up_left:

                            // Add the path between steparent and insertion point
                            for (std::ptrdiff_t i = 2; i < pos; i++) {
                                tmpVector_C.push_back(tmpVector_B.at(i));
                            }

                            // Add target node
                            tmpVector_C.push_back(tmpVector_B.at(pos));

                            // Add source and sourceparent nodes
                            tmpVector_C.push_back(tmpVector_B.at(0));
                            tmpVector_C.push_back(tmpVector_B.at(1));

                            // Add path from insertion point to root
                            for (std::ptrdiff_t i = pos + 1; i < tmpVector_B.size(); i++) {
                                tmpVector_C.push_back(tmpVector_B.at(i));
                            }

                            break;

                        case MoveDirections::up:

//                            if(pos+1 == tmpVector_B.size()){
//
//                                // Add
//                                tmpVector_C.push_back(tmpVector_B.at(pos - 1));
//
//                                // Add the beginning of the vector B to vector C
//                                for (std::ptrdiff_t i = 1; i < pos - 1; i++) {
//                                    tmpVector_C.push_back(tmpVector_B.at(i));
//                                }
//
//                            }else {

                                //tmpVector_C.push_back(tmpVector_B.at(pos + 1));
                                tmpVector_C.push_back(tmpVector_B.at(pos + 1));

                                // Add the beginning of the vector B to vector C
                                for (std::ptrdiff_t i = 1; i < pos; i++) {
                                    tmpVector_C.push_back(tmpVector_B.at(i));
                                }

                                // Add the beginning of the vector B to vector C
                                for (std::ptrdiff_t i = pos + 2; i < tmpVector_B.size(); i++) {
                                    tmpVector_C.push_back(tmpVector_B.at(i));
                                }
//                            }

                            break;

                        case MoveDirections::undef:
                            DLOG(WARNING) << "[tshlib::TreeRearrangment::updatePathBetweenNodes] The current type of move is not yet supported!";
                            break;
                    }


                    // Return the newly ordered vector
                    updatedNodesInPath = tmpVector_C;
                }


                outNodePath = updatedNodesInPath;
                break;

            case MoveType::NNI:
            case MoveType::FNNI:
            case MoveType::VFNNI:
                if (move->moveDirection_ != MoveDirections::up) {

                    tmpVector_B = inPath;

                    // Remove the first element of the array since it is not a likelihood component (source node)
                    tmpVector_B.erase(tmpVector_B.begin());

                    // Find the position of the target node (it is not necessarily the end of the vector)


                    if (move->moveDirection_ == MoveDirections::up_left || move->moveDirection_ == MoveDirections::up_right) {
                        pos = std::distance(tmpVector_B.begin(), std::find(tmpVector_B.begin(), tmpVector_B.end(), move->getTargetNode()));
                    } else {
                        pos = std::distance(tmpVector_B.begin(), std::find(tmpVector_B.begin(), tmpVector_B.end(), move->getSourceNode()));
                    }

                    if (pos <= tmpVector_B.size()) {

                        // Copy the reference to the pointers of the node starting from the target node to the end of the vector (root)
                        for (std::ptrdiff_t i = pos; i < tmpVector_B.size(); i++) {
                            tmpVector_C.push_back(tmpVector_B.at(i));
                        }

                        // Revert the order of the elements in the vector
                        std::reverse(std::begin(tmpVector_C), std::end(tmpVector_C));

                        // Add the beginning of the vector B to vector C
                        for (std::ptrdiff_t i = 0; i < pos; i++) {
                            tmpVector_C.push_back(tmpVector_B.at(i));
                        }

                        // Return the newly ordered vector
                        updatedNodesInPath = tmpVector_C;
                    }

                    outNodePath = updatedNodesInPath;

                } else {

                    tmpVector_B = inPath;

                    tmpVector_B.erase(std::find(tmpVector_B.begin(), tmpVector_B.end(), move->getTargetNode()));
                    tmpVector_B.erase(std::find(tmpVector_B.begin(), tmpVector_B.end(), move->getSourceNode()));

                    outNodePath = tmpVector_B;
                }

                break;

            case MoveType::TBR:
            case MoveType::undef:
                DLOG(WARNING) << "[tshlib::TreeRearrangment::updatePathBetweenNodes] The current type of move is not yet supported!";

                break;
        }

        return outNodePath;
    }

    void TreeRearrangment::initialize() {

        // Define the radius for pruning and regrafting the input tree.
        switch (trTreeCoverage) {

            case tshlib::TreeRearrangmentOperations::classic_NNI:
                trSearchRadius_min = 3;
                trSearchRadius_max = 3;
                break;

            case tshlib::TreeRearrangmentOperations::classic_SPR:
                trSearchRadius_min = 4;
                trSearchRadius_max = utree_->getMaxNodeDistance() / 2;
                break;

            case tshlib::TreeRearrangmentOperations::classic_TBR:
                trSearchRadius_min = 5;
                trSearchRadius_max = utree_->getMaxNodeDistance() / 2;
                break;

            case tshlib::TreeRearrangmentOperations::classic_Mixed:
                trSearchRadius_min = 3;  // Minimum radius for an NNI move is 3 nodes
                trSearchRadius_max = utree_->getMaxNodeDistance(); // Full tree traversing from any nodeInterface of the tree
                break;

        }

        VLOG(1) << "[tsh] Defined coverage radius as [" << trSearchRadius_min << ";" << trSearchRadius_max << "] on a max of ["
                << utree_->getMaxNodeDistance() << "]";

        trInitialized_ = true;
    }

    bool TreeRearrangment::_applySPR(Move *move, Utree &_utree__topology) {

        bool execstatus = false;

        VirtualNode *sourceNode = _utree__topology.getNode(move->getSourceNode());
        VirtualNode *targetNode = _utree__topology.getNode(move->getTargetNode());

        if (sourceNode == targetNode) {

            LOG(FATAL) << "[tshlib::_applySPR] The source and the target nodes must be different!";

        } else if (targetNode == nullptr) {

            LOG(FATAL) << "[tshlib::_applySPR] The target node is empty";


        } else {

            // References
            VirtualNode *moveStepChild = nullptr;
            VirtualNode *moveStepParent = nullptr;
            VirtualNode *parentTarget = targetNode->getNodeUp();
            VirtualNode *parentSource = sourceNode->getNodeUp();

            std::string parentS_RN = parentSource->isPseudoRootNode() ? "yes" : "no";
            std::string parentT_RN = parentTarget->isPseudoRootNode() ? "yes" : "no";

            // Assign flags to nodes

            switch (move->getMoveDirection()) {
                case MoveDirections::up:
                case MoveDirections::up_left:
                case MoveDirections::up_right:

                    // Assign flags to nodes
                    moveStepChild = sourceNode->getSiblingNode();
                    moveStepParent = parentSource->getNodeUp();

                    // Disconnections
                    parentSource->disconnectNode();
                    moveStepChild->disconnectNode();
                    targetNode->disconnectNode();

                    // Connect step-child and step-parent
                    if (!moveStepParent->getNodeUp() && !moveStepParent->getNodeLeft() && !moveStepParent->getNodeRight()) {
                        moveStepChild->connectNode(moveStepParent);
                    } else {
                        moveStepParent->connectNode(moveStepChild);
                    }

                    break;

                case MoveDirections::down_left:
                case MoveDirections::down_right:

                    // Rotations //TODO: Check for rotation direction (it could be the way around)
                    if (move->getMoveDirection() == MoveDirections::down_left) {
                        sourceNode->rotateCounterClockwise();
                        //sourceNode->rotateClockwise();
                    } else {
                        //sourceNode->rotateCounterClockwise();
                        sourceNode->rotateClockwise();
                    }
                    // Resolve parent-source node
                    parentSource = sourceNode->getNodeUp();

                    // Assign flags to nodes
                    moveStepChild = parentSource->getNodeRight();
                    moveStepParent = parentSource->getNodeLeft();

                    // Disconnections
                    moveStepChild->disconnectNode();
                    moveStepParent->disconnectNode();
                    targetNode->disconnectNode();

                    // Re-solve
                    VirtualNodeUtils::rotateNodeClockwise(parentSource);

                    // Connect step-child and step-parent
                    moveStepChild->_bidirectionalUpwardConnection(moveStepParent);

                    break;

                default:
                    LOG(ERROR) << "[thslib::_applySPR] Move #" << move->getUID() << " has an undefined direction.";
                    break;
            }

            // Re-Connections
            parentTarget->connectNode(parentSource);
            parentSource->connectNode(targetNode);


            VLOG(2) << "[tshlib::_applySPR] Debug: [MOVE#" << move->getUID() << "] pS: " << parentSource->getNodeName() << "(" << parentS_RN << ")" <<
                    " pT: " << parentTarget->getNodeName() << "(" << parentT_RN << ")";
            VLOG(2) << "[tshlib::_applySPR] Debug: [MOVE#" << move->getUID() << "] stepParent: " << moveStepParent->getNodeName() << " stepChild: "
                    << moveStepChild->getNodeName();


            // Set StepParent and StepChild in move description
            move->setStepChildNode(moveStepChild->getVnode_id());
            move->setStepParentNode(moveStepParent->getVnode_id());

            execstatus = true;

        }
        return execstatus;
    }

    bool TreeRearrangment::_revertSPR(Move *move, Utree &_utree__topology) {

        bool execstatus = false;

        VirtualNode *sourceNode = _utree__topology.getNode(move->getSourceNode());
        VirtualNode *parentSourceNode = _utree__topology.getNode(move->getSourceNode())->getNodeUp();

        VirtualNode *targetNode = _utree__topology.getNode(move->getTargetNode());
        VirtualNode *parentTargetNode = _utree__topology.getNode(move->getSourceNode())->getNodeUp()->getNodeUp();

        VirtualNode *moveStepParent = _utree__topology.getNode(move->getStepParentNode());
        VirtualNode *moveStepChild = _utree__topology.getNode(move->getStepChildNode());


        if (sourceNode == targetNode) {

            LOG(FATAL) << "[tshlib::sprNode_revert] The source and the target nodes must be different!";

        } else if (targetNode == nullptr) {

            LOG(FATAL) << "[tshlib::sprNode_revert] The target node is empty";


        } else {

            // 1. Disconnections
            moveStepChild->disconnectNode();
            parentSourceNode->disconnectNode();
            targetNode->disconnectNode();

            // 2. Revert rotations on the parentSourceNode (resolve)
            bool rotatedcase = false;
            switch (sourceNode->getNodeRotation()) {
                case NodeRotation::counterclockwise:
                    sourceNode->rotateClockwise(true);
                    VirtualNodeUtils::rotateNodeCounterClockwise(parentSourceNode);
                    rotatedcase = true;
                    break;
                case NodeRotation::clockwise:
                    sourceNode->rotateCounterClockwise(true);
                    VirtualNodeUtils::rotateNodeCounterClockwise(parentSourceNode);
                    rotatedcase = true;
                    break;
                case NodeRotation::undef:
                    break;
            }

            // 4. Reconnections
            if (rotatedcase) {
                parentSourceNode->connectNode(moveStepParent);
            } else {
                if (!moveStepParent->getNodeLeft() && !moveStepParent->getNodeRight()) {
                    moveStepParent->_setNodeUp(parentSourceNode);
                    parentSourceNode->_setNodeUp(moveStepParent);
                } else {
                    moveStepParent->connectNode(parentSourceNode);
                }
            }

            parentSourceNode->connectNode(moveStepChild);

            parentTargetNode->connectNode(targetNode);

            execstatus = true;
        }
        return execstatus;
    }

}



