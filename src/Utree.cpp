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
 * @file Utree.cpp
 * @author Lorenzo Gatti
 * @date 26 10 2017
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

#include <string>
#include <random>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <glog/logging.h>
#include <map>
#include "Utree.hpp"


#include "Utree.hpp"

using namespace tshlib;

VirtualNode *UtreeUtils::getPseudoRoot(VirtualNode *vn) {

    VirtualNode *vtemp = vn;
    while (!vtemp->isRootNode()) {
        vtemp = vtemp->getNodeUp();
    }

    return vtemp;
}

namespace tshlib {

    VirtualNode::VirtualNode() {
        // Initialise a new VirtualNode completely disconnected from the tree
        vnode_right = nullptr;
        vnode_left = nullptr;
        vnode_up = nullptr;
        vnode_branchlength = 0;
        vnode_leaf = false;
        vnode_rotated = NodeRotation::undef;


    };


    void VirtualNode::connectNode(VirtualNode *inVNode) {

        // Connect this node to the parent node
        inVNode->_setNodeUp(this);
        //inNode->_oneway_connectNode(this);
        _oneway_connectNode(inVNode);
    }


    bool VirtualNode::isTerminalNode() {

        return vnode_leaf;
    }


    bool VirtualNode::isRootNode() {

        return getNodeUp() == nullptr;
    }


    void VirtualNode::setNodeName(const std::string s) {

        vnode_name = s;
    }


    void VirtualNode::_setNodeRight(VirtualNode *inVNode) {

        vnode_right = inVNode;

    }


    void VirtualNode::_setNodeLeft(VirtualNode *inVNode) {

        vnode_left = inVNode;

    }


    void VirtualNode::_setNodeUp(VirtualNode *inVNode) {

        vnode_up = inVNode;

    }


    const std::string VirtualNode::getNodeName() {

        return vnode_name;
    }


    VirtualNode *VirtualNode::getNodeUp() {

        return vnode_up ?: nullptr;
    }


    void VirtualNode::rotateClockwise() {

        _recursive_cw_rotation(this, false);
    }


    void VirtualNode::rotateClockwise(bool revertRotations) {

        _recursive_cw_rotation(this, revertRotations);
    }


    void VirtualNode::rotateCounterClockwise() {

        _recursive_ccw_rotation(this, false);
    }


    void VirtualNode::rotateCounterClockwise(bool revertRotations) {

        _recursive_ccw_rotation(this, revertRotations);
    }


    VirtualNode *VirtualNode::getNodeLeft() {

        return vnode_left ?: nullptr;
    }


    VirtualNode *VirtualNode::getNodeRight() {

        return vnode_right ?: nullptr;
    }


    void VirtualNode::_traverseVirtualNodeTree() {

        std::cout << "[Name] " << vnode_name << std::endl;

        if (vnode_up) {
            std::cout << "[Up name] " << getNodeUp()->vnode_name << std::endl;
        } else {
            std::cout << "[Up name] " << " ' ' " << std::endl;
        }

        std::cout << "[Node branch length] " << vnode_branchlength << std::endl;

        if (isTerminalNode()) {
            std::cout << std::endl;
        } else {

            std::cout << "[Left name] " << getNodeLeft()->vnode_name << std::endl;
            std::cout << "[Right name] " << getNodeRight()->vnode_name << std::endl;
            std::cout << std::endl;

            getNodeLeft()->_traverseVirtualNodeTree();
            getNodeRight()->_traverseVirtualNodeTree();
        }


    }


    double VirtualNode::computeTotalTreeLength() {

        if (isTerminalNode()) {
            return vnode_branchlength;
        } else {
            return vnode_branchlength +
                   getNodeLeft()->computeTotalTreeLength() +
                   getNodeRight()->computeTotalTreeLength();
        }
    }


    void VirtualNode::clearChildren() {
        vnode_left = nullptr;
        vnode_right = nullptr;
    }



    void Utree::addMember(VirtualNode *inVNode, bool isStartNode) {

        listVNodes.push_back(inVNode);
        if (isStartNode) {
            startVNodes.push_back(inVNode);
        }

        // Add pointer reference to the utree map
        utreeNodeIdsMap_[inVNode->getVnode_id()] = inVNode;

    }

    std::vector<int> Utree::getChildrenId(int nodeID) const {
        std::vector<int> sonsId(2,-2);

        VirtualNode *node = &getNode(nodeID);

        sonsId[0] = node->getNodeLeft()->getVnode_id();
        sonsId[1] = node->getNodeRight()->getVnode_id();

        return sonsId;
    }

    std::vector<VirtualNode *> Utree::findPseudoRoot(VirtualNode *inVNode, bool fixPseudoRootOnNextSubtree) {

        // This method returns a vector which contains all the nodes in the path
        // from the start node to the pseudoroot
        std::vector<VirtualNode *> path2root;

        VirtualNode *CurrentNode = inVNode;

        // Add the first element of the path to the root
        //path2root.push_back(inVNode);

        // Add all the other nodes traversing the tree in post-order
        while (CurrentNode != nullptr) {

            path2root.push_back(CurrentNode);

            if (CurrentNode == CurrentNode->getNodeUp()->getNodeUp()) { break; }


            if (CurrentNode->getNodeUp() != nullptr) {
                CurrentNode = CurrentNode->getNodeUp();
            } else {
                CurrentNode = nullptr;
            }

        }

        if (fixPseudoRootOnNextSubtree) {

            if (CurrentNode->getNodeUp() != nullptr) {
                path2root.push_back(CurrentNode->getNodeUp());
            }

        }

        return path2root;
    }

    void Utree::myAddRoot() {

        this->rootnode->_setNodeUp(nullptr);
        this->rootnode->_setNodeLeft(this->startVNodes.at(0));
        this->startVNodes.at(0)->_setNodeUp(this->rootnode);
        this->rootnode->_setNodeRight(this->startVNodes.at(1));
        this->startVNodes.at(1)->_setNodeUp(this->rootnode);

    }

    void Utree::myRemoveRoot() {

        this->startVNodes.at(0)->_setNodeUp(this->startVNodes.at(1));
        this->startVNodes.at(1)->_setNodeUp(this->startVNodes.at(0));

        this->rootnode->_setNodeLeft(nullptr);
        this->rootnode->_setNodeRight(nullptr);

    }

    std::string Utree::printTreeNewick(bool showInternalNodeNames, bool updateStartNodes) {

        try {
            std::string s, terminator;

            if(updateStartNodes) _updateStartNodes();

            s += "(";
            for (unsigned long i = 0; i < 2; i++) {
                if (i == 2 - 1) {
                    terminator = ");";
                } else {
                    terminator = ",";
                }
                s += _recursiveFormatNewick(startVNodes.at(i), showInternalNodeNames) + terminator;


            }
            return s;
        } catch (const std::exception &e) {

            LOG(FATAL) << "[Utree::printTreeNewick]" << e.what();

        }


    }


    std::string Utree::_recursiveFormatNewick(VirtualNode *vnode, bool showInternalNodeNames) {

        try {
            std::stringstream newick;

            if (vnode->isTerminalNode()) {

                newick << vnode->vnode_name << ":" << vnode->vnode_branchlength;

            } else {

                newick << "(";
                newick << _recursiveFormatNewick(vnode->getNodeLeft(), showInternalNodeNames);
                newick << ",";
                newick << _recursiveFormatNewick(vnode->getNodeRight(), showInternalNodeNames);
                newick << ")";

                if (showInternalNodeNames) {
                    newick << vnode->vnode_name;
                }

                newick << ":" << vnode->vnode_branchlength;
            }

            return newick.str();
        } catch (const std::exception &e) {

            LOG(FATAL) << "[Utree::printTreeNewick]" << e.what();

        }

    }


    Utree::~Utree() {

        for (std::vector<VirtualNode *>::reverse_iterator i = listVNodes.rbegin(); i < listVNodes.rend(); i++) {
            VirtualNode *vnode = *i;
            delete vnode;
        }
        delete rootnode;
        std::vector<VirtualNode *>().swap(listVNodes);

    }


    double Utree::computeTotalTreeLength() {
        double t_length;

        t_length = 0.0;
        for (unsigned long i = 0; i < listVNodes.size(); i++) {
            t_length += listVNodes.at(i)->vnode_branchlength;
        }

        return t_length;
    }


    void Utree::_printUtree() {

        VirtualNode *vn;

        for (unsigned long i = 0; i < listVNodes.size(); i++) {
            vn = listVNodes.at(i);

            std::cout << "[Node name] " << vn->vnode_name << std::endl;

            std::cout << "[Up child node name] " << vn->getNodeUp()->vnode_name << std::endl;

            if (!vn->isTerminalNode()) {
                std::cout << "[L. child node name] " << vn->getNodeLeft()->vnode_name << std::endl;
                std::cout << "[R. child node name] " << vn->getNodeRight()->vnode_name << std::endl;
            }
            std::cout << "[Node branch length] " << vn->vnode_branchlength << std::endl;

            std::cout << std::endl;
        }

    }

    Utree::Utree() {

        initialized_treeDepth = false;

        // Added virtual root to utree
        auto root = new VirtualNode;
        root->vnode_name = "root";
        root->vnode_id = -1;
        rootnode = root;

        // Add pointer reference to the utree map
        utreeNodeIdsMap_[rootnode->getVnode_id()] = rootnode;
    };


    void Utree::_updateStartNodes() {
        try {
            // Find the pseudoroot node on the first side of the tree
            bool fixPseudoRootOnNextSubtree = true;
            std::vector<VirtualNode *> sideA = findPseudoRoot(listVNodes.at(0), fixPseudoRootOnNextSubtree);
            VirtualNode *sideA_node = sideA.back();

            // Find the pseudoroot node on the opposite side of the tree
            fixPseudoRootOnNextSubtree = false;
            std::vector<VirtualNode *> sideB = findPseudoRoot(listVNodes.at(0), fixPseudoRootOnNextSubtree);
            VirtualNode *sideB_node = sideB.back();

            // Reset the original utree attribute
            startVNodes.clear();

            // Update the utree attribute
            startVNodes.push_back(sideA_node);
            startVNodes.push_back(sideB_node);

        } catch (const std::exception &e) {

            LOG(FATAL) << "[Utree::_updateStartNodes]" << e.what();

        }

    }


    void Utree::_testReachingPseudoRoot() {

        std::string strpath;
        std::vector<VirtualNode *> vnodes2root;

        for (auto &i : listVNodes) {

            vnodes2root = findPseudoRoot(i, true);

            for (auto &tnode: vnodes2root) {
                strpath += tnode->vnode_name;

                if (tnode->vnode_rotated != NodeRotation::undef) {

                    strpath += "\033[1;34m(" + std::to_string(static_cast<int>(tnode->vnode_rotated)) + ")\033[0m";

                } else {

                    strpath += "(" + std::to_string(static_cast<int>(tnode->vnode_rotated)) + ")";
                }

                strpath += ">";
            }
            strpath.pop_back();

            VLOG(1) << strpath << std::endl;
            vnodes2root.clear();
            strpath.clear();

        }

    }


    void Utree::saveTreeOnFile(std::string outfilepath) {

        std::ofstream outfile;

        outfile.open(outfilepath, std::ios_base::app);
        outfile << printTreeNewick(true) << std::endl;

    }


    int Utree::getMaxNodeDistance() {

        long int max_distance = 0;

        for (auto &vnode:listVNodes) {

            if (vnode->isTerminalNode()) {

                if (max_distance < findPseudoRoot(vnode).size()) {
                    max_distance = findPseudoRoot(vnode).size();
                };

            }
        }
        return max_distance * 2;
    }

    int Utree::mygetMaxNodeDistance() {

        long int max_distance = 0;

        int rootId1 = startVNodes.at(0)->vnode_id;
        int rootId2 = startVNodes.at(1)->vnode_id;

        for (auto &vnode:listVNodes) {

            if (vnode->isTerminalNode()) {

                VirtualNode *node = vnode;
                int currdist = 1;
                while(node->vnode_id != rootId1 && node->vnode_id != rootId2){
                    node = node->getNodeUp();
                    currdist++;
                }

                if (max_distance < currdist) {
                    max_distance = currdist;
                };

            }
        }

        return max_distance * 2;
    }

    std::string VirtualNode::printNeighbours() {

        std::stringstream description;

        if (!isTerminalNode()) {

            if (getNodeUp()) {
                description << vnode_name << " (^" << getNodeUp()->vnode_name << ";";
            } else {
                description << vnode_name << " (^NULL;";
            }

            description << "<" << getNodeLeft()->vnode_name << ";";
            description << getNodeRight()->vnode_name << ">)";

        } else {
            description << vnode_name << " (^" << getNodeUp()->vnode_name << "; ";
            description << "<-;->)";
        }

        return description.str();
    }


    NodePosition VirtualNode::indexOf() {
        NodePosition node_position;
        VirtualNode *parent = getNodeUp();
        if (parent->getNodeLeft() == this) {

            node_position = NodePosition::left;

        } else if (parent->getNodeRight() == this) {

            node_position = NodePosition::right;

        } else if (parent->getNodeUp() == this) {

            node_position = NodePosition::up;

        } else {
            //perror("The pointer index of the current node is unknown");
            node_position = NodePosition::undef;
        }
        return node_position;
    }

//    void VirtualNode::moveTheRoot(VirtualNode *pnode, VirtualNode *qnode, std::vector<VirtualNode *> &startVNodes) {
//
//        VirtualNode *new_root  = nullptr;
//        VirtualNode *tmp = nullptr;
//
//        if(pnode->isPseudoRootNode()){
//            tmp = pnode;
//            pnode = qnode;
//            qnode = tmp;
//        }
//
//        tmp = qnode->getNodeUp();
//
//        new_root = qnode->getNodeLeft();
//        new_root->_setNodeUp(qnode);
//        qnode->_setNodeLeft(tmp);
//        qnode->_setNodeUp(new_root);
//
//        tmp->_setNodeUp(qnode);
//
//        startVNodes.at(0)=new_root;
//        startVNodes.at(1)=new_root->getNodeUp();
//
//    }
//
//    void VirtualNode::NNImoveUP(VirtualNode *pnode, VirtualNode *qnode, std::vector<VirtualNode *> &startVNodes) {
//        VirtualNode * n1 = nullptr;
//        VirtualNode * n2 = nullptr;
//
//        n1 = pnode->getNodeUp();
//        n2 = qnode->getNodeUp();
//
//        if(pnode->vnode_id==n1->getNodeLeft()->vnode_id){
//            n1->_setNodeLeft(nullptr);
//            pnode->_setNodeUp(nullptr);
//        }else{
//            n1->_setNodeRight(nullptr);
//            pnode->_setNodeUp(nullptr);
//        }
//
//        if(qnode->vnode_id==n2->getNodeLeft()->vnode_id){
//            n2->_setNodeLeft(nullptr);
//            qnode->_setNodeUp(nullptr);
//        }else{
//            n2->_setNodeRight(nullptr);
//            qnode->_setNodeUp(nullptr);
//        }
//
//        if(!n2->getNodeLeft()){
//            n2->_setNodeLeft(pnode);
//            pnode->_setNodeUp(n2);
//        }else{
//            n2->_setNodeRight(pnode);
//            pnode->_setNodeUp(n2);
//        }
//
//        if(!n1->getNodeLeft()){
//            n1->_setNodeLeft(qnode);
//            qnode->_setNodeUp(n1);
//        }else{
//            n1->_setNodeRight(qnode);
//            qnode->_setNodeUp(n1);
//        }
//
//    }
//
//    void VirtualNode::NNImoveUP_LR(VirtualNode *pnode, VirtualNode *qnode, std::vector<VirtualNode *> &startVNodes) {
//        VirtualNode * n1 = nullptr;
//        VirtualNode * n2 = nullptr;
//        VirtualNode * tmp1 = nullptr;
//        VirtualNode * tmp2 = nullptr;
//        int id;
//        double bl;
//        std::string name;
//
//        n1 = pnode->getNodeUp();
//        n2 = n1->getNodeUp();
//
//        tmp1 = pnode->getSiblingNode();
//        tmp2 = n1->getSiblingNode();
//
//        if(tmp1->vnode_id==n1->getNodeLeft()->vnode_id){
//            tmp1->getNodeUp()->_setNodeLeft(nullptr);
//            tmp1->_setNodeUp(nullptr);
//        }else{
//            tmp1->getNodeUp()->_setNodeRight(nullptr);
//            tmp1->_setNodeUp(nullptr);
//        }
//
//        if(tmp2->vnode_id==n2->getNodeLeft()->vnode_id){
//            tmp2->getNodeUp()->_setNodeLeft(nullptr);
//            tmp2->_setNodeUp(nullptr);
//        }else{
//            tmp2->getNodeUp()->_setNodeRight(nullptr);
//            tmp2->_setNodeUp(nullptr);
//        }
//
//        if(!n2->getNodeLeft()){
//            n2->_setNodeLeft(tmp1);
//            tmp1->_setNodeUp(n2);
//        }else{
//            n2->_setNodeRight(tmp1);
//            tmp1->_setNodeUp(n2);
//        }
//
//        if(!n1->getNodeLeft()){
//            n1->_setNodeLeft(tmp2);
//            tmp2->_setNodeUp(n1);
//        }else{
//            n1->_setNodeRight(tmp2);
//            tmp2->_setNodeUp(n1);
//        }
//
//        id = n1->vnode_id;
//        bl = n1->vnode_branchlength;
//        name = n1->vnode_name;
//        n1->vnode_id = n2->vnode_id;
//        n1->vnode_branchlength = n2->vnode_branchlength;
//        n1->vnode_name = n2->vnode_name;
//        n2->vnode_id = id;
//        n2->vnode_branchlength = bl;
//        n2->vnode_name = name;
//
//    }
//
//    bool VirtualNode::swapNode(VirtualNode *targetNode, MoveDirections move_direction, bool revertRotations,std::vector<VirtualNode *> &startVNodes) {
//
//        bool execstatus = false;
//
//        if (this == targetNode) {
//
//            LOG(FATAL) << "[tshlib::swapNode] The source and the target nodes must be different!";
//
//
//        } else if (targetNode == nullptr) {
//
//            LOG(FATAL) << "[tshlib::swapNode] The target node is empty";
//
//
//        } else {
//
//            VirtualNode *pnode = this;
//            VirtualNode *qnode = targetNode;
//
//            switch (move_direction) {
//
//                case MoveDirections::up:
//
//                    if(pnode->isPseudoRootNode() || qnode->isPseudoRootNode()){
//
//                        moveTheRoot(pnode,qnode,startVNodes);
//
//                        NNImoveUP_LR(pnode,qnode,startVNodes);
//
//                    }else{
//
//                        NNImoveUP(pnode,qnode,startVNodes);
//
//                    }
//
//                    break;
//
//                case MoveDirections::down_left:
//                case MoveDirections::down_right:
//                case MoveDirections::up_left:
//                case MoveDirections::up_right:
//
//                    NNImoveUP_LR(pnode,qnode,startVNodes);
//
//                    break;
//
//                default:
//                    std::cout << "I cannot move the node" << std::endl;
//                    exit(EXIT_FAILURE);
//            }
//
//            execstatus = true;
//        }
//
//        return execstatus;
//    }


    void VirtualNode::disconnectNode() {

        // 1. Get the location where this node is connected on the parent node and disconnect
        switch (indexOf()) {
            case NodePosition::left:
                // The node is connected on the left side
                getNodeUp()->_setNodeLeft(nullptr);
                break;
            case NodePosition::right:
                // The node is connected on the right side
                getNodeUp()->_setNodeRight(nullptr);
                break;
            case NodePosition::up:
                // The node is connect on the up side (pseudoroot)
                getNodeUp()->_setNodeUp(nullptr);
                break;
            default:
                break;
        }

        // 2. Disconnect the parent node
        _setNodeUp(nullptr);
    }


    void VirtualNode::_oneway_connectNode(VirtualNode *inVNode) {

        // Check which direction is still available (either left or right)
        if (!getNodeLeft()) {
            _setNodeLeft(inVNode);
            return;
        } else if (!getNodeRight()) {
            _setNodeRight(inVNode);
            return;
        } else if (!getNodeUp()) {
            _setNodeUp(inVNode);
            return;
        } else {
            LOG(ERROR) << "[VirtualNode::_oneway_connectNode] No direction available in the VirtualNode to add a new child";
        }

    }

    void VirtualNode::_bidirectionalUpwardConnection(VirtualNode *inNode) {

        _setNodeUp(inNode);
        inNode->_setNodeUp(this);

    }

    bool VirtualNode::isParent(VirtualNode *inVNode) {

        bool status = false;

        VirtualNode *CurrentNode = inVNode;

        do {
            CurrentNode = CurrentNode->getNodeUp();
            if (CurrentNode == this) {
                status = true;
                break;
            };
        } while (CurrentNode != CurrentNode->getNodeUp()->getNodeUp());

        return status;
    }

    void Utree::printAllNodesNeighbors() {
        std::string token = "";
        for (auto &node:listVNodes) {

            if (node->isPseudoRootNode()) token = "***";
            VLOG(3) << node->printNeighbours() << token;
            token = "";
        }

    }

    void Utree::addVirtualRootNode() {

        if (rootnode->getNodeRight() == nullptr && rootnode->getNodeLeft() == nullptr) {

            // Update starting nodes in the startVNodes array
            _updateStartNodes();

            // Clear children nodes in the virtual root
            rootnode->clearChildren();

            // Disconnect the bidirectional connection in the pseudoroot nodes
            startVNodes.at(0)->disconnectNode();

            // Reconnect nodes to virtual root node
            rootnode->connectNode(startVNodes.at(0));
            rootnode->connectNode(startVNodes.at(1));
        }

    }

    void Utree::removeVirtualRootNode() {

        if (rootnode->getNodeRight() != nullptr && rootnode->getNodeLeft() != nullptr) {
            // Get reference children nodes from virtual root
            VirtualNode *left = rootnode->getNodeLeft();
            VirtualNode *right = rootnode->getNodeRight();

            // Disconnect nodes from virtual root node
            left->disconnectNode();
            right->disconnectNode();

            // Establish bidirectional connection
            left->_setNodeUp(right);
            right->_setNodeUp(left);

            // Clear residual connection with virtual root node
            rootnode->clearChildren();

        }
    }

    std::vector<int> Utree::computePathBetweenNodes2(VirtualNode *vnode_1, VirtualNode *vnode_2) {
        std::vector<int> vect1;
        std::vector<int> vect2;
        VirtualNode *v1 = vnode_1;
        VirtualNode *v2 = vnode_2;

        vect1.push_back(v1->vnode_id);
        while(!v1->isPseudoRootNode()){
            v1 = v1->getNodeUp();
            vect1.push_back(v1->vnode_id);
        }

        vect2.push_back(v2->vnode_id);
        while(!v2->isPseudoRootNode()){
            v2 = v2->getNodeUp();
            vect2.push_back(v2->vnode_id);
        }

        if(v1->vnode_id != v2->vnode_id){
            vect1.insert( vect1.end(), vect2.begin(), vect2.end() );
        }else{

            std::reverse(vect1.begin(), vect1.end());
            std::reverse(vect2.begin(), vect2.end());

            if(vect1.size() <= vect2.size()){
                int idx;
                for(idx=0;idx<vect1.size();idx++){
                    if(vect1.at(idx)!=vect2.at(idx)){
                        break;
                    }
                }
                vect1.insert( vect1.end(), vect2.begin()+idx, vect2.end() );

                return vect1;
            }else{
                int idx;
                for(idx=0;idx<vect2.size();idx++){
                    if(vect1.at(idx)!=vect2.at(idx)){
                        break;
                    }
                }
                vect2.insert( vect2.end(), vect1.begin()+idx, vect1.end() );

                return vect2;
            }

        }

        return vect1;
    }

    std::vector<int> Utree::computePathBetweenNodes(VirtualNode *vnode_1, VirtualNode *vnode_2) {

        std::vector<VirtualNode *> list_vnode_to_root;

        std::vector<VirtualNode *> path2root_1 = findPseudoRoot(vnode_1, false);
        std::vector<VirtualNode *> path2root_2 = findPseudoRoot(vnode_2, false);

        std::reverse(path2root_1.begin(), path2root_1.end());
        std::reverse(path2root_2.begin(), path2root_2.end());

        // Test map

        std::map<VirtualNode *, int> tmpMap;

        for (int i = 0; i < path2root_1.size(); i++) {
            tmpMap.insert(std::pair<VirtualNode *, int>(path2root_1.at(i), i));
        }

        for (int i = 0; i < path2root_2.size(); i++) {

            if (tmpMap.find(path2root_2.at(i)) == tmpMap.end()) {

                tmpMap.insert(std::pair<VirtualNode *, int>(path2root_2.at(i), tmpMap.size()));

            }
        }

        //auto newvector = std::vector<VirtualNode *>(tmpMap.size());
        auto newvector = std::vector<int>(tmpMap.size());
        // Create a map iterator and point to beginning of map
        //std::map<VirtualNode *,int>::iterator it = tmpMap.begin();

        // Iterate over the map using c++11 range based for loop
        for (std::pair<VirtualNode *, int> element : tmpMap) {

            // Accessing KEY from element
            VirtualNode *node = element.first;

            // Accessing VALUE from element.
            int order = element.second;

            newvector[order] = node->getVnode_id();

        }

        std::reverse(newvector.begin(), newvector.end());


        return newvector;
    }

    std::vector<VirtualNode *> Utree::_unique(std::vector<VirtualNode *> &list_nodes_n1, std::vector<VirtualNode *> &list_nodes_n2) {

        std::vector<VirtualNode *> list_nodes;
        VirtualNode *n1;
        VirtualNode *n2;

        while (list_nodes_n1.size() > 0 && list_nodes_n2.size() > 0) {
            n1 = list_nodes_n1.at(list_nodes_n1.size() - 1);
            n2 = list_nodes_n2.at(list_nodes_n2.size() - 1);
            if (n1 == n2) {
                list_nodes.push_back(n1);
                list_nodes_n1.pop_back();
                list_nodes_n2.pop_back();
            } else {
                break;
            }
        }

        while (list_nodes_n1.size() > 0) {
            n1 = list_nodes_n1.at(list_nodes_n1.size() - 1);
            list_nodes.push_back(n1);
            list_nodes_n1.pop_back();
        }

        while (list_nodes_n2.size() > 0) {
            n2 = list_nodes_n2.at(list_nodes_n2.size() - 1);
            list_nodes.push_back(n2);
            list_nodes_n2.pop_back();
        }

        return list_nodes;
    }

    std::vector<VirtualNode *> Utree::getPostOrderNodeList() {

        bool removeRoot = false;
        std::vector<VirtualNode *> rlist;
        if (rootnode->getNodeRight() == nullptr && rootnode->getNodeLeft() == nullptr) {
            addVirtualRootNode();
            removeRoot = true;
        }
        _getPostOrderNodeList(rlist, rootnode);


        if (removeRoot) {
            removeVirtualRootNode();
        }

        return rlist;
    }

    std::vector<VirtualNode *> Utree::getPostOrderNodeList(VirtualNode *startNode) {

        bool removeRoot = false;
        std::vector<VirtualNode *> rlist;

        if (rootnode->getNodeRight() == nullptr && rootnode->getNodeLeft() == nullptr) {
            addVirtualRootNode();
            removeRoot = true;
        }
        _getPostOrderNodeList(rlist, startNode);

        if (removeRoot) {
            removeVirtualRootNode();
        }
        return rlist;
    }

    std::vector<int> Utree::getPostOrderNodeList(int NodeID){
        bool removeRoot = false;
        std::vector<int> rlist;
        if (rootnode->getNodeRight() == nullptr && rootnode->getNodeLeft() == nullptr) {
            addVirtualRootNode();
            removeRoot = true;
        }

        _getPostOrderNodeList(rlist, NodeID);

        if (removeRoot) {
            removeVirtualRootNode();
        }
        return rlist;

    }

    void Utree::_getPostOrderNodeList(std::vector<int> &rlist, int nodeID){

        if (!getNode(nodeID)->isTerminalNode()) {

            _getPostOrderNodeList(rlist, getNode(nodeID)->getNodeLeft()->getVnode_id());
            _getPostOrderNodeList(rlist, getNode(nodeID)->getNodeRight()->getVnode_id());

        }

        rlist.push_back(nodeID);

    }


    void Utree::_getPostOrderNodeList(std::vector<VirtualNode *> &rlist, VirtualNode *node) {


        if (!node->isTerminalNode()) {

            _getPostOrderNodeList(rlist, node->getNodeLeft());
            _getPostOrderNodeList(rlist, node->getNodeRight());

        }

        rlist.push_back(node);

    }

    void Utree::computeTreeDepth() {
        for (auto &node : listVNodes) {
            computeNodeDepth(node);
        }
        initialized_treeDepth = true;
    }

    int Utree::getTreeDepthAtNode(VirtualNode *vnode) {
        int level = 0;

        if (initialized_treeDepth) {
            level = vnode->getNodeLevel();
        } else {
            LOG(ERROR) << "No depth on this tree! You should call setTreeDepth before accessing the level of node: " << vnode->getNodeName();
        }
        return level;
    }

    void Utree::computeNodeDepth(VirtualNode *vnode) {

        int level = (int) (findPseudoRoot(vnode)).size();

        vnode->setNodeLevel(level);

    }

    const std::map<int, VirtualNode *> &Utree::getNodeIdsMap() const {
        return utreeNodeIdsMap_;
    }


    void VirtualNode::resetNodeDirections(bool revertRotations) {

        if (vnode_rotated != NodeRotation::undef) {
            switch (vnode_rotated) {
                case NodeRotation::clockwise:
                    rotateCounterClockwise(revertRotations);
                    break;
                case NodeRotation::counterclockwise:
                    rotateClockwise(revertRotations);
                    break;
                default:
                    break;
            }
        }

    }

    void VirtualNode::_recursive_cw_rotation(VirtualNode *vnode, bool revertRotations) {

        bool continueTraverse = true;
        NodePosition curr_node_position = NodePosition::undef;

        if (revertRotations and vnode->vnode_rotated == NodeRotation::undef) { return; }

        if (revertRotations) {

            continueTraverse = true;

        } else {
            if (vnode->getNodeUp() != nullptr) {
                if (vnode->getNodeUp()->getNodeUp() == vnode) {
                    continueTraverse = false;
                }
            }
            curr_node_position = vnode->indexOf();
        }

        // rotate the pointers only if the node is not a leaf
        if (!vnode->isTerminalNode()) {

            VirtualNode *n_up;
            VirtualNode *curr_vn_up = vnode->getNodeUp();

            VirtualNodeUtils::rotateNodeClockwise(vnode);

            if (revertRotations) {
                n_up = vnode->getNodeUp();

                switch (n_up->vnode_rotated) {
                    case NodeRotation::clockwise:
                        curr_node_position = NodePosition::left;
                        break;
                    case NodeRotation::counterclockwise:
                        curr_node_position = NodePosition::right;
                        break;
                    default:
                        curr_node_position = NodePosition::undef;
                }

            } else {
                n_up = curr_vn_up;
            }

            if (continueTraverse) {
                if (n_up != nullptr) {
                    switch (curr_node_position) {
                        case NodePosition::left:
                            _recursive_ccw_rotation(n_up, revertRotations);
                            break;
                        case NodePosition::right:
                            _recursive_cw_rotation(n_up, revertRotations);
                            break;
                        case NodePosition::up:
                            _recursive_cw_rotation(n_up, revertRotations);
                        default:
                            return;
                    }
                }


            } else {

                //std::cout << "This node is the pseudoroot" << std::endl;

            }
        }

    }

    void VirtualNode::_recursive_ccw_rotation(VirtualNode *vnode, bool revertRotations) {

        bool continueTraverse = true;
        NodePosition curr_node_position = NodePosition::undef;

        if (revertRotations and vnode->vnode_rotated == NodeRotation::undef) { return; }

        if (revertRotations) {

            continueTraverse = true;

        } else {
            if (vnode->getNodeUp() != nullptr) {
                if (vnode->getNodeUp()->getNodeUp() == vnode) {

                    continueTraverse = false;
                }
            }
            curr_node_position = vnode->indexOf();
        }

        // rotate the pointers only if the node is not a leaf
        if (!vnode->isTerminalNode()) {

            VirtualNode *n_up;
            VirtualNode *curr_vn_up = vnode->getNodeUp();

            VirtualNodeUtils::rotateNodeCounterClockwise(vnode);

            if (revertRotations) {
                n_up = vnode->getNodeUp();

                switch (n_up->vnode_rotated) {
                    case NodeRotation::clockwise: //clockwise
                        curr_node_position = NodePosition::left; // force to rotate counterclockwise
                        break;
                    case NodeRotation::counterclockwise: //counterclockwise
                        curr_node_position = NodePosition::right; // force to rotate clockwise
                        break;
                    default:
                        curr_node_position = NodePosition::undef;
                }

            } else {
                n_up = curr_vn_up;
            }

            if (continueTraverse) {
                if (n_up != nullptr) {
                    switch (curr_node_position) {
                        case NodePosition::left:
                            _recursive_ccw_rotation(n_up, revertRotations);
                            break;
                        case NodePosition::right:
                            _recursive_cw_rotation(n_up, revertRotations);
                            break;
                        case NodePosition::up:
                            _recursive_ccw_rotation(n_up, revertRotations);
                        default:
                            return;
                    }
                }

            } else {

                //std::cout << "This node is the pseudoroot" << std::endl;

            }

        }

    }

    void VirtualNodeUtils::rotateNodeCounterClockwise(VirtualNode *vnode) {

        VirtualNode *curr_vn_up;
        VirtualNode *curr_vn_left;
        VirtualNode *curr_vn_right;

        curr_vn_up = vnode->getNodeUp();
        curr_vn_left = vnode->getNodeLeft();
        curr_vn_right = vnode->getNodeRight();

        // Revert pointers
        vnode->_setNodeUp(curr_vn_left);
        vnode->_setNodeLeft(curr_vn_right);
        vnode->_setNodeRight(curr_vn_up);

        // Store the information about the rotation
        if (vnode->vnode_rotated != NodeRotation::undef) {
            vnode->vnode_rotated = NodeRotation::undef;
        } else {
            vnode->vnode_rotated = NodeRotation::counterclockwise;
        }
    }


    void VirtualNodeUtils::rotateNodeClockwise(VirtualNode *vnode) {


        VirtualNode *curr_vn_up;
        VirtualNode *curr_vn_left;
        VirtualNode *curr_vn_right;

        curr_vn_up = vnode->getNodeUp();
        curr_vn_left = vnode->getNodeLeft();
        curr_vn_right = vnode->getNodeRight();

        vnode->_setNodeUp(curr_vn_right);
        vnode->_setNodeLeft(curr_vn_up);
        vnode->_setNodeRight(curr_vn_left);

        // Store the information about the rotation
        if (vnode->vnode_rotated != NodeRotation::undef) {
            vnode->vnode_rotated = NodeRotation::undef;
        } else {
            vnode->vnode_rotated = NodeRotation::clockwise;
        }

    }


    bool VirtualNode::isPseudoRootNode() {
        return getNodeUp()->getNodeUp() == this;
    }

    int VirtualNode::getNodeLevel() {

        return vnode_depth;
    }

    void VirtualNode::setNodeLevel(int level) {

        vnode_depth = level;

    }

    void VirtualNode::setBranchLength(double blength) {

        vnode_branchlength = blength;

    }

    int VirtualNode::getVnode_id() const {
        return vnode_id;
    }

    void VirtualNode::setVnode_id(int in_vnode_id) {
        VirtualNode::vnode_id = in_vnode_id;
    }

    VirtualNode *VirtualNode::getSiblingNode() {
        VirtualNode *siblingNode = nullptr;
        if (indexOf() == NodePosition::left) {
            siblingNode = getNodeUp()->getNodeRight();
        } else if (indexOf() == NodePosition::right) {
            siblingNode = getNodeUp()->getNodeLeft();
        } else if (isTerminalNode()) {
            LOG(WARNING) << "[tshlib::VirtualNode::getSiblingNode] Node " << getNodeName() << " has no siblings since it is a leaf.";
        } else {
            //siblingNode = getNodeUp();
        }
        return siblingNode;
    }

    VirtualNode::~VirtualNode() = default;

}
