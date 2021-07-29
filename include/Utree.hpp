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
 * @file Utree.hpp
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

#ifndef TSHLIB_UTREE_HPP
#define TSHLIB_UTREE_HPP

#include <string>
#include <vector>
#include <map>

#include "Utilities.hpp"

namespace tshlib {

/*
 * VirtualNode contains the virtual directions for traversing the tree. Each VirtualNode represents
 * either an internal node or a leaf.
 *
 *         [VN]
 *          ||
 *          ||
 *  +----------------+
 *  |       up       |
 *  |     /    \     |  Node of a real phylogenetic tree
 *  | left  -  right |
 *  +----------------+
 *    //          \\
 *   //            \\
 * [VN]            [VN]
 */

    class VirtualNode {

    protected:
        VirtualNode *vnode_up;                          // NodeUp - This is the pointer to the VirtualNode above
        VirtualNode *vnode_left;                        // NodeLeft  -  This is the pointer to the VirtualNode on the leftside
        VirtualNode *vnode_right;                       // NodeRight - This is the pointer to the VirtualNode on the rightside

    public:

        int vnode_id;                                   // Node ID - Useful in case of parallel independent executions
        std::string vnode_name;                         // Node Name - Useful in case of parallel independent executions
        double vnode_branchlength;                      // Branch length connecting the node the parent node
        int vnode_depth;                                // Depth level of the node in the tree
        bool vnode_leaf;                                // Flag: terminal node in the tree listVNodes
        int vnode_move_direction;                       // Int: This attribute is used to perform the correct rotation of the p-node w.r.t q-node.
        NodeRotation vnode_rotated;


        // Flag: if node was rotaded during a tree rearrangement move
        int vnode_seqid;                                // Needed to associate the seqID to the nodeID

        /*!
         *  Standard constructor
         */
        VirtualNode();

        VirtualNode(const VirtualNode &inNode){

            vnode_id = inNode.vnode_id;
            vnode_seqid = inNode.vnode_seqid;
            vnode_name = inNode.vnode_name;
            vnode_depth = inNode.vnode_depth;
            vnode_leaf = inNode.vnode_leaf;
            vnode_move_direction = inNode.vnode_move_direction;
            vnode_branchlength = inNode.vnode_branchlength;
            vnode_leaf = inNode.vnode_leaf;
            vnode_rotated = inNode.vnode_rotated;

            vnode_right = inNode.vnode_right;
            vnode_left = inNode.vnode_left;
            vnode_up = inNode.vnode_up;

        };

        VirtualNode &operator=(const VirtualNode &inNode) {


            vnode_id = inNode.vnode_id;
            vnode_seqid = inNode.vnode_seqid;
            vnode_name = inNode.vnode_name;
            vnode_depth = inNode.vnode_depth;
            vnode_leaf = inNode.vnode_leaf;
            vnode_move_direction = inNode.vnode_move_direction;
            vnode_branchlength = inNode.vnode_branchlength;
            vnode_leaf = inNode.vnode_leaf;
            vnode_rotated = inNode.vnode_rotated;

            vnode_right = inNode.vnode_right;
            vnode_left = inNode.vnode_left;
            vnode_up = inNode.vnode_up;

        };

        /*!
         * Virtual deconstructor
         */
        ~VirtualNode();

        /*!
         * @brief This function connects the current node to another one. It automatically performs a bidirectional connection
         * @param inVNode Target node to apply the connection to
         */
        void connectNode(VirtualNode *inVNode);


        void setBranchLength(double blength);


        /*!
         * @brief This function disconnect the current node from any other one above it. The function is bidirectional.
         */
        void disconnectNode();

        /*!
         * @brief This function perfomrms a one-step clockwise rotation of the virtualnode pointers
         */
        void rotateClockwise();

        void rotateClockwise(bool revertRotations);

        /*!
         * @brief This function perfomrms a one-step counter-clockwise rotation of the virtualnode pointers
         */
        void rotateCounterClockwise();

        void rotateCounterClockwise(bool revertRotations);

        /*!
         * @brief This function resets any rotation previously performed on the node
         */
        void resetNodeDirections(bool revertRotations);

        /*!
         * @brief The function prints the neighborhood of the node in the format <^nodeUp;nodeLeft;nodeRight >
         * @return std::string with the node neighborhood
         */
        std::string printNeighbours();

        void setNodeName(std::string s);

        const std::string getNodeName();

        VirtualNode *getNodeUp();

        VirtualNode *getNodeLeft();

        VirtualNode *getNodeRight();

        void _traverseVirtualNodeTree();

        double computeTotalTreeLength();

        VirtualNode *getSiblingNode();

        int getNodeLevel();

        void setNodeLevel(int level);

        void clearChildren();

        int getVnode_id() const;

        void setVnode_id(int vnode_id);

        NodeRotation getVnode_rotated() const;

        /*!
         * @brief This function returns true if the node is terminal
         * @return boolean value (true or false)
         */
        bool isTerminalNode();

        /*!
         * @brief This function returns true if the node is root
         * @return boolean value (true or false)
         */
        bool isRootNode();

        /*!
         * @brief This function returns true if the node is pseudoroot
         * @return boolean value (true or false)
         */
        bool isPseudoRootNode();

        /*!
         * @brief This function return the index of node as seen from the parent immediate above it.
         * @return
         */
        NodePosition indexOf();

        /*!
         * @brief The function checks if the current node is a parent of another node using a recursive structure
         * @param inVNode VirtualNode pointer
         * @return False if the node passed is not parent of the current one, True otherwise
         * @deprecated
         */
        bool isParent(VirtualNode *inVNode);

        void _setNodeRight(VirtualNode *inVNode);

        void _setNodeLeft(VirtualNode *inVNode);

        void _setNodeUp(VirtualNode *inVNode);

        void _bidirectionalUpwardConnection(VirtualNode *inNode);

        NodeRotation getNodeRotation() const { return vnode_rotated; }

    protected:

        void _oneway_connectNode(VirtualNode *inVNode);

        void _recursive_cw_rotation(VirtualNode *vnode, bool revertRotations);

        void _recursive_ccw_rotation(VirtualNode *vnode, bool revertRotations);


    };


    class Utree {
    private:
        bool initialized_treeDepth;

    public:

        std::vector<VirtualNode *> listVNodes;
        std::vector<VirtualNode *> startVNodes;
        std::map<int, VirtualNode *> utreeNodeIdsMap_;

        VirtualNode *rootnode;

        Utree();

        Utree(const Utree &inTree) {

            initialized_treeDepth = false;
            // copy nodes of the tree
            for(auto &_node__intree:inTree.listVNodes){

                // Deep copy node attributes
                auto _node__utree = new VirtualNode(*_node__intree);
                _node__utree->_setNodeLeft(nullptr);
                _node__utree->_setNodeRight(nullptr);
                _node__utree->_setNodeUp(nullptr);

                // add copied node to the tree
                addMember(_node__utree, false);

            }

            // update nodes connections with new pointers
            for(auto &_node__utree:listVNodes){

                VirtualNode &_node__intree = inTree.getNode(_node__utree->getVnode_id());

                if(_node__intree.getNodeLeft())
                    _node__utree->_setNodeLeft(getNode(_node__intree.getNodeLeft()->getVnode_id()));

                if(_node__intree.getNodeRight())
                    _node__utree->_setNodeRight(getNode(_node__intree.getNodeRight()->getVnode_id()));

                if(_node__intree.getNodeUp())
                    _node__utree->_setNodeUp(getNode(_node__intree.getNodeUp()->getVnode_id()));

            }

            for(auto &_startnode__intree:inTree.startVNodes){
                startVNodes.push_back(getNode(_startnode__intree->getVnode_id()));
            }

            // Add rootnode
            rootnode = new VirtualNode(*inTree.rootnode);

            // Add pointer reference to the utree map
            utreeNodeIdsMap_[rootnode->getVnode_id()] = rootnode;

        }

        Utree &operator=(const Utree &inTree) {};

        ~Utree();

        /*!
         * @brief This function finds the pseudoroot traversing the tree from starting node until a bidirectional connection is found
         * @param inVNode starting node
         * @return std::vector of VirtualNode pointers from the starting point until the pseudoroot
         */
        std::vector<VirtualNode *> findPseudoRoot(VirtualNode *inVNode, bool fixPseudoRootOnNextSubtree = false);

        /*!
         * @brief This function adds a complete VirtualNode to the listVNodes attribute of the class Utree
         * @param inVNode       Pointer to the VirtualNode to add
         * @param isStartNode
         */
        void addMember(VirtualNode *inVNode, bool isStartNode = false);

        //=======================
        void myAddRoot();
        void myRemoveRoot();
        //=======================

        /*!
         * @brief
         * @param showInternalNodeNames
         * @return
         */
        std::string printTreeNewick(bool showInternalNodeNames, bool updateStartNodes = true);

        /*!
         * @brief
         * @return
         */
        int getMaxNodeDistance();

        //===========================
        int mygetMaxNodeDistance();
        //===========================


        /*!
         * @brief
         * @param outfilepath
         */
        void saveTreeOnFile(std::string outfilepath);

        /*!
         * @brief
         */
        void printAllNodesNeighbors();

        /*!
         * @brief This function create a virtual root on the utree object. It breaks the link between two pseudoroot virtualnodes
         */
        void addVirtualRootNode();

        /*!
         * @brief This function remove the virtual root node added with the companion function addVirtualRootNode()
         */
        void removeVirtualRootNode();

        void _testReachingPseudoRoot();

        double computeTotalTreeLength();

        void _printUtree();

        std::vector<int> computePathBetweenNodes2(VirtualNode *vnode_1, VirtualNode *vnode_2);

        std::vector<int> computePathBetweenNodes(VirtualNode *vnode_1, VirtualNode *vnode_2);

        std::vector<VirtualNode *> _unique(std::vector<VirtualNode *> &list_nodes_n1, std::vector<VirtualNode *> &list_nodes_n2);

        std::vector<VirtualNode *> getPostOrderNodeList();

        std::vector<VirtualNode *> getPostOrderNodeList(VirtualNode *startNode);

        std::vector<int> getPostOrderNodeList(int NodeID);

        void _getPostOrderNodeList(std::vector<VirtualNode *> &rlist, VirtualNode *node);

        void _getPostOrderNodeList(std::vector<int> &rlist, int nodeID) ;

        void computeTreeDepth();

        int getTreeDepthAtNode(VirtualNode *vnode);

        void computeNodeDepth(VirtualNode *vnode);

        const std::map<int, VirtualNode *> &getNodeIdsMap() const;

        VirtualNode *getNode(int nodeID) {
            if(nodeID == -1)
                return rootnode;

            return utreeNodeIdsMap_[nodeID];
        }

        VirtualNode & getNode(int nodeID) const {
            if(nodeID == -1)
                return (*rootnode);

            return (*utreeNodeIdsMap_.at(nodeID));
        }

        std::vector<int> getChildrenId(int nodeID) const;


    protected:
        std::string _recursiveFormatNewick(VirtualNode *vnode, bool showInternalNodeNames);


        void _updateStartNodes();

    private:

    };


    namespace VirtualNodeUtils {

        void rotateNodeClockwise(VirtualNode *vnode);

        void rotateNodeCounterClockwise(VirtualNode *vnode);
    }


    namespace UtreeUtils {
        using namespace tshlib;

        VirtualNode *getPseudoRoot(VirtualNode *vn);

    }
}
#endif //TSHLIB_UTREE_HPP
