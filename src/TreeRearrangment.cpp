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
#include <string>
#include <iostream>
#include <numeric>
#include <limits>
#include <iomanip>
#include <iterator>
#include <chrono>
#include <algorithm>
#include <glog/logging.h>
#include <map>
#include "TreeRearrangment.hpp"

#include <algorithm>
#include <random>

namespace tshlib {
    TreeRearrangment::TreeRearrangment() = default;

    TreeRearrangment::~TreeRearrangment() {

        for (std::vector<Move *>::reverse_iterator i = trMoveSet.rbegin(); i < trMoveSet.rend(); i++) {
            Move *move = *i;
            delete move;
        }
        std::vector<Move *>().swap(trMoveSet);


    };

    std::string TreeRearrangment::getStrategy() const {
        std::string rtToken;
        switch (trStrategy) {
            case TreeSearchHeuristics::swap:
                rtToken = "Swap";
                break;
            case TreeSearchHeuristics::phyml:
                rtToken = "phyML";
                break;
            case TreeSearchHeuristics::mixed:
                rtToken = "Mixed(Swap+phyML)";
                break;
            case TreeSearchHeuristics::nosearch:
                rtToken = "no-search";
                break;
        }
        return rtToken;
    }

    std::string TreeRearrangment::getRearrangmentCoverageDescription() const {
        std::string rtToken;
        switch (trTreeCoverage) {
            case TreeRearrangmentOperations::classic_NNI:
                rtToken = "NNI-like";
                break;
            case TreeRearrangmentOperations::classic_SPR:
                rtToken = "SPR-like";
                break;
            case TreeRearrangmentOperations::classic_TBR:
                rtToken = "TBR-like";
                break;
            case TreeRearrangmentOperations::classic_Mixed:
                rtToken = "Complete";
                break;
        }
        return rtToken;
    }

    void TreeRearrangment::getNodesInRadiusDown(VirtualNode *node_source,VirtualNode *node_target, int radius_min, int radius_curr, int radius_max, bool includeSelf,
                                                MoveDirections direction, bool allowDuplicatedMoves) {

        VirtualNode *node = node_target;

        if (!includeSelf) {
            includeSelf = true;
        } else {
            if (radius_curr <= (radius_max - radius_min) && radius_curr >= 0) {
                //addMove(radius_max - radius_curr, node, direction, allowDuplicatedMoves);

                addMove(node_source,node, direction,radius_max - radius_curr);
            }
        }

        if (radius_curr < 0) {
            return;
        }

        if (!node->isTerminalNode()) {
            radius_curr--;
            // Direction is required to know where to look at in the tree when the move is performed
            getNodesInRadiusDown(node_source,node->getNodeLeft(), radius_min, radius_curr, radius_max, includeSelf, direction, allowDuplicatedMoves);
            getNodesInRadiusDown(node_source,node->getNodeRight(), radius_min, radius_curr, radius_max, includeSelf, direction, allowDuplicatedMoves);
        }

    }

    void TreeRearrangment::getNodesInRadiusUp(VirtualNode *node_source,VirtualNode *node_target, int radius_min, int radius_curr, int radius_max, NodePosition traverse_direction,
                                         bool allowDuplicatedMoves) {

        VirtualNode *vnode;

        NodePosition idx;
        MoveDirections moving_direction;

        vnode = node_target;

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

            addMove(node_source, vnode, moving_direction, radius_max - radius_curr);
        }

        if (radius_curr > 0) {
            if (traverse_direction == NodePosition::left) {
                radius_curr--;
                if (vnode->getNodeUp() != nullptr) {

                    idx = vnode->indexOf();
                    getNodesInRadiusUp(node_source,vnode->getNodeUp(), radius_min, radius_curr, radius_max, idx, allowDuplicatedMoves);
                }

                getNodesInRadiusDown(node_source,vnode->getNodeRight(), radius_min, radius_curr, radius_max, true, MoveDirections::up, allowDuplicatedMoves);

            } else if (traverse_direction == NodePosition::right) {
                radius_curr--;
                if (vnode->getNodeUp() != nullptr) {

                    idx = vnode->indexOf();
                    getNodesInRadiusUp(node_source,vnode->getNodeUp(), radius_min, radius_curr, radius_max, idx, allowDuplicatedMoves);
                }

                getNodesInRadiusDown(node_source,vnode->getNodeLeft(), radius_min, radius_curr, radius_max, true, MoveDirections::up, allowDuplicatedMoves);

            } else if (traverse_direction == NodePosition::up) {

                // If the traverse_direction is 2, we are moving across the pseudoroot, therefore no need to decrease the radious
                if (vnode->getNodeUp() != nullptr) {

                    // Get the nodes in the radius from the node we reached after crossing the pseudoroot
                    getNodesInRadiusDown(node_source,vnode, radius_min, radius_curr - 1, radius_max - 1, false, MoveDirections::up, allowDuplicatedMoves);

                }

            }
        }
    }

    Move * TreeRearrangment::setMove(tshlib::Move *m,int idx) {
        trMoveSet.at(idx) = m;
    }

    void TreeRearrangment::sortData(std::vector<int> &array,std::vector<int> &indeces){

        double swap;

        int idx;

        for (int i=0;i<array.size()-1;i++){
            for (int j=0;j<array.size()-i-1;j++){
                if (array[j] < array[j+1]){
                    swap=array[j];
                    array[j]=array[j+1];
                    array[j+1]=swap;

                    idx=indeces[j];
                    indeces[j]=indeces[j+1];
                    indeces[j+1]=idx;
                }
            }
        }

    }

    std::vector<int> TreeRearrangment::sortByDepth(std::vector<int> &path,tshlib::Utree *utree){

        VirtualNode *node = nullptr;

        std::vector<int> path_len(path.size());
        std::vector<int> indeces(path.size());
        std::vector<int> path_sorted(path.size());
        for(int node_i=0;node_i<path.size();node_i++){
            node = utree->getNode(path.at(node_i));
            path_len.at(node_i)=0;
            indeces.at(node_i)=node_i;
            path_sorted.at(node_i)=-1;
            while(!node->isPseudoRootNode()){
                path_len.at(node_i)=path_len.at(node_i)+1;
                node=node->getNodeUp();
            }
            path_len.at(node_i)=path_len.at(node_i)+1;
        }

        sortData(path_len,indeces);

        for(int node_i=0;node_i<path.size();node_i++){
            path_sorted.at(node_i)=path.at(indeces.at(node_i));
        }

        return path_sorted;
    }

    std::vector<int> TreeRearrangment::myPathBetweenNodes(tshlib::VirtualNode *_node__source,
                                                          tshlib::VirtualNode *_node__target,
                                                          tshlib::Utree *utree) {

        std::vector<int> path_between_nodes;
        int idx;
        VirtualNode *node = nullptr;

        std::vector<int> path_source;
        node = _node__source;
        while(!node->isPseudoRootNode()){
            path_source.push_back(node->vnode_id);
            node=node->getNodeUp();
        }
        path_source.push_back(node->vnode_id);

        std::reverse(path_source.begin(),path_source.end());

        std::vector<int> path_target;
        node = _node__target;
        while(!node->isPseudoRootNode()){
            path_target.push_back(node->vnode_id);
            node=node->getNodeUp();
        }
        path_target.push_back(node->vnode_id);

        std::reverse(path_target.begin(),path_target.end());

        if(path_source.size()>path_target.size()){

            for(idx=0;idx<path_target.size();idx++){
                if(path_target.at(idx)==path_source.at(idx)){
                    path_between_nodes.push_back(path_target.at(idx));
                }else{
                    break;
                }
            }
            for(int i=idx;i<path_target.size();i++){
                path_between_nodes.push_back(path_target.at(i));
            }
            for(int i=idx;i<path_source.size();i++){
                path_between_nodes.push_back(path_source.at(i));
            }

        }else{

            for(idx=0;idx<path_source.size();idx++){
                if(path_target.at(idx)==path_source.at(idx)){
                    path_between_nodes.push_back(path_target.at(idx));
                }else{
                    break;
                }
            }
            for(int i=idx;i<path_source.size();i++){
                path_between_nodes.push_back(path_source.at(i));
            }
            for(int i=idx;i<path_target.size();i++){
                path_between_nodes.push_back(path_target.at(i));
            }

        }

        return sortByDepth(path_between_nodes,utree);;
    }

    void TreeRearrangment::mydefineMoves(tshlib::Utree *utree,std::default_random_engine &random_engine,bool shuffle) {

        utree->removeVirtualRootNode();

        int num_nodes = utree->listVNodes.size();
        int up = 1;
        int up_lr = 2;

        std::vector<std::vector <int> > Mat(num_nodes);
        std::vector<std::vector <int> > Dir(num_nodes);
        for(int i=0;i<num_nodes;i++){
            Mat[i].resize(num_nodes);
            Dir[i].resize(num_nodes);
            std::fill(Mat[i].begin(),Mat[i].end(),0);
            std::fill(Dir[i].begin(),Dir[i].end(),0);
        }

        std::vector<std::vector<int>> paths_to_root;

        for (int i=0;i < num_nodes; i++) {
            std::vector<int> v;
            VirtualNode *node = utree->listVNodes.at(i);
            while(!node->isPseudoRootNode()){
                v.push_back(node->vnode_id);
                node=node->getNodeUp();
            }
            v.push_back(node->vnode_id);
            paths_to_root.push_back(v);
        }

        int ix,iy;
        int dim1, dim2;
        int id1,id2;
        int count = 0;
        int d;
        int k;
        bool is_same;
        std::vector<int> *v1 = nullptr;
        std::vector<int> *v2 = nullptr;
        for(int i=0;i<num_nodes-1;i++){
            v1 = &(paths_to_root.at(i));
            dim1 = v1->size();
            for(int j=i+1;j<num_nodes;j++){
                v2 = &(paths_to_root.at(j));
                dim2 = v2->size();
                id1 = v1->at(dim1-1);
                id2 = v2->at(dim2-1);

                if(id1!=id2){
                    d = dim1 + dim2 - 1;
                    Mat[i][j] = d;
                    Dir[i][j] = up;
                    if(d>2){
                        count++;
                    }
                }else{

                    if(dim1>dim2){

                        // source:i ; target:j
                        k = 0;
                        ix=dim1;
                        iy=dim2;
                        is_same = true;
                        for(k=0;k<dim2;k++){
                            ix--;
                            iy--;
                            if(v1->at(ix) != v2->at(iy)){
                                is_same = false;
                                break;
                            }
                        }

                        if(is_same){
                            Dir[i][j] = up_lr;
                        }else{
                            Dir[i][j] = up;
                        }

                        d = dim1 -k + dim2 -k;
                        Mat[i][j] = d;

                        if(d>2){
                            count++;
                        }

                    }else{
                        // source:j ; target:i
                        k = 0;
                        ix=dim1;
                        iy=dim2;
                        is_same = true;
                        for(k=0;k<dim1;k++){
                            ix--;
                            iy--;
                            if(v1->at(ix) != v2->at(iy)){
                                is_same = false;
                                break;
                            }
                        }

                        if(is_same){
                            Dir[j][i] = up_lr;
                        }else{
                            Dir[j][i] = up;
                        }

                        d = dim1 -k + dim2 -k;
                        Mat[j][i] = d;

                        if(d>2){
                            count++;
                        }

                    }

                }
            }
        }

        trMoveSet.resize(count);

        std::vector<int> indeces;
        if(shuffle){
            indeces.resize(count);
            for(int i=0;i<count;i++){
                indeces.at(i)=i;
            }

            std::shuffle(std::begin(indeces), std::end(indeces), random_engine);
        }

        int id=0;
        for(int i=0;i<num_nodes;i++) {
            for (int j = 0; j < num_nodes; j++) {
                if(Mat[i][j]==3){

                    Move *move = new Move;
                    move->initMove();
                    move->setSourceNode(utree->listVNodes.at(i)->vnode_id);
                    move->setTargetNode(utree->listVNodes.at(j)->vnode_id);
                    if(Dir[i][j]==up){
                        move->setDirection(tshlib::MoveDirections::up);
                    }else if(Dir[i][j]==up_lr){
                        move->setDirection(tshlib::MoveDirections::up_left);
                    }
                    move->setRadius(Mat[i][j]);
                    move->moveClassDescription_ = "NNI";
                    move->moveUID_ = id;

                    if(shuffle){
                        trMoveSet.at(indeces.at(id)) = move;
                    }else{
                        trMoveSet.at(id) = move;
                    }

                    id++;

                }else if(Mat[i][j]>3){

                    Move *move = new Move;
                    move->initMove();
                    move->setSourceNode(utree->listVNodes.at(i)->vnode_id);
                    move->setTargetNode(utree->listVNodes.at(j)->vnode_id);
                    if(Dir[i][j]==up){
                        move->setDirection(tshlib::MoveDirections::up);
                    }else if(Dir[i][j]==up_lr){
                        move->setDirection(tshlib::MoveDirections::up_left);
                    }
                    move->setRadius(Mat[i][j]);
                    move->moveClassDescription_ = "SPR";
                    move->moveUID_ = id;

                    if(shuffle){
                        trMoveSet.at(indeces.at(id)) = move;
                    }else{
                        trMoveSet.at(id) = move;
                    }

                    id++;
                }
            }
        }

        trCandidateMovesFound_ = trMoveSet.size();

        utree->addVirtualRootNode();

    }

    void TreeRearrangment::mydefineMoves_new(tshlib::Utree *utree,int source_id) {

        //utree->removeVirtualRootNode();
        utree->myRemoveRoot();

        int num_nodes = utree->listVNodes.size();
        int up = 1;
        int up_lr = 2;

        std::vector <int> Mat(num_nodes);
        std::vector <int> Dir(num_nodes);
        std::vector <int> Rev_S_T(num_nodes);
        std::fill(Mat.begin(),Mat.end(),0);
        std::fill(Dir.begin(),Dir.end(),0);
        std::fill(Rev_S_T.begin(),Rev_S_T.end(),0);

        std::vector<std::vector<int>> paths_to_root(num_nodes);
        //std::vector<int> paths_to_root_source;
        for (int i=0;i < num_nodes; i++) {
            std::vector<int> v;

            VirtualNode *node = utree->listVNodes.at(i);
            int node_id = node->vnode_id;

            while(!node->isPseudoRootNode()){
                v.push_back(node->vnode_id);
                node=node->getNodeUp();
            }
            v.push_back(node->vnode_id);
            //if(node_id!=source_id){
            paths_to_root.at(node_id)=v;
            //}else{
            //    paths_to_root_source = v;
            //}
        }

        int ix,iy;
        int dim1, dim2;
        int id1,id2;
        int count = 0;
        int d;
        int k;
        bool is_same;
        std::vector<int> *v1 = nullptr;
        std::vector<int> *v2 = nullptr;
        //for(int i=0;i<num_nodes-1;i++){
        v1 = &(paths_to_root.at(source_id));;
        dim1 = v1->size();
        for(int target_id=0; target_id < paths_to_root.size(); target_id++){
            if(source_id!=target_id){

                v2 = &(paths_to_root.at(target_id));
                dim2 = v2->size();
                id1 = v1->at(dim1-1);
                id2 = v2->at(dim2-1);

                if(id1!=id2){
                    d = dim1 + dim2 - 1;
                    Mat[target_id] = d;
                    Dir[target_id] = up;
                    if(d>2){
                        count++;
                    }
                }else{

                    if(dim1>dim2){

                        // source:i ; target:target_id
                        k = 0;
                        ix=dim1;
                        iy=dim2;
                        is_same = true;
                        for(k=0;k<dim2;k++){
                            ix--;
                            iy--;
                            if(v1->at(ix) != v2->at(iy)){
                                is_same = false;
                                break;
                            }
                        }

                        if(is_same){
                            Dir[target_id] = up_lr;
                        }else{
                            Dir[target_id] = up;
                        }

                        d = dim1 -k + dim2 -k;
                        Mat[target_id] = d;

                        if(d>2){
                            count++;
                        }

                    }else{
                        // source:target_id ; target:i
                        k = 0;
                        ix=dim1;
                        iy=dim2;
                        is_same = true;
                        for(k=0;k<dim1;k++){
                            ix--;
                            iy--;
                            if(v1->at(ix) != v2->at(iy)){
                                is_same = false;
                                break;
                            }
                        }

                        Rev_S_T[target_id] = 1;


                        if(is_same){
                            Dir[target_id] = up_lr;
                        }else{
                            Dir[target_id] = up;
                        }

                        d = dim1 -k + dim2 -k;
                        Mat[target_id] = d;

                        if(d>2){
                            count++;
                        }

                    }

                }
            }
        }

        trMoveSet.resize(count);

        int id=0;
        //for(int i=0;i<num_nodes;i++) {
        for (int target_id = 0; target_id < num_nodes; target_id++) {
            //if(target_id != source_id){

            if(Mat[target_id] == 3){

                Move *move = new Move;
                move->initMove();
                if(Rev_S_T[target_id]==0){
                    move->setSourceNode(source_id);
                    move->setTargetNode(target_id);
                }else{
                    move->setSourceNode(target_id);
                    move->setTargetNode(source_id);
                }
                if(Dir[target_id] == up){
                    move->setDirection(tshlib::MoveDirections::up);
                }else if(Dir[target_id] == up_lr){
                    move->setDirection(tshlib::MoveDirections::up_left);
                }
                move->setRadius(Mat[target_id]);
                move->moveClassDescription_ = "NNI";
                move->moveUID_ = id;

                trMoveSet.at(id) = move;

                id++;

            }else if(Mat[target_id] > 3){

                Move *move = new Move;
                move->initMove();
                if(Rev_S_T[target_id]==0){
                    move->setSourceNode(source_id);
                    move->setTargetNode(target_id);
                }else{
                    move->setSourceNode(target_id);
                    move->setTargetNode(source_id);
                }
                if(Dir[target_id] == up){
                    move->setDirection(tshlib::MoveDirections::up);
                }else if(Dir[target_id] == up_lr){
                    move->setDirection(tshlib::MoveDirections::up_left);
                }
                move->setRadius(Mat[target_id]);
                move->moveClassDescription_ = "SPR";
                move->moveUID_ = id;

                trMoveSet.at(id) = move;

                id++;
            }
            //}
        }

        trCandidateMovesFound_ = trMoveSet.size();

        //utree->addVirtualRootNode();
        utree->myAddRoot();

    }

    void TreeRearrangment::defineMoves(VirtualNode *sourceNode,bool includeSelf, bool allowDuplicatedMoves) {
        // Flag the nodes according to their position on the tree (left or right or above the source node -- p node).
        // For each node within the radius extremities, define a move and add it to TreeRearrangment.
        // Start from the children of the current starting node (if any)

        trCandidateMovesFound_ = 0;

        if (!sourceNode->isTerminalNode()) {

//            getNodesInRadiusDown(sourceNode->getNodeLeft(),
//                                 trSearchRadius_min,
//                                 trSearchRadius_max - 1,
//                                 trSearchRadius_max,
//                                 includeSelf,
//                                 MoveDirections::down_left,
//                                 allowDuplicatedMoves);
            getNodesInRadiusDown(sourceNode,
                                 sourceNode->getNodeLeft(),
                                 trSearchRadius_min,
                                 trSearchRadius_max - 1,
                                 trSearchRadius_max,
                                 includeSelf,
                                 MoveDirections::up_left,
                                 allowDuplicatedMoves);
//            getNodesInRadiusDown(sourceNode->getNodeRight(),
//                                 trSearchRadius_min,
//                                 trSearchRadius_max - 1,
//                                 trSearchRadius_max,
//                                 includeSelf,
//                                 MoveDirections::down_right,
//                                 allowDuplicatedMoves);
            getNodesInRadiusDown(sourceNode,
                                 sourceNode->getNodeRight(),
                                 trSearchRadius_min,
                                 trSearchRadius_max - 1,
                                 trSearchRadius_max,
                                 includeSelf,
                                 MoveDirections::up_right,
                                 allowDuplicatedMoves);
        }
        // If the node is a leaf, then go up
        if (nullptr != sourceNode->getNodeUp()) {

            getNodesInRadiusUp(sourceNode,
                               sourceNode->getNodeUp(),
                               trSearchRadius_min,
                               trSearchRadius_max - 1,
                               trSearchRadius_max,
                               sourceNode->indexOf(),
                               allowDuplicatedMoves);
        }


        //VLOG(1) << "[TSH Cycle]  Found " << trCandidateMovesFound_ << " candidate moves for node " << utree_->getNode(trSourceNode_)->getNodeName();
    }
    /*
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
    */
    /*
    void TreeRearrangment::addMove(VirtualNode *sourceNode,VirtualNode *targetNode, MoveDirections moveDirection,int radius) {

        bool allowDuplicatedMoves = false;

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
            //moveInstance->setSourceNode(trSourceNode_);
            moveInstance->setDirection(moveDirection);
            moveInstance->setRadius(radius);
            //moveInstance->setTargetNode(targetNode->getVnode_id());
            moveInstance->setSourceNode(targetNode->vnode_id);
            moveInstance->setTargetNode(sourceNode->vnode_id);
            moveInstance->setClass(ts_strategy, (sourceNode->isPseudoRootNode() || targetNode->isPseudoRootNode()));


            bool storeMove = true;

            if (!allowDuplicatedMoves) {
                for (auto &query:trMoveSet) {

                    if (query->getTargetNode() == targetNode->getVnode_id()
                        && query->getSourceNode() == sourceNode->vnode_id
                        && query->getMoveStrategy() == ts_strategy) {

                        storeMove = false;
                    }

                    if (targetNode->getVnode_id() == query->getSourceNode()
                        && sourceNode->vnode_id == query->getTargetNode()
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
    */

    void TreeRearrangment::addMove(VirtualNode *sourceNode,VirtualNode *targetNode, MoveDirections moveDirection,int radius) {

        std::vector<TreeSearchHeuristics> strategies;

//        if (trStrategy == TreeSearchHeuristics::mixed) {
//            strategies.push_back(TreeSearchHeuristics::swap);
//            strategies.push_back(TreeSearchHeuristics::phyml);
//        } else {
            strategies.push_back(trStrategy);
        //}

        for (auto &ts_strategy:strategies) {

            // Skip SPR-like moves if the radius is insufficient
            if (ts_strategy == TreeSearchHeuristics::phyml && radius < 4) continue;

            auto moveInstance = new Move;
            moveInstance->initMove();

            //moveInstance->setSourceNode(trSourceNode_);
            //moveInstance->setTargetNode(targetNode->getVnode_id());
            moveInstance->setSourceNode(targetNode->vnode_id);
            moveInstance->setTargetNode(sourceNode->vnode_id);

            moveInstance->setDirection(moveDirection);

            moveInstance->setRadius(radius);
            moveInstance->setClass(ts_strategy, (sourceNode->isPseudoRootNode() || targetNode->isPseudoRootNode()));


            //bool storeMove = true;

//            if (!allowDuplicatedMoves) {
//                for (auto &query:trMoveSet) {
//
//                    if (query->getTargetNode() == targetNode->getVnode_id()
//                        && query->getSourceNode() == trSourceNode_
//                        && query->getMoveStrategy() == ts_strategy) {
//
//                        storeMove = false;
//                    }
//
//                    if (targetNode->getVnode_id() == query->getSourceNode()
//                        && trSourceNode_ == query->getTargetNode()
//                        && query->getMoveStrategy() == ts_strategy) {
//                        storeMove = false;
//                    }
//
//                }
//            }

            //if (storeMove) {
            moveInstance->moveUID_ = (int) trMoveSet.size();
            trMoveSet.push_back(moveInstance);
            trCandidateMovesFound_++;
//            } else {
//                delete moveInstance;
//            }

        }
    }

    void TreeRearrangment::printMoves() {

        VLOG(2) << "[set " << trUID_ << "] " << getStrategy() << " strategy" << std::endl;
        VLOG(2) << "[class]\t(P\t; Q)" << std::endl;
//        for (auto &nmove: trMoveSet) {
//            VLOG(2) << "[" << nmove->moveClassDescription_ << "]\t" << nmove->moveRadius_
//                    << "\t(" << utree_->getNode(trSourceNode_)->vnode_name << "\t; "
//                    << utree_->getNode(nmove->getTargetNode())->vnode_name << ")\t"
//                    << static_cast<int>(nmove->moveDirection_) << std::endl;
//        }

    }

    unsigned long TreeRearrangment::getNumberOfMoves() {

        return trMoveSet.size();
    }


    //================================================================================
    bool TreeRearrangment::checkVectors(std::vector<int> &v1,std::vector<int> &v2){

        bool flag;

        if(v1.size() != v2.size()){
            return false;
        }

        for(int i=0;i<v1.size();i++){
            flag=false;
            for(int j=0;j<v2.size();j++){
                if(v1.at(i)==v2.at(j)){
                    flag=true;
                    break;
                }
            }
            if(!flag){
                return false;
            }
        }

        return true;
    }

    void TreeRearrangment::myCheckTree_rec(VirtualNode *node) {

        if(node->isTerminalNode()){
            return;
        }else{
            if(!node->getNodeLeft()){
                LOG(FATAL) << "ERROR";
            }
            if(node->getNodeLeft()->getNodeUp()!=node){
                LOG(FATAL) << "ERROR";
            }
            if(!node->getNodeRight()){
                LOG(FATAL) << "ERROR";
            }
            if(node->getNodeRight()->getNodeUp()!=node){
                LOG(FATAL) << "ERROR";
            }
            myCheckTree_rec(node->getNodeLeft());
            myCheckTree_rec(node->getNodeRight());
        }

    }

    void TreeRearrangment::myPrintTree(VirtualNode *node,std::vector<VirtualNode *> &startnodes) {

        if(node->isPseudoRootNode()){
            std::cout<<std::setprecision(18)<<node->vnode_name<<" |"<<node->vnode_branchlength<<"| "<<" ==== "<<node->getNodeUp()->vnode_name<<std::endl;
            std::cout<<std::setprecision(18)<<startnodes.at(0)->vnode_name<<" == [0] ==== [1] == "<<startnodes.at(1)->vnode_name<<std::endl;
        }
        if(node->isTerminalNode()) {
            std::cout<<std::setprecision(18)<<node->vnode_name<<" {"<<node->vnode_id<<"} "<<" |"<<node->vnode_branchlength<<"| "<<" ["<<node->getNodeUp()->vnode_name<<"] "<<std::endl;
        }else{
            std::cout<<std::setprecision(18)<<node->vnode_name<<" {"<<node->vnode_id<<"} "<<" |"<<node->vnode_branchlength<<"| "<<" : ("<<node->getNodeLeft()->vnode_name<<","<<node->getNodeRight()->vnode_name<<")"<<"  ["<<node->getNodeUp()->vnode_name<<"]"<<std::endl;
            myPrintTree(node->getNodeLeft(),startnodes);
            myPrintTree(node->getNodeRight(),startnodes);
        }

    }

    void TreeRearrangment::myCheckTreeCountLeaves_rec(VirtualNode *node,int &n){

        if(node->isTerminalNode()){
            n++;
        }else{
            myCheckTreeCountLeaves_rec(node->getNodeLeft(),n);
            myCheckTreeCountLeaves_rec(node->getNodeRight(),n);
        }

    }

    bool TreeRearrangment::myCheckTreeCountLeaves(tshlib::Utree *utree){

        int tot;
        int num_leaves_L = 0;
        int num_leaves_R = 0;
        myCheckTreeCountLeaves_rec(utree->startVNodes.at(0),num_leaves_L);
        myCheckTreeCountLeaves_rec(utree->startVNodes.at(1),num_leaves_R);
        tot=num_leaves_L+num_leaves_R;
        if( tot*2-2 != utree->listVNodes.size() ){
            return false;
        }

        return true;
    }

    void TreeRearrangment::myCheckTree(VirtualNode *root1,VirtualNode *root2) {

        // check the roots
        if(root1->getNodeUp() != root2){
            LOG(FATAL) << "ERROR";
        }
        if(root2->getNodeUp() != root1){
            LOG(FATAL) << "ERROR";
        }

        // check that each leaf points to the parent
        myCheckTree_rec(root1);
        myCheckTree_rec(root2);
    }

    void TreeRearrangment::myCheckTreeBranchLength(tshlib::Utree *utree) {

        for(int i=0;i<utree->listVNodes.size();i++){
            if(abs(utree->listVNodes.at(i)->vnode_branchlength)<0.000001){
                LOG(FATAL) << "ERROR: myCheckTreeBranchLength";
                exit(EXIT_FAILURE);
            }
        }

    }

    void TreeRearrangment::myCheckTreeNodeName(tshlib::Utree *utree) {

        for(auto &node: utree->listVNodes){
            if(node->vnode_name=="new_root"){
                LOG(FATAL) << "ERROR: myCheckTreeNodeName";
                exit(EXIT_FAILURE);
            }
            if(node->vnode_name=="new_root"){
                LOG(FATAL) << "ERROR: myCheckTreeNodeName";
                exit(EXIT_FAILURE);
            }
            if(node->vnode_name=="left_node"){
                LOG(FATAL) << "ERROR: myCheckTreeNodeName";
                exit(EXIT_FAILURE);
            }
            if(node->vnode_name=="left_leaf"){
                LOG(FATAL) << "ERROR: myCheckTreeNodeName";
                exit(EXIT_FAILURE);
            }
            if(node->vnode_name=="right_node"){
                LOG(FATAL) << "ERROR: myCheckTreeNodeName";
                exit(EXIT_FAILURE);
            }
            if(node->vnode_name=="right_leaf"){
                LOG(FATAL) << "ERROR: myCheckTreeNodeName";
                exit(EXIT_FAILURE);
            }
        }

    }

    void TreeRearrangment::checkMoveNodeName(Move* move, tshlib::Utree *utree) {

        VirtualNode *node = nullptr;

        node = utree->getNode(move->getSourceNode());

        if(node->vnode_name=="new_root"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="new_root"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="left_node"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="left_leaf"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="right_node"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="right_leaf"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }

        node = utree->getNode(move->getTargetNode());

        if(node->vnode_name=="new_root"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="new_root"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="left_node"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="left_leaf"){
            LOG(FATAL) << "ERROR: myCheckTreeNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="right_node"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
        if(node->vnode_name=="right_leaf"){
            LOG(FATAL) << "ERROR: checkMoveNodeName";
            exit(EXIT_FAILURE);
        }
    }
    //================================================================================

    void TreeRearrangment::moveTheRootLeaf(VirtualNode *pnode, VirtualNode *qnode, std::vector<VirtualNode *> &startVNodes) {

        VirtualNode *new_root  = nullptr;
        VirtualNode *root_2 = nullptr;

        root_2 = qnode->getNodeUp();

        if(pnode->getNodeUp()->vnode_id==root_2->getNodeLeft()->vnode_id){
            new_root = root_2->getNodeRight();
            root_2->_setNodeRight(qnode);
            qnode->_setNodeUp(root_2);
            new_root->_setNodeUp(root_2);
            root_2->_setNodeUp(new_root);
        }else if(pnode->getNodeUp()->vnode_id==root_2->getNodeRight()->vnode_id){
            new_root = root_2->getNodeLeft();
            root_2->_setNodeLeft(qnode);
            qnode->_setNodeUp(root_2);
            new_root->_setNodeUp(root_2);
            root_2->_setNodeUp(new_root);
        }else{
            LOG(FATAL) << "ERROR: moveTheRootLeaf";
            exit(EXIT_FAILURE);
        }

        startVNodes.at(0)=new_root;
        startVNodes.at(1)=root_2;

    }

    void TreeRearrangment::moveTheRoot(VirtualNode *pnode, VirtualNode *qnode, std::vector<VirtualNode *> &startVNodes) {

        VirtualNode *new_root  = nullptr;
        VirtualNode *tmp = nullptr;

        tmp = qnode->getNodeUp();

        new_root = qnode->getNodeLeft();
        new_root->_setNodeUp(qnode);
        qnode->_setNodeLeft(tmp);
        qnode->_setNodeUp(new_root);

        tmp->_setNodeUp(qnode);

        startVNodes.at(0)=new_root;
        startVNodes.at(1)=new_root->getNodeUp();

    }

    void TreeRearrangment::NNImoveUP(tshlib::Move *move,VirtualNode *pnode, VirtualNode *qnode, std::vector<VirtualNode *> &startVNodes) {
        VirtualNode * n1 = nullptr;
        VirtualNode * n2 = nullptr;

        n1 = pnode->getNodeUp();
        n2 = qnode->getNodeUp();

        if(n1->getNodeUp()->vnode_id!=n2->vnode_id){
            VirtualNode * n = n1;
            n1 = n2;
            n2 = n;
            if(n1->getNodeUp()->vnode_id!=n2->vnode_id){
                LOG(FATAL) << "ERROR in NNImoveUP";
                exit(EXIT_FAILURE);
            }
            n = pnode;
            pnode = qnode;
            qnode = n;
        }

        move->node2Opt.push_back(pnode->vnode_id);
        move->node2Opt.push_back(n1->vnode_id);
        move->node2Opt.push_back(n2->vnode_id);
        move->node2Opt.push_back(qnode->vnode_id);

        if(pnode->vnode_id==n1->getNodeLeft()->vnode_id){
            n1->_setNodeLeft(nullptr);
            pnode->_setNodeUp(nullptr);
        }else if(pnode->vnode_id==n1->getNodeRight()->vnode_id){
            n1->_setNodeRight(nullptr);
            pnode->_setNodeUp(nullptr);
        }else{
            LOG(FATAL) << "ERROR in NNImoveUP";
            exit(EXIT_FAILURE);
        }

        if(qnode->vnode_id==n2->getNodeLeft()->vnode_id){
            n2->_setNodeLeft(nullptr);
            qnode->_setNodeUp(nullptr);
        }else if(qnode->vnode_id==n2->getNodeRight()->vnode_id){
            n2->_setNodeRight(nullptr);
            qnode->_setNodeUp(nullptr);
        }else{
            LOG(FATAL) << "ERROR in NNImoveUP";
            exit(EXIT_FAILURE);
        }

        if(!n2->getNodeLeft()){
            n2->_setNodeLeft(pnode);
            pnode->_setNodeUp(n2);
        }else{
            n2->_setNodeRight(pnode);
            pnode->_setNodeUp(n2);
        }

        if(!n1->getNodeLeft()){
            n1->_setNodeLeft(qnode);
            qnode->_setNodeUp(n1);
        }else{
            n1->_setNodeRight(qnode);
            qnode->_setNodeUp(n1);
        }

    }

    void TreeRearrangment::checkMove(Move* move, tshlib::Utree *utree) {

        if(move->moveClassDescription_=="NNI"){

            VirtualNode *pnode = utree->getNode(move->getSourceNode());
            VirtualNode *qnode = utree->getNode(move->getTargetNode());

            if(move->moveDirection_==MoveDirections::up){
                VirtualNode *n1 = pnode->getNodeUp();
                VirtualNode *n2 = qnode->getNodeUp();

                if(n1->getNodeUp()->vnode_id!=n2->vnode_id){
                    VirtualNode * n = n1;
                    n1 = n2;
                    n2 = n;
                    if(n1->getNodeUp()->vnode_id!=n2->vnode_id){
                        LOG(FATAL) << "ERROR in NNImoveUP";
                        exit(EXIT_FAILURE);
                    }
                    move->setSourceNode(qnode->vnode_id);
                    move->setTargetNode(pnode->vnode_id);
                }

            }else if(move->moveDirection_==MoveDirections::up_right || move->moveDirection_==MoveDirections::up_left){

                if(pnode->getNodeUp()->getNodeUp()->getNodeUp()->vnode_id!=qnode->vnode_id){
                    move->setSourceNode(qnode->vnode_id);
                    move->setTargetNode(pnode->vnode_id);
                    pnode = utree->getNode(move->getSourceNode());
                    qnode = utree->getNode(move->getTargetNode());
                    if(pnode->getNodeUp()->getNodeUp()->getNodeUp()->vnode_id!=qnode->vnode_id){


                        std::cout<<utree->printTreeNewick(true,false)<<std::endl;


                        LOG(FATAL) << "ERROR in NNImoveUP_LR";
                        exit(EXIT_FAILURE);
                    }
                }

            }

        }else if(move->moveClassDescription_=="SPR"){

            VirtualNode *source = utree->getNode(move->getSourceNode());
            VirtualNode *target = utree->getNode(move->getTargetNode());

            if(source->isPseudoRootNode() && target->isPseudoRootNode()){
                LOG(FATAL) << "ERROR in checkMove, both nood are pseudoroot";
                exit(EXIT_FAILURE);
            }

            if(source->isPseudoRootNode()){
                move->setSourceNode(target->vnode_id);
                move->setTargetNode(source->vnode_id);
            }

            if(move->moveDirection_==MoveDirections::up){

                VirtualNode *node = source;
                bool flag1 = false;
                while(!node->isPseudoRootNode()){
                    node = node->getNodeUp();
                    if(node->vnode_id==target->vnode_id){
                        flag1=true;
                        break;
                    }
                }
                node = target;
                bool flag2 = false;
                while(!node->isPseudoRootNode()){
                    node = node->getNodeUp();
                    if(node->vnode_id==source->vnode_id){
                        flag2=true;
                        break;
                    }
                }
                if(flag1==true || flag2==true){
                    move->setDirection(MoveDirections::up_left); // either up_left or up_right
                }

            }else if(move->moveDirection_==MoveDirections::up_right || move->moveDirection_==MoveDirections::up_left){


                VirtualNode *node = source;
                bool flag = false;
                while(!node->isPseudoRootNode()){
                    node = node->getNodeUp();
                    if(node->vnode_id==target->vnode_id){
                        flag=true;
                        break;
                    }
                }
                if(!flag){
                    move->setSourceNode(target->vnode_id);
                    move->setTargetNode(source->vnode_id);
                    source = utree->getNode(move->getSourceNode());
                    target = utree->getNode(move->getTargetNode());
                    node = source;
                    flag = false;
                    while(!node->isPseudoRootNode()){
                        node = node->getNodeUp();
                        if(node->vnode_id==target->vnode_id){
                            flag=true;
                            break;
                        }
                    }
                    if(!flag){
                        LOG(FATAL) << "ERROR in checkMove";
                        exit(EXIT_FAILURE);
                    }
                }
            }

        }


    }

    void TreeRearrangment::NNImoveUP_LR(tshlib::Move *move,VirtualNode *pnode, VirtualNode *qnode, std::vector<VirtualNode *> &startVNodes) {
        VirtualNode * n1 = nullptr;
        VirtualNode * n2 = nullptr;
        VirtualNode * sibling = nullptr;

        if(pnode->getNodeUp()->getNodeUp()->getNodeUp()->vnode_id!=qnode->vnode_id){
            VirtualNode * n = pnode;
            pnode = qnode;
            qnode = n;
            if(pnode->getNodeUp()->getNodeUp()->getNodeUp()->vnode_id!=qnode->vnode_id){
                LOG(FATAL) << "ERROR in NNImoveUP_LR";
                exit(EXIT_FAILURE);
            }
        }

        n1 = pnode->getNodeUp();
        n2 = n1->getNodeUp();


        move->node2Opt.push_back(pnode->vnode_id);
        move->node2Opt.push_back(n1->vnode_id);
        move->node2Opt.push_back(n2->vnode_id);
        move->node2Opt.push_back(qnode->vnode_id);


        if(n1->getNodeUp()->vnode_id!=n2->vnode_id){
            VirtualNode * n = n1;
            n1 = n2;
            n2 = n;
            if(n1->getNodeUp()->vnode_id!=n2->vnode_id){
                LOG(FATAL) << "ERROR in NNImoveUP_LR";
                exit(EXIT_FAILURE);
            }
        }

        sibling = pnode->getSiblingNode();

        if(n1->getNodeLeft()->vnode_id==pnode->vnode_id){
            n1->_setNodeLeft(nullptr);
            pnode->_setNodeUp(nullptr);
        }else if(n1->getNodeRight()->vnode_id==pnode->vnode_id){
            n1->_setNodeRight(nullptr);
            pnode->_setNodeUp(nullptr);
        }else{
            LOG(FATAL) << "ERROR in NNImoveUP_LR";
            exit(EXIT_FAILURE);
        }

        if(n2->getNodeLeft()->vnode_id==n1->vnode_id){
            n2->_setNodeLeft(nullptr);
            n1->_setNodeUp(nullptr);
        }else if(n2->getNodeRight()->vnode_id==n1->vnode_id){
            n2->_setNodeRight(nullptr);
            n1->_setNodeUp(nullptr);
        }else{
            LOG(FATAL) << "ERROR in NNImoveUP_LR";
            exit(EXIT_FAILURE);
        }

        if(qnode->getNodeLeft()->vnode_id==n2->vnode_id){
            qnode->_setNodeLeft(nullptr);
            n2->_setNodeUp(nullptr);
        }else if(qnode->getNodeRight()->vnode_id==n2->vnode_id){
            qnode->_setNodeRight(nullptr);
            n2->_setNodeUp(nullptr);
        }else{
            LOG(FATAL) << "ERROR in NNImoveUP_LR";
            exit(EXIT_FAILURE);
        }

        if(!qnode->getNodeLeft()){
            qnode->_setNodeLeft(n1);
            n1->_setNodeUp(qnode);
        }else{
            qnode->_setNodeRight(n1);
            n1->_setNodeUp(qnode);
        }

        if(!n1->getNodeLeft()){
            n1->_setNodeLeft(n2);
            n2->_setNodeUp(n1);
        }else{
            n1->_setNodeRight(n2);
            n2->_setNodeUp(n1);
        }

        if(!n2->getNodeLeft()){
            n2->_setNodeLeft(pnode);
            pnode->_setNodeUp(n2);
        }else{
            n2->_setNodeRight(pnode);
            pnode->_setNodeUp(n2);
        }

    }

    /*
    void TreeRearrangment::NNImoveUP_LR(VirtualNode *pnode, VirtualNode *qnode, std::vector<VirtualNode *> &startVNodes) {
        VirtualNode * n1 = nullptr;
        VirtualNode * n2 = nullptr;
        VirtualNode * tmp1 = nullptr;
        VirtualNode * tmp2 = nullptr;
        int id;
        double bl;
        std::string name;

        if(pnode->getNodeUp()->getNodeUp()->getNodeUp()->vnode_id!=qnode->vnode_id){
            VirtualNode * n = pnode;
            pnode = qnode;
            qnode = n;
            if(pnode->getNodeUp()->getNodeUp()->getNodeUp()->vnode_id!=qnode->vnode_id){
                LOG(FATAL) << "ERROR in NNImoveUP_LR";
                exit(EXIT_FAILURE);
            }
        }

        n1 = pnode->getNodeUp();
        n2 = n1->getNodeUp();

        if(n1->getNodeUp()->vnode_id!=n2->vnode_id){
            VirtualNode * n = n1;
            n1 = n2;
            n2 = n;
            if(n1->getNodeUp()->vnode_id!=n2->vnode_id){
                LOG(FATAL) << "ERROR in NNImoveUP_LR";
                exit(EXIT_FAILURE);
            }
        }

        tmp1 = pnode->getSiblingNode();
        tmp2 = n1->getSiblingNode();

        if(tmp1->vnode_id==n1->getNodeLeft()->vnode_id){
            tmp1->getNodeUp()->_setNodeLeft(nullptr);
            tmp1->_setNodeUp(nullptr);
        }else if(tmp1->vnode_id==n1->getNodeRight()->vnode_id){
            tmp1->getNodeUp()->_setNodeRight(nullptr);
            tmp1->_setNodeUp(nullptr);
        }else{
            LOG(FATAL) << "ERROR in NNImoveUP_LR";
            exit(EXIT_FAILURE);
        }

        if(tmp2->vnode_id==n2->getNodeLeft()->vnode_id){
            tmp2->getNodeUp()->_setNodeLeft(nullptr);
            tmp2->_setNodeUp(nullptr);
        }else if(tmp2->vnode_id==n2->getNodeRight()->vnode_id){
            tmp2->getNodeUp()->_setNodeRight(nullptr);
            tmp2->_setNodeUp(nullptr);
        }else{
            LOG(FATAL) << "ERROR in NNImoveUP_LR";
            exit(EXIT_FAILURE);
        }

        if(!n2->getNodeLeft()){
            n2->_setNodeLeft(tmp1);
            tmp1->_setNodeUp(n2);
        }else{
            n2->_setNodeRight(tmp1);
            tmp1->_setNodeUp(n2);
        }

        if(!n1->getNodeLeft()){
            n1->_setNodeLeft(tmp2);
            tmp2->_setNodeUp(n1);
        }else{
            n1->_setNodeRight(tmp2);
            tmp2->_setNodeUp(n1);
        }

        id = n1->vnode_id;
        bl = n1->vnode_branchlength;
        name = n1->vnode_name;
        n1->vnode_id = n2->vnode_id;
        n1->vnode_branchlength = n2->vnode_branchlength;
        n1->vnode_name = n2->vnode_name;
        n2->vnode_id = id;
        n2->vnode_branchlength = bl;
        n2->vnode_name = name;

    }
    */
    void TreeRearrangment::_applyNNI(tshlib::Move *move,VirtualNode *source, VirtualNode *target, MoveDirections move_direction, std::vector<VirtualNode *> &startVNodes) {

        if (source == target) {
            LOG(FATAL) << "[tshlib::swapNode] The source and the target nodes must be different!";
        } else if (target == nullptr) {
            LOG(FATAL) << "[tshlib::swapNode] The target node is empty";
        } else {
            switch (move_direction) {

                case MoveDirections::up:

                    if(source->isPseudoRootNode() || target->isPseudoRootNode()){

                        if(source->isPseudoRootNode()){

                            if(source->isTerminalNode()) {
                                moveTheRootLeaf(target, source, startVNodes);
                                NNImoveUP(move,target, source, startVNodes);
                            }else{
                                moveTheRoot(target,source,startVNodes);
                                NNImoveUP_LR(move,target, source, startVNodes);
                            }

                        }else{

                            if(target->isTerminalNode()){
                                moveTheRootLeaf(source,target,startVNodes);
                                NNImoveUP(move,source,target,startVNodes);
                            }else{
                                moveTheRoot(source,target,startVNodes);
                                NNImoveUP_LR(move,source, target, startVNodes);
                            }

                        }

                    }else{
                        NNImoveUP(move,source,target,startVNodes);
                    }

                    break;

                case MoveDirections::down_left:
                case MoveDirections::down_right:
                case MoveDirections::up_left:
                case MoveDirections::up_right:

                    NNImoveUP_LR(move,source, target, startVNodes);

                    break;

                default:
                    //std::cout << "I cannot move the node" << std::endl;
                    exit(EXIT_FAILURE);
            }

            //execstatus = true;
        }

        //return execstatus;
    }

    void TreeRearrangment::applyMove(Move *move, Utree &_utree__topology,VirtualNode *source,VirtualNode *target) {

        if(move->moveClassDescription_=="NNI"){
            _applyNNI(move,source, target, move->moveDirection_, _utree__topology.startVNodes);
        }else if(move->moveClassDescription_=="SPR"){
            _applySPR(move,source,target,_utree__topology.startVNodes);
        }else{
            LOG(FATAL) << "[tshlib::TreeRearrangment::applyMove] Something went wrong during the application of the move";
            exit(EXIT_FAILURE);
        }

    }

    /*
    bool TreeRearrangment::revertMove(Move *move, Utree &_utree__topology) {

        bool outcomeExecutionMove = false;

        int pnode = move->getTargetNode();
        int qnode = move->getSourceNode();

        if(move->moveClassDescription_=="NNI"){
            outcomeExecutionMove = _applyNNI(_utree__topology.getNode(pnode),_utree__topology.getNode(qnode),MoveDirections::up,_utree__topology.startVNodes);
        }else if(move->moveClassDescription_=="SPR"){
            outcomeExecutionMove = _revertSPR(move, _utree__topology);
        }else{
            LOG(FATAL) << "[tshlib::TreeRearrangment::revertMove] Something went wrong during the application of the move";
            exit(EXIT_FAILURE);
        }

        return outcomeExecutionMove;
    }
    */

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

    //void TreeRearrangment::commitMove(int moveID, Utree &_utree__topology) {
    void TreeRearrangment::commitMove(Move *move, Utree &_utree__topology,VirtualNode *source,VirtualNode *target){

        // apply move
        //applyMove(moveID, _utree__topology);
        this->applyMove(move,_utree__topology,source,target);

        // reset node rotations
        //================================================================================
        // m@x
//        for (auto &nodeID:_utree__topology.getNodeIdsMap()) {
//            _utree__topology.getNode(nodeID.first)->vnode_rotated = NodeRotation::undef;
//        }
        //================================================================================

    }

    void TreeRearrangment::storeMove(Move *inMove) {

        inMove->moveUID_ = (int) trMoveSet.size();
        trMoveSet.push_back(inMove);

    }

//    void TreeRearrangment::setTree(Utree *inTree) {
//        utree_ = inTree;
//    }

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
//            VLOG(2) << "[test  move]\t" << getMove(idMove)->moveClassDescription_ << "." << std::setfill('0') << std::setw(3) << idMove
//                    << " | " << start_col_line << getMove(idMove)->moveScore_ << end_col_line << "\t"
//                    << " | (" << _utree__topology.getNode(getMove(idMove)->getSourceNode())->vnode_name << "->"
//                    << _utree__topology.getNode(getMove(idMove)->getTargetNode())->vnode_name << ")"
                    //<< "\t[" << getMove(idMove)->moveRadius_ << "] | " << getTree()->printTreeNewick(true);
        } else {
            VLOG(2) << "[test  move]\t" << getMove(idMove)->moveClassDescription_ << "." << std::setfill('0') << std::setw(3) << idMove
                    << " | " << start_col_line << getMove(idMove)->moveScore_ << end_col_line << "\t"
                    << " | (" << _utree__topology.getNode(getMove(idMove)->getSourceNode())->vnode_name
                    << "->" << _utree__topology.getNode(getMove(idMove)->getTargetNode())->vnode_name << ")"
                    << "\t[" << getMove(idMove)->moveRadius_ << "]";
        }

    }

    void TreeRearrangment::removeMoveDuplicates(int num_nodes){

        std::vector<std::vector <bool> > Mat(num_nodes);
        for(int i =0;i<num_nodes;i++){
            Mat[i].resize(num_nodes);
            std::fill(Mat[i].begin(),Mat[i].end(),false);
        }

        int num_elements=0;
        for(auto &move: trMoveSet){
            int source_id = move->getSourceNode();
            int target_id = move->getTargetNode();
            move->is_duplicate = false;
            if(Mat[source_id][target_id]==false){
                Mat[source_id][target_id]=true;
                Mat[target_id][source_id]=true;
                num_elements++;
            }else{
                move->is_duplicate = true;
            }
        }

        std::vector<Move *> moves;
        moves.reserve(num_elements);
        for(auto &move: trMoveSet) {
            if(!move->is_duplicate){
                moves.push_back(move);
            }
        }
        trMoveSet.clear();
        trMoveSet=moves;
    }

    /*
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
    */

    const std::vector<int> TreeRearrangment::updatePathBetweenNodes(unsigned long moveID, std::vector<int> inPath) {

        // Retrieve move object
        Move *move = getMove(moveID);

        // Declaring temporary objects
        std::vector<int> tmpVector_B, tmpVector_C, updatedNodesInPath, outNodePath;
        std::ptrdiff_t pos;

        // According to the different kind of move, apply the specific method to compute the list of nodes
        if(move->moveClassDescription_=="NNI"){



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


        }else if(move->moveClassDescription_=="SPR"){


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


        }else{
            LOG(FATAL) << "[tshlib::TreeRearrangment::updatePathBetweenNodes] Something went wrong during the update";
            exit(EXIT_FAILURE);
        }



        return outNodePath;
    }

    void TreeRearrangment::initialize(int n) {

        // Define the radius for pruning and regrafting the input tree.
        switch (trTreeCoverage) {

            case tshlib::TreeRearrangmentOperations::classic_NNI:
                trSearchRadius_min = 3;
                trSearchRadius_max = 3;
                break;

            case tshlib::TreeRearrangmentOperations::classic_SPR:
                trSearchRadius_min = 4;
                trSearchRadius_max = n / 2;
                break;

            case tshlib::TreeRearrangmentOperations::classic_TBR:
                trSearchRadius_min = 5;
                trSearchRadius_max = n / 2;
                break;

            case tshlib::TreeRearrangmentOperations::classic_Mixed:
                trSearchRadius_min = 3;  // Minimum radius for an NNI move is 3 nodes
                trSearchRadius_max = n; // Full tree traversing from any nodeInterface of the tree
                break;

        }

//        VLOG(1) << "[tsh] Defined coverage radius as [" << trSearchRadius_min << ";" << trSearchRadius_max << "] on a max of ["
//                << utree_->getMaxNodeDistance() << "]";

        trInitialized_ = true;
    }

    void TreeRearrangment::_applySPR(tshlib::Move *move,VirtualNode *source,VirtualNode *target,std::vector<VirtualNode *> &startVNodes) {

        if(source->isPseudoRootNode() && target->isPseudoRootNode()){
            LOG(FATAL) << "ERROR in _applySPR, both nood are pseudoroot";
            exit(EXIT_FAILURE);
        }

        if(source->isPseudoRootNode()){
            tshlib::VirtualNode *tmp = source;
            source = target;
            target = tmp;
        }

        tshlib::VirtualNode *new_root = new VirtualNode();
        new_root->setNodeName("new_root");
        new_root->setVnode_id(-2);

        tshlib::VirtualNode *left_node = new VirtualNode();
        left_node->setNodeName("left_node");
        left_node->setVnode_id(-3);

        tshlib::VirtualNode *left_leaf = new VirtualNode();
        left_leaf->setNodeName("left_leaf");
        left_leaf->setVnode_id(-4);

        tshlib::VirtualNode *right_node = new VirtualNode();
        right_node->setNodeName("right_node");
        right_node->setVnode_id(-5);

        tshlib::VirtualNode *right_leaf = new VirtualNode();
        right_leaf->setNodeName("right_leaf");
        right_leaf->setVnode_id(-6);

        tshlib::VirtualNode *pseudo_root_0 = startVNodes.at(0);
        tshlib::VirtualNode *pseudo_root_1 = startVNodes.at(1);

        new_root->_setNodeUp(nullptr);
        new_root->_setNodeLeft(left_node);
        new_root->_setNodeRight(right_node);

        left_node->_setNodeUp(new_root);
        left_node->_setNodeLeft(pseudo_root_0);
        pseudo_root_0->_setNodeUp(left_node);
        left_node->_setNodeRight(left_leaf);

        left_leaf->_setNodeUp(left_node);
        left_leaf->_setNodeLeft(nullptr);
        left_leaf->_setNodeRight(nullptr);

        right_node->_setNodeUp(new_root);
        right_node->_setNodeLeft(right_leaf);
        right_node->_setNodeRight(pseudo_root_1);
        pseudo_root_1->_setNodeUp(right_node);

        right_leaf->_setNodeUp(right_node);
        right_leaf->_setNodeLeft(nullptr);
        right_leaf->_setNodeRight(nullptr);

        startVNodes.at(0) = nullptr;
        startVNodes.at(1) = nullptr;

        VirtualNode *siblingSource = nullptr;
        VirtualNode *parentSource = nullptr;
        VirtualNode *grandparentSource = nullptr;
        VirtualNode *parentTarget = nullptr;
        double bl = 0.0;


        parentSource = source->getNodeUp();
        if(parentSource->getNodeUp()->vnode_id == source->vnode_id){
            LOG(FATAL) << "ERROR in _applySPR";
            exit(EXIT_FAILURE);
        }

        if(parentSource->getNodeLeft()->vnode_id == source->vnode_id){
            siblingSource = parentSource->getNodeRight();
        }else{
            siblingSource = parentSource->getNodeLeft();
        }

        grandparentSource = parentSource->getNodeUp();
        parentTarget = target->getNodeUp();

        move->node2Opt.push_back(siblingSource->vnode_id);
        move->node2Opt.push_back(target->vnode_id);
        move->node2Opt.push_back(parentSource->vnode_id);

        bl = siblingSource->vnode_branchlength + parentSource->vnode_branchlength;
        //bl = parentSource->vnode_branchlength;

        if(grandparentSource->getNodeLeft()->vnode_id == parentSource->vnode_id){
            grandparentSource->_setNodeLeft(siblingSource);
            siblingSource->_setNodeUp(grandparentSource);
            parentSource->_setNodeUp(nullptr);
            if(parentSource->getNodeLeft()->vnode_id == source->vnode_id){
                parentSource->_setNodeRight(nullptr);
            }else{
                parentSource->_setNodeLeft(nullptr);
            }
        }else{
            grandparentSource->_setNodeRight(siblingSource);
            siblingSource->_setNodeUp(grandparentSource);
            parentSource->_setNodeUp(nullptr);
            if(parentSource->getNodeLeft()->vnode_id == source->vnode_id){
                parentSource->_setNodeRight(nullptr);
            }else{
                parentSource->_setNodeLeft(nullptr);
            }
        }

        if(parentTarget->getNodeLeft()->vnode_id == target->vnode_id){
            parentTarget->_setNodeLeft(parentSource);
            parentSource->_setNodeUp(parentTarget);
            target->_setNodeUp(nullptr);
        }else{
            parentTarget->_setNodeRight(parentSource);
            parentSource->_setNodeUp(parentTarget);
            target->_setNodeUp(nullptr);
        }

        if(parentSource->getNodeRight() == nullptr){
            parentSource->_setNodeRight(target);
            target->_setNodeUp(parentSource);
        }else{
            parentSource->_setNodeLeft(target);
            target->_setNodeUp(parentSource);
        }

        siblingSource->vnode_branchlength = bl;

        bl = target->vnode_branchlength / 2.0;
        //bl = target->vnode_branchlength;

        target->vnode_branchlength = bl;
        parentSource->vnode_branchlength = bl;

        VirtualNode *r1 = new_root->getNodeLeft()->getNodeLeft();
        VirtualNode *r2 = new_root->getNodeRight()->getNodeRight();
        r1->_setNodeUp(r2);
        r2->_setNodeUp(r1);
        startVNodes.at(0)=r1;
        startVNodes.at(1)=r2;
        delete new_root;
        delete left_node;
        delete left_leaf;
        delete right_node;
        delete right_leaf;

    }

    /*
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
    */

}



