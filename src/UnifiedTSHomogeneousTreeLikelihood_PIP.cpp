/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2019 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of Castor
 *
 * Castor is a computer program whose purpose is to infer phylogentic trees
 * under indel-aware and indel-non-aware substitution models for nucleotide,
 * protein, and codon datasets
 *
 * This software is based and extends the following libraries:
 *
 * - the Bio++ libraries
 *   developed by the Bio++ Development Team <http://biopp.univ-montp2.fr>
 *
 * - The Tree Search Heuristic Library (TSH-LIB)
 *   developed by L. Gatti & M. Maiolo <http://bit.ly/tsh-lib>
 *
 * Castor is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Castor is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Castor. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file UnifiedTSHomogeneousTreeLikelihood_PIP.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 04 2018
 * @version 1.0.10
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
 * @see For more information visit: https://bitbucket.org/lorenzogatti89/castor/wiki/Home
 */

#include <glog/logging.h>
#include "UnifiedTSHomogeneousTreeLikelihood_PIP.hpp"


using namespace bpp;

/*
 * Implementation for the interface  likelihood under tree search engines (all mixed models)
 */

UnifiedTSHomogeneousTreeLikelihood_PIP::UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                                                               const SiteContainer &data,
                                                                               TransitionModel *model,
                                                                               DiscreteDistribution *rDist,
                                                                               tshlib::Utree *utree,
                                                                               UtreeBppUtils::treemap *tm,
                                                                               bool optNumericalDerivatives,
                                                                               std::map <std::string, std::string> &params,
                                                                               const std::string &suffix,
                                                                               bool checkRooted,
                                                                               bool verbose,
                                                                               bool usePatterns) :
        RHomogeneousTreeLikelihood_PIP(tree, utree, data, model, rDist, tm, checkRooted, verbose, usePatterns), utree_(utree) {


    setOptimiser(static_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(this), optNumericalDerivatives, params, suffix, true, verbose, 0);

}

UnifiedTSHomogeneousTreeLikelihood_PIP::UnifiedTSHomogeneousTreeLikelihood_PIP(const Tree &tree,
                                                                               TransitionModel *model,
                                                                               DiscreteDistribution *rDist,
                                                                               tshlib::Utree *utree,
                                                                               UtreeBppUtils::treemap *tm,
                                                                               bool optNumericalDerivatives,
                                                                               std::map <std::string, std::string> &params,
                                                                               const std::string &suffix,
                                                                               bool checkRooted,
                                                                               bool verbose,
                                                                               bool usePatterns) :
        RHomogeneousTreeLikelihood_PIP(tree, utree, model, rDist, tm, checkRooted, verbose, usePatterns), utree_(utree) {


    setOptimiser(static_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(this), optNumericalDerivatives, params, suffix, true, verbose, 0);

}


UnifiedTSHomogeneousTreeLikelihood_PIP::~UnifiedTSHomogeneousTreeLikelihood_PIP() {}


void UnifiedTSHomogeneousTreeLikelihood_PIP::init_(bool usePatterns) {}


void UnifiedTSHomogeneousTreeLikelihood_PIP::addTestLikelihoodData(int idxThread) {

    // allocate thread memory
    for (auto &nodeID:tree_->getNodesId()) {

        // Instantiation empty objects
        std::vector < std::vector < std::vector < double >> > sub;
        std::vector < std::vector < std::vector < double >> > empty;
        std::vector<int> descCount;
        std::vector<bool> setA;

        // Allocate memory
        descCount.resize(nbDistinctSites_);
        setA.resize(nbDistinctSites_);

        VectorTools::resize3(sub, nbDistinctSites_, nbClasses_, nbStates_);
        VectorTools::resize3(empty, 1, nbClasses_, nbStates_);

        // Store references
        testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID] = sub;
        testVectorLikelihoodData_[LKDataClass::empty][idxThread][nodeID] = empty;
        tsTemp_descCountData_[idxThread][nodeID] = descCount;
        tsTemp_setAData_[idxThread][nodeID] = setA;
        tsTemp_node_data_origin[idxThread][nodeID] = false;

    }

}

void UnifiedTSHomogeneousTreeLikelihood_PIP::removeTestLikelihoodData(int idxThread) {

    // deallocate thread memory
    for (auto &nodeID:tree_->getNodesId()) {

        size_t dim1 = testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID].size();
        size_t dim2 = testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID][0].size();

        for (size_t i = 0; i < dim1; ++i) {
            for (size_t j = 0; j < dim2; ++j)
                testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID][i][j].resize(0);

            testVectorLikelihoodData_[LKDataClass::sub][idxThread][nodeID][i].resize(0);
        }
        testVectorLikelihoodData_[LKDataClass::empty][idxThread][nodeID][0].resize(0);

        tsTemp_descCountData_[idxThread][nodeID].resize(0);
        tsTemp_setAData_[idxThread][nodeID].resize(0);

    }

}


void UnifiedTSHomogeneousTreeLikelihood_PIP::fireTopologyChange(std::vector<int> nodeList) {

    // Recompute the value of the FV 3D arrays
    //computeSubtreeLikelihood(likelihoodDataTest_, likelihoodEmptyDataTest_);
    computeSubtreeLikelihood(nodeList);
    // Compute the insertion histories set (recompute the desc_count and set A)
    setInsertionHistories(*data_, nodeList, &descCountData_, &setAData_);
}


void UnifiedTSHomogeneousTreeLikelihood_PIP::fireTopologyChange(std::vector<int> nodeList, std::map<int, VVVdouble> *ts_lkdata,
                                                                std::map<int, VVVdouble> *ts_lkemptydata,
                                                                std::map<int, std::vector<int>> *ts_desccount,
                                                                std::map<int, std::vector<bool>> *ts_setadata,
                                                                std::map<int, bool> *ts_node__data_origin,
                                                                tshlib::Utree &_utree__topology) {


    //========================================================
    // m@x
    //????????????????????????????

    //computeAllTransitionProbabilities();
    for(int i=0;i<_utree__topology.listVNodes.size();i++){
        double l = _utree__topology.listVNodes.at(i)->vnode_branchlength;
        int id = _utree__topology.listVNodes.at(i)->vnode_id;
        DiscreteDistribution *rDist_ = this->getRateDistribution();
        // Computes all pxy once for all:
        for (size_t c = 0; c < nbClasses_; c++){
            VVdouble* pxy__c = &pxy_[id][c];
            RowMatrix<double> Q = model_->getPij_t(l * rDist_->getCategory(c));
            for (size_t x = 0; x < nbStates_; x++){
                Vdouble* pxy__c_x = &(*pxy__c)[x];
                for (size_t y = 0; y < nbStates_; y++){
                    (*pxy__c_x)[y] = Q(x, y);
                }
            }
        }
    }
    //========================================================

    // Recompute the value of the FV 3D arrays
    //computeSubtreeLikelihood(likelihoodDataTest_, likelihoodEmptyDataTest_);
    computeSubtreeLikelihood(ts_lkdata, ts_lkemptydata, nodeList, ts_node__data_origin, _utree__topology);
    // Compute the insertion histories set (recompute the desc_count and set A)
    setInsertionHistories(*data_, nodeList, ts_desccount, ts_setadata, ts_node__data_origin, _utree__topology);

    //========================================================
    // m@x
    this->tau_ = _utree__topology.computeTotalTreeLength();
    computeNu();
    setAllIotas(&_utree__topology);
    setAllBetas(&_utree__topology);
    //========================================================

}

double UnifiedTSHomogeneousTreeLikelihood_PIP::updateLikelihoodOnTreeRearrangement(std::vector<int> &nodeList,
                                                                                   tshlib::Utree &_utree__topology,
                                                                                   int idxThread,
                                                                                   std::vector<int> &nodeInvolved) {

    //fetch temporary arrays
    std::map<int, VVVdouble> *ts_lkdata = &testVectorLikelihoodData_[LKDataClass::sub][idxThread];
    std::map<int, VVVdouble> *ts_lkemptydata = &testVectorLikelihoodData_[LKDataClass::empty][idxThread];
    std::map<int, std::vector<int>> *ts_desccount = &tsTemp_descCountData_[idxThread];
    std::map<int, std::vector<bool>> *ts_setadata = &tsTemp_setAData_[idxThread];
    std::map<int, bool> *ts_node__data_origin = &tsTemp_node_data_origin[idxThread];


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(false){

    // Add root to the utree structure
    //**************************************************************
    //_utree__topology.addVirtualRootNode();
    //if (_utree__topology.rootnode->getNodeRight() == nullptr && _utree__topology.rootnode->getNodeLeft() == nullptr) {
        //bpp::ApplicationTools::displayMessage("\nStep2a");
        //----------------------------------------------------------
        //_utree__topology._updateStartNodes();
        try {
            //bpp::ApplicationTools::displayMessage("\nStep2a_1");
            bool fixPseudoRootOnNextSubtree = true;
            //bpp::ApplicationTools::displayMessage("\nStep2a_2");
            //......................................................
            //std::vector<tshlib::VirtualNode *> sideA = _utree__topology.findPseudoRoot(_utree__topology.listVNodes.at(0), fixPseudoRootOnNextSubtree);
            std::vector<tshlib::VirtualNode *> sideA;
            tshlib::VirtualNode *CurrentNode = _utree__topology.listVNodes.at(0);
            //bpp::ApplicationTools::displayResult("Node:",CurrentNode->getVnode_id());
            while (CurrentNode != nullptr) {
                sideA.push_back(CurrentNode);
                if (CurrentNode == CurrentNode->getNodeUp()->getNodeUp()) {
                    //bpp::ApplicationTools::displayMessage("exiting while loop\n");
                    break;
                }
                if (CurrentNode->getNodeUp() != nullptr) {
                    CurrentNode = CurrentNode->getNodeUp();
                    //bpp::ApplicationTools::displayResult("Node:",CurrentNode->getVnode_id());
                } else {
                    CurrentNode = nullptr;
                    //bpp::ApplicationTools::displayMessage("nullptr\n");
                }
            }
            if (fixPseudoRootOnNextSubtree) {

                //bpp::ApplicationTools::displayMessage("fixPseudoRootOnNextSubtree\n");
                if (CurrentNode->getNodeUp() != nullptr) {
                    sideA.push_back(CurrentNode->getNodeUp());
                }
            }
            //......................................................
            //bpp::ApplicationTools::displayMessage("\nStep2a_3");
            tshlib::VirtualNode *sideA_node = sideA.back();
            //bpp::ApplicationTools::displayMessage("\nStep2a_4");
            fixPseudoRootOnNextSubtree = false;
            //bpp::ApplicationTools::displayMessage("\nStep2a_5");
            std::vector<tshlib::VirtualNode *> sideB = _utree__topology.findPseudoRoot(_utree__topology.listVNodes.at(0), fixPseudoRootOnNextSubtree);
            tshlib::VirtualNode *sideB_node = sideB.back();
            //bpp::ApplicationTools::displayMessage("\nStep2a_6");
            _utree__topology.startVNodes.clear();
            //bpp::ApplicationTools::displayMessage("\nStep2a_7");
            _utree__topology.startVNodes.push_back(sideA_node);
            //bpp::ApplicationTools::displayMessage("\nStep2a_8");
            _utree__topology.startVNodes.push_back(sideB_node);
            //bpp::ApplicationTools::displayMessage("\nStep2a_9");
        } catch (const std::exception &e) {
            LOG(FATAL) << "[Utree::_updateStartNodes]" << e.what();
        }
        //-----------------------------------------------------------
        //bpp::ApplicationTools::displayMessage("\nStep2b");
        _utree__topology.rootnode->clearChildren();
        //bpp::ApplicationTools::displayMessage("\nStep2c");
        _utree__topology.startVNodes.at(0)->disconnectNode();
        //bpp::ApplicationTools::displayMessage("\nStep2d");
        _utree__topology.rootnode->connectNode(_utree__topology.startVNodes.at(0));
        _utree__topology.rootnode->connectNode(_utree__topology.startVNodes.at(1));
        //bpp::ApplicationTools::displayMessage("\nStep2e");
//    }else{
//        //bpp::ApplicationTools::displayMessage("\nStep2z");
//    }
    //**************************************************************

    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //===============================================================
    // m@x
    std::vector<int> dist_to_root(_utree__topology.listVNodes.size());
    std::vector<int> indeces(_utree__topology.listVNodes.size());
    for(int i=0;i<_utree__topology.listVNodes.size();i++){
        indeces.at(i)=i;
        //int id = utree_->listVNodes.at(i);
        dist_to_root.at(i)=0;
        tshlib::VirtualNode *node =  _utree__topology.listVNodes.at(i);
        if(node != nullptr){
            while (node->getNodeUp() != nullptr){
                dist_to_root.at(i) += 1;
                node=node->getNodeUp();
            }
        }
    }
    this->sortData(dist_to_root,indeces);
    std::vector<int> _affected__nodes(_utree__topology.listVNodes.size());
    for(int i=0;i<_utree__topology.listVNodes.size();i++){
        _affected__nodes.at(i)=_utree__topology.listVNodes.at(indeces.at(i))->vnode_id;
    }
    // add root node
    //_affected__nodes.push_back(utree_->listVNodes.size());
    //===============================================================

    //================================================================
    // 0. convert the list of tshlib::VirtualNodes into bpp::Node
    //std::vector<int> _affected__nodes = remapVirtualNodeLists(nodeList);
    // m@x
    //std::vector<int> _affected__nodes = remapVirtualNodeLists(full_node_list);
    //================================================================


    // 1. Flag nodes to read from reference
    //for (std::vector<int>::iterator it = _affected__nodes.begin(); it != _affected__nodes.end(); ++it)
    /*
    for (auto it = _affected__nodes.begin(); it != _affected__nodes.end(); ++it){
        (*ts_node__data_origin)[(*it)] = true;
    }
    */

    for (auto it = ts_node__data_origin->begin(); it != ts_node__data_origin->end(); ++it){
        it->second = true;
    }

    //================================================================
//    for(int i=0;i<nodeInvolved.size();i++) {
//
//
//        VVVdouble *pxy__node1 = &pxy_[nodeInvolved.at(i)];
//
//        bpp::Node *node = tree_->getNode(nodeInvolved.at(i), true);
//        computeTransitionProbabilitiesForNode(node);
//
//
//        VVVdouble *pxy__node2 = &pxy_[nodeInvolved.at(i)];
//
//        int kkk=1;
//    }

//    ParameterList parameters;
//    // For each node involved in the move, get the corrisponding branch parameter (no root)
//    for(int i=0;i<nodeInvolved.size();i++) {
//        bpp::Node *bnode = tree_->getNode(nodeInvolved.at(i), true);
//        if (bnode->hasFather()) {
//            Parameter brLen = lk_->getParameter("BrLen" + TextTools::toString(bnode->getId()));
//            brLen.setName("BrLen" + TextTools::toString(bnode->getId()));
//            parameters.addParameter(brLen);
//        }
//    }
//
//    // set parameters on the likelihood function (inherited)
//    lk_->setParameters(parameters);
    //================================================================



    // 2. Fire topology change
    fireTopologyChange(_affected__nodes, ts_lkdata, ts_lkemptydata, ts_desccount, ts_setadata, ts_node__data_origin, _utree__topology);




    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    /*
    map<int, std::vector<bool>>::iterator itr;
    std::vector<unsigned long> idx = this->rootPatternLinksInverse_;
    for(int i=0;i<idx.size();i++){
        std::cout<<idx.at(i)<<",";
    }
    std::cout<<endl;
    int k=0;
    for (itr = (*ts_setadata).begin(); itr != (*ts_setadata).end(); ++itr) {
        int id = itr->first;
        std::vector<bool> bo = itr->second;
        if(id==_utree__topology.listVNodes.size()){
            std::cout<<"root";
        }else{
            std::cout<<_utree__topology.getNode(id)->vnode_name;
        }
        std::cout<<" : ";
        for(int j=0;j<bo.size();j++){
            std::cout<<bo.at(j);
            std::cout<<",";
        }
        std::cout<<endl;
        k++;
    }
     */
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










    // 3. Compute loglikelihood
    double logLk = getLogLikelihoodOnTreeRearrangement(_affected__nodes,
                                                       ts_lkdata,
                                                       ts_lkemptydata,
                                                       ts_desccount,
                                                       ts_setadata,
                                                       ts_node__data_origin,
                                                       _utree__topology);

    // 4. Remove root node from the utree structure
    //_utree__topology.removeVirtualRootNode();

    // 5. Remove flags
    for (auto it = _affected__nodes.begin(); it != _affected__nodes.end(); ++it) {
        (*ts_node__data_origin)[(*it)] = false;
        std::fill((*ts_setadata)[(*it)].begin(), (*ts_setadata)[(*it)].end(), false);
        std::fill((*ts_desccount)[(*it)].begin(), (*ts_desccount)[(*it)].end(), 0);
    }

    return logLk;


//#endif FORCE_FULL_LK_COMPUTATION

}

void UnifiedTSHomogeneousTreeLikelihood_PIP::sortData(std::vector<int> &array,std::vector<int> &indeces) const{

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

double UnifiedTSHomogeneousTreeLikelihood_PIP::getLogLikelihoodOnTreeRearrangement(const std::vector<int> &_ts__node_list,
                                                                                   std::map<int, VVVdouble> *_ts__likelihoods,
                                                                                   std::map<int, VVVdouble> *_ts__likelihoods_empty,
                                                                                   std::map<int, std::vector<int>> *_ts__desccount,
                                                                                   std::map<int, std::vector<bool>> *_ts__setadata,
                                                                                   std::map<int, bool> *_ts__node_data_origin,
                                                                                   tshlib::Utree &_utree__topology) const {

    // 1. Initialise variables and contenitors
    double logLK;
    std::vector<double> lk_sites(nbDistinctSites_);

    // 2. Compute the lk of the empty column
    double lk_site_empty = computeLikelihoodWholeAlignmentEmptyColumn(_ts__likelihoods_empty, _ts__node_data_origin, _utree__topology);



    //!!!!!!!!!!!!!!!!!!!!
    //std::cout<<std::setprecision(18)<<"lk[0]:"<<lk_site_empty<<std::endl;
    //!!!!!!!!!!!!!!!!!!!!

    // 3. Compute the likelihood of each site
    const std::vector<unsigned int> *_root__weights = &likelihoodData_->getWeights();

    //===============================================================
    // m@x
    std::vector<int> dist_to_root(_utree__topology.listVNodes.size());
    std::vector<int> indeces(_utree__topology.listVNodes.size());
    for(int i=0;i<_utree__topology.listVNodes.size();i++){
        indeces.at(i)=i;
        dist_to_root.at(i)=0;
        tshlib::VirtualNode *node =  _utree__topology.listVNodes.at(i);
        if(node != nullptr){
            while (node->getNodeUp() != nullptr){
                dist_to_root.at(i) += 1;
                node=node->getNodeUp();
            }
        }
    }
    this->sortData(dist_to_root,indeces);
    std::vector<int> full_node_list(_utree__topology.listVNodes.size());
    for(int i=0;i<_utree__topology.listVNodes.size();i++){
        full_node_list.at(i)=_utree__topology.listVNodes.at(indeces.at(i))->vnode_id;
    }
    full_node_list.push_back(_utree__topology.listVNodes.size());
    //===============================================================

    for (unsigned long i = 0; i < nbDistinctSites_; i++) {

        // Extend rearranged-node-list including all the nodes in the setA for each site
        //std::vector<int> _node__list;
        //_extendNodeListOnSetA(_ts__node_list.back(), i, _node__list, _ts__setadata, _utree__topology);


        // call to function which retrieves the lk value for each site
        //==============================================================================================================
        // m@x
        /*
        lk_sites[i] = log(computeLikelihoodForASite(i,
                                                    _ts__likelihoods,
                                                    _ts__likelihoods_empty,
                                                    _ts__setadata,
                                                    _node__list,
                                                    _ts__node_data_origin,
                                                    _utree__topology)) * _root__weights->at(i);
        */
        /*
        lk_sites[i] = log(computeLikelihoodForASite(i,
                                                    _ts__likelihoods,
                                                    _ts__likelihoods_empty,
                                                    _ts__setadata,
                                                    const_cast<std::vector<int> &>(_ts__node_list),
                                                    _ts__node_data_origin,
                                                    _utree__topology)) * _root__weights->at(i);
        */
        lk_sites[i] = log(computeLikelihoodForASite(i,
                                            _ts__likelihoods,
                                            _ts__likelihoods_empty,
                                            _ts__setadata,
                                            const_cast<std::vector<int> &>(full_node_list),
                                            _ts__node_data_origin,
                                            _utree__topology)) * _root__weights->at(i);
        //==============================================================================================================
        //=========================
        // m@x
        //!!!!!!!!!!!!!!!!!!!!
        //std::cout<<std::setprecision(18)<<"lk["<<i+1<<"]:"<<lk_sites[i]<<std::endl;
        //!!!!!!!!!!!!!!!!!!!!
        //=========================

        DVLOG(2) << "site log_lk[" << i << "]=" << std::setprecision(18) << lk_sites[i] << std::endl;
    }

    // Sum all the values stored in the lk vector
    logLK = MatrixBppUtils::sumVector(&lk_sites);
    DVLOG(2) << "LK Sites [BPP] " << std::setprecision(18) << logLK;

    // compute PHi
    double log_phi_value = computePhi(lk_site_empty);
    DVLOG(2) << "PHI [BPP] " << std::setprecision(18) << log_phi_value;

    logLK += log_phi_value;


    return logLK;
}


void UnifiedTSHomogeneousTreeLikelihood_PIP::topologyChangeSuccessful(std::vector<int> listNodes) {

    // Update BPP tree using the structure in Utree
    topologyCommitTree();

    // Add virtual root to compute the likelihood
    utree_->addVirtualRootNode();

    // Fire topology change
    std::vector<int> ponl = getNodeListPostOrder(tree_->getRootNode()->getId());

    setLikelihoodNodes(ponl);

    fireTopologyChange(ponl);

    DVLOG(1) << "loglikelihood after commit" << TextTools::toString(getLogLikelihood(), 15);

    //========================================================================================================
    //m@x
    // Optimise branches involved in the tree rearrangement
    //fireBranchOptimisation(UtreeBppUtils::remapNodeLists(listNodes, tree_, treemap_));
    //========================================================================================================

    // Remove the virtual root to allow for further tree topology improvements
    utree_->removeVirtualRootNode();

}

void UnifiedTSHomogeneousTreeLikelihood_PIP::commitBranchLength(int id_Bpp,double bl) {

    if(bl<MIN_BRANCH_LEN){
        bl=MIN_BRANCH_LEN;
    }

    this->nodes_.at(id_Bpp)->setDistanceToFather(bl);

}

void UnifiedTSHomogeneousTreeLikelihood_PIP::commitBranchLength(tshlib::Utree *utree) {

    for(int i=0;i<utree->listVNodes.size();i++){
        int id_Bpp = this->treemap_.right.at(utree->listVNodes.at(i)->vnode_id);
        double bl = this->nodes_.at(id_Bpp)->getDistanceToFather();

        if(bl<MIN_BRANCH_LEN){
            bl=MIN_BRANCH_LEN;
        }

        utree->listVNodes.at(i)->vnode_branchlength = bl;
    }

}

void UnifiedTSHomogeneousTreeLikelihood_PIP::commitBranchLengths() {

    for(int i=0;i<this->nodes_.size();i++){
        int id = this->nodes_.at(i)->getId();
        double bl = this->nodes_.at(i)->getDistanceToFather();

        if(bl<MIN_BRANCH_LEN){
            bl=MIN_BRANCH_LEN;
        }

        tree_->getNode(id)->setDistanceToFather(bl);
    }

}

void UnifiedTSHomogeneousTreeLikelihood_PIP::topologyCommitTree() {

    std::vector < tshlib::VirtualNode * > nodelist;
    nodelist = utree_->listVNodes;

    std::map < int, bpp::Node * > tempMap;
    std::map<int, double> tempDistanceToFather;

    double bl;

    // reset inBtree
    for (auto &bnode:tree_->getNodes()) {

        tempMap.insert(std::pair<int, bpp::Node *>(bnode->getId(), bnode));
        // Empty array of sons on the node
        bnode->removeSons();

        if (bnode->hasFather()) {

            tempDistanceToFather.insert(std::pair<int, double>(bnode->getId(), bnode->getDistanceToFather()));
            // Empty father connection
            bnode->removeFather();
        }
    }

    for (auto &vnode:nodelist) {

        if (!vnode->isTerminalNode()) {

            // get corresponding sons in inBTree
            bpp::Node *leftBNode = tempMap[treemap_.right.at(vnode->getNodeLeft()->getVnode_id())];
            bpp::Node *rightBNode = tempMap[treemap_.right.at(vnode->getNodeRight()->getVnode_id())];

            // get corresponding parent in inBTree
            bpp::Node *pNode = tempMap[treemap_.right.at(vnode->getVnode_id())];

            leftBNode->setFather(pNode);
            rightBNode->setFather(pNode);

            bl=tempDistanceToFather[leftBNode->getId()];
            if(bl<MIN_BRANCH_LEN){
                bl=MIN_BRANCH_LEN;
            }
            leftBNode->setDistanceToFather(bl);

            bl=tempDistanceToFather[rightBNode->getId()];
            if(bl<MIN_BRANCH_LEN){
                bl=MIN_BRANCH_LEN;
            }
            rightBNode->setDistanceToFather(bl);

            //Add new sons
            pNode->setSon(0, leftBNode);
            pNode->setSon(1, rightBNode);

            pNode->setDistanceToFather(tempDistanceToFather[pNode->getId()]);

        }

        // in case the current vnode is also the pseudo-root
        if (vnode == vnode->getNodeUp()->getNodeUp()) {

            bpp::Node *leftBNode = tempMap[treemap_.right.at(vnode->getVnode_id())];
            bpp::Node *rightBNode = tempMap[treemap_.right.at(vnode->getNodeUp()->getVnode_id())];

            tree_->getRootNode()->removeSons();

            leftBNode->setFather(tree_->getRootNode());
            rightBNode->setFather(tree_->getRootNode());

            bl=tempDistanceToFather[leftBNode->getId()];
            if(bl<MIN_BRANCH_LEN){
                bl=MIN_BRANCH_LEN;
            }
            leftBNode->setDistanceToFather(bl);

            bl=tempDistanceToFather[rightBNode->getId()];
            if(bl<MIN_BRANCH_LEN){
                bl=MIN_BRANCH_LEN;
            }
            rightBNode->setDistanceToFather(bl);

            tree_->getRootNode()->setSon(0, leftBNode);
            tree_->getRootNode()->setSon(1, rightBNode);

        }

    }

}
