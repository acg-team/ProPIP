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

    // Recompute the value of the FV 3D arrays
    //computeSubtreeLikelihood(likelihoodDataTest_, likelihoodEmptyDataTest_);
    computeSubtreeLikelihood(ts_lkdata, ts_lkemptydata, nodeList, ts_node__data_origin, _utree__topology);
    // Compute the insertion histories set (recompute the desc_count and set A)
    setInsertionHistories(*data_, nodeList, ts_desccount, ts_setadata, ts_node__data_origin, _utree__topology);
}


double UnifiedTSHomogeneousTreeLikelihood_PIP::updateLikelihoodOnTreeRearrangement(std::vector<int> &nodeList,
                                                                                   tshlib::Utree &_utree__topology,
                                                                                   int idxThread) {
    //fetch temporary arrays
    std::map<int, VVVdouble> *ts_lkdata = &testVectorLikelihoodData_[LKDataClass::sub][idxThread];
    std::map<int, VVVdouble> *ts_lkemptydata = &testVectorLikelihoodData_[LKDataClass::empty][idxThread];
    std::map<int, std::vector<int>> *ts_desccount = &tsTemp_descCountData_[idxThread];
    std::map<int, std::vector<bool>> *ts_setadata = &tsTemp_setAData_[idxThread];
    std::map<int, bool> *ts_node__data_origin = &tsTemp_node_data_origin[idxThread];


    // Add root to the utree structure
    _utree__topology.addVirtualRootNode();

    // 0. convert the list of tshlib::VirtualNodes into bpp::Node
    std::vector<int> _affected__nodes = remapVirtualNodeLists(nodeList);

    // 1. Flag nodes to read from reference
    //for (std::vector<int>::iterator it = _affected__nodes.begin(); it != _affected__nodes.end(); ++it)
    for (auto it = _affected__nodes.begin(); it != _affected__nodes.end(); ++it)
        (*ts_node__data_origin)[(*it)] = true;

    // 2. Fire topology change
    fireTopologyChange(_affected__nodes, ts_lkdata, ts_lkemptydata, ts_desccount, ts_setadata, ts_node__data_origin, _utree__topology);

    // 3. Compute loglikelihood
    double logLk = getLogLikelihoodOnTreeRearrangement(_affected__nodes,
                                                       ts_lkdata,
                                                       ts_lkemptydata,
                                                       ts_desccount,
                                                       ts_setadata,
                                                       ts_node__data_origin,
                                                       _utree__topology);

    // 4. Remove root node from the utree structure
    _utree__topology.removeVirtualRootNode();



    // 5. Remove flags
    for (auto it = _affected__nodes.begin(); it != _affected__nodes.end(); ++it) {
        (*ts_node__data_origin)[(*it)] = false;
        std::fill((*ts_setadata)[(*it)].begin(), (*ts_setadata)[(*it)].end(), false);
        std::fill((*ts_desccount)[(*it)].begin(), (*ts_desccount)[(*it)].end(), 0);
    }

    return logLk;

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

    // 3. Compute the likelihood of each site
    const std::vector<unsigned int> *_root__weights = &likelihoodData_->getWeights();

    for (unsigned long i = 0; i < nbDistinctSites_; i++) {

        // Extend rearranged-node-list including all the nodes in the setA for each site
        std::vector<int> _node__list;
        _extendNodeListOnSetA(_ts__node_list.back(), i, _node__list, _ts__setadata, _utree__topology);

        // call to function which retrieves the lk value for each site
        lk_sites[i] = log(computeLikelihoodForASite(i,
                                                    _ts__likelihoods,
                                                    _ts__likelihoods_empty,
                                                    _ts__setadata,
                                                    _node__list,
                                                    _ts__node_data_origin,
                                                    _utree__topology)) * _root__weights->at(i);

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

    // Optimise branches involved in the tree rearrangement
    fireBranchOptimisation(UtreeBppUtils::remapNodeLists(listNodes, tree_, treemap_));

    // Remove the virtual root to allow for further tree topology improvements
    utree_->removeVirtualRootNode();

}


void UnifiedTSHomogeneousTreeLikelihood_PIP::topologyCommitTree() {

    std::vector < tshlib::VirtualNode * > nodelist;
    nodelist = utree_->listVNodes;

    std::map < int, bpp::Node * > tempMap;
    std::map<int, double> tempDistanceToFather;
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

            leftBNode->setDistanceToFather(tempDistanceToFather[leftBNode->getId()]);
            rightBNode->setDistanceToFather(tempDistanceToFather[rightBNode->getId()]);

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

            leftBNode->setDistanceToFather(tempDistanceToFather[leftBNode->getId()]);
            rightBNode->setDistanceToFather(tempDistanceToFather[rightBNode->getId()]);

            tree_->getRootNode()->setSon(0, leftBNode);
            tree_->getRootNode()->setSon(1, rightBNode);

        }

    }

}
