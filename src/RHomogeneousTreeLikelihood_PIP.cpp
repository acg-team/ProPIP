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
 * @file RHomogeneousTreeLikelihood_PIP.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 11 01 2018
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

#include <vector>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <glog/logging.h>
#include "Utils.hpp"

using namespace bpp;

#include "RHomogeneousTreeLikelihood_PIP.hpp"

RHomogeneousTreeLikelihood_PIP::RHomogeneousTreeLikelihood_PIP(
        const Tree &tree,
        tshlib::Utree *utree,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        UtreeBppUtils::treemap *tm,
        bool checkRooted,
        bool verbose,
        bool usePatterns)
throw(Exception) :
        AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
        likelihoodData_(0),
        minusLogLik_(-1.),
        utree_(utree),
        treemap_(*tm) {


    init_(usePatterns);

}

RHomogeneousTreeLikelihood_PIP::RHomogeneousTreeLikelihood_PIP(
        const Tree &tree,
        tshlib::Utree *utree,
        const SiteContainer &data,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        UtreeBppUtils::treemap *tm,
        bool checkRooted,
        bool verbose,
        bool usePatterns)
throw(Exception) :
        AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
        likelihoodData_(0),
        minusLogLik_(-1.),
        utree_(utree),
        treemap_(*tm) {

    init_(usePatterns);
    setData(data);


}

RHomogeneousTreeLikelihood_PIP::RHomogeneousTreeLikelihood_PIP(
        const RHomogeneousTreeLikelihood_PIP &lik) :
        AbstractHomogeneousTreeLikelihood(lik),
        likelihoodData_(0),
        minusLogLik_(lik.minusLogLik_) {
    likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData *>(lik.likelihoodData_->clone());
    likelihoodData_->setTree(tree_);
}


RHomogeneousTreeLikelihood_PIP &RHomogeneousTreeLikelihood_PIP::operator=(const RHomogeneousTreeLikelihood_PIP &lik) {
    AbstractHomogeneousTreeLikelihood::operator=(lik);
    if (likelihoodData_) delete likelihoodData_;
    likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData *>(lik.likelihoodData_->clone());
    likelihoodData_->setTree(tree_);
    minusLogLik_ = lik.minusLogLik_;
    return *this;
}


RHomogeneousTreeLikelihood_PIP::~RHomogeneousTreeLikelihood_PIP() {
    delete likelihoodData_;
}


void RHomogeneousTreeLikelihood_PIP::init_(bool usePatterns) noexcept(false) {
    // This call initialises the data structure to compute the partial likelihoods of the nodes (it allows for ASRV distributions to be added on
    // top of the s.m.)
    likelihoodData_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);
    likelihoodEmptyData_ = new DRASRTreeLikelihoodData(tree_, rateDistribution_->getNumberOfCategories(), usePatterns);
}


void RHomogeneousTreeLikelihood_PIP::setData(const SiteContainer &sites) throw(Exception) {

    if (data_)
        delete data_;

    data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());

    if (verbose_)
        ApplicationTools::displayTask("Initializing data structure");

    likelihoodData_->initLikelihoods(*data_, *model_);

    // Initialise gap column
    auto gaps = new std::vector<int>(sites.getNumberOfSequences(), sites.getAlphabet()->getGapCharacterCode());
    auto emptySite = new Site(*gaps, sites.getAlphabet());

    auto emptyContainer = new VectorSiteContainer(sites.getSequencesNames(), sites.getAlphabet());
    emptyContainer->addSite(*emptySite, false);

    // Prepare FV empty
    likelihoodEmptyData_->initLikelihoods(*emptyContainer, *model_);

    delete gaps;
    delete emptySite;
    delete emptyContainer;


    if (verbose_)
        ApplicationTools::displayTaskDone();

    nbSites_ = likelihoodData_->getNumberOfSites();
    nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
    nbStates_ = likelihoodData_->getNumberOfStates();

    //auto rootPatternLinksInverse_ = new std::vector<unsigned long>(nbDistinctSites_);
    rootPatternLinksInverse_.resize(nbDistinctSites_);

    // fill inverse map patternlinks
    auto rootPatternLinks = likelihoodData_->getRootArrayPositions();
    for (int i = 0; i < nbSites_; i++) {
        rootPatternLinksInverse_.at(rootPatternLinks.at(i)) = static_cast<unsigned long>(i);
    }


    // Initialise iotas, betas, evolutionary events and weighted evolutionary events maps
    for (auto &nodeID:tree_->getNodesId()) {
        iotasData_.insert(std::make_pair(nodeID, 0));
        betasData_.insert(std::make_pair(nodeID, 0));
        if (tree_->isLeaf(nodeID))
            indicatorFun_[nodeID].resize(nbDistinctSites_);

    }

    // Initialise vectors for storing insertion histories values
    initialiseInsertionHistories();

    // Set indicator function on leafs
    setIndicatorFunction(*data_);

    // Set the likelihoodNodeList to the postorder one
    std::vector<int> ponl = getNodeListPostOrder(tree_->getRootNode()->getId());
    setLikelihoodNodes(ponl);

    // Set initial Insertion Histories
    setInsertionHistories(*data_, ponl, &descCountData_, &setAData_);

    if (verbose_)
        ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));

    initialized_ = false;

}


void RHomogeneousTreeLikelihood_PIP::setParameters(const ParameterList &parameters)
throw(ParameterNotFoundException, ConstraintException) {
    setParametersValues(parameters);
}


double RHomogeneousTreeLikelihood_PIP::getValue() const
throw(Exception) {
    if (!isInitialized()) throw Exception("RHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
    return minusLogLik_;
}

void RHomogeneousTreeLikelihood_PIP::_computeHadamardFVSons(std::vector<VVVdouble *> inFVSons, VVVdouble *outFVParent) const {

    size_t t_distSites = (*outFVParent).size();
    size_t t_classes = (*outFVParent)[0].size();
    size_t t_states = (*outFVParent)[0][0].size();

    for (int i = 0; i < t_distSites; i++) {
        for (int c = 0; c < t_classes; c++) {
            for (int x = 0; x < t_states; x++) {
                (*outFVParent)[i][c][x] = 1.0;

                for (int l = 0; l < inFVSons.size(); l++)
                    (*outFVParent)[i][c][x] *= (*inFVSons[l])[i][c][x];

            }
        }
    }
}

Vdouble RHomogeneousTreeLikelihood_PIP::_SingleRateCategoryHadamardMultFvSons(int nodeID,
                                                                              unsigned long site,
                                                                              unsigned long rate,
                                                                              std::vector<VVVdouble *> lk_sons) const {

    std::vector<double> fv_out(nbStates_, 0);

    for (int x = 0; x < nbStates_; x++) {
        fv_out[x] = 1.0;

        for (int i = 0; i < 2; i++)
            fv_out[x] *= (*lk_sons[i])[site][rate][x];

    }

    return fv_out;
}

void RHomogeneousTreeLikelihood_PIP::_computePrTimesFv(VVVdouble *pxy__node, VVVdouble *_likelihoods_node) const {
    size_t t_distSites = (*_likelihoods_node).size();
    size_t t_classes = (*_likelihoods_node)[0].size();
    size_t t_states = (*_likelihoods_node)[0][0].size();

    for (int c = 0; c < t_classes; c++) {
        //VVdouble *pxy__node_c = &(*pxy__node)[c];
        for (int i = 0; i < t_distSites; i++) {
            // this vector stores the old values before recomputing them
            Vdouble reference = (*_likelihoods_node)[i][c];
            for (int x = 0; x < t_states; x++) {
                (*_likelihoods_node)[i][c][x] = 0.0;
                for (int y = 0; y < t_states; y++)
                    (*_likelihoods_node)[i][c][x] += (*pxy__node)[c][x][y] * reference[y];
            }
        }
    }

}

void RHomogeneousTreeLikelihood_PIP::_computePrTimesIndicator(VVVdouble *pxy__node, VVdouble *indicator_node, VVVdouble *_likelihoods_node) const {
    size_t t_distSites = (*_likelihoods_node).size();
    size_t t_classes = (*_likelihoods_node)[0].size();
    size_t t_states = (*_likelihoods_node)[0][0].size();

    for (int c = 0; c < t_classes; c++) {
        for (int i = 0; i < t_distSites; i++) {
            for (int x = 0; x < t_states; x++) {
                // Reset value for the current state
                (*_likelihoods_node)[i][c][x] = 0.0;

                for (int y = 0; y < t_states; y++)
                    (*_likelihoods_node)[i][c][x] += (*pxy__node)[c][x][y] * (*indicator_node)[i][y];

            }
        }
    }
}

void RHomogeneousTreeLikelihood_PIP::_computePrTimesIndicatorEmpty(VVVdouble *pxy__node, VVdouble *indicator_node, VVVdouble *_likelihoods_node)
const {

    for (int c = 0; c < nbClasses_; c++) {
        for (int x = 0; x < nbStates_; x++) {
            (*_likelihoods_node)[0][c][x] = (*pxy__node)[c][x][nbStates_ - 1];
        }
    }
}


std::vector<int> RHomogeneousTreeLikelihood_PIP::_getMappedNodeChildren(int nodeID, tshlib::Utree &_utree__topology) const {

    std::vector<int> sonsIDs(2);
    // TODO: Check if this condition is needed during the remapping
//    if (_utree__topology.getNode(nodeID)->isRootNode()){
//    //if (tree_->isRoot(nodeID)) {
//
//        sonsIDs[0] = treemap_.right.at(_utree__topology.rootnode->getNodeLeft()->getVnode_id());
//        sonsIDs[1] = treemap_.right.at(_utree__topology.rootnode->getNodeRight()->getVnode_id());
//
//    } else {
        sonsIDs[0] = treemap_.right.at(_utree__topology.getNode(treemap_.left.at(nodeID))->getNodeLeft()->getVnode_id());
        sonsIDs[1] = treemap_.right.at(_utree__topology.getNode(treemap_.left.at(nodeID))->getNodeRight()->getVnode_id());
//    }

    //sonsIDs[0] = treemap_.right.at(candUTree.getNodeIdsMap().at(treemap_.left.at(nodeID))->getNodeLeft()->getVnode_id());
    //sonsIDs[1] = treemap_.right.at(candUTree.getNodeIdsMap().at(treemap_.left.at(nodeID))->getNodeRight()->getVnode_id());

    return sonsIDs;
}

void RHomogeneousTreeLikelihood_PIP::computeSubtreeLikelihood() const {

}

void RHomogeneousTreeLikelihood_PIP::computeSubtreeLikelihood(const std::vector<int> &nodeList) {

    for (auto &nodeID:nodeList) {

        // Get references to data-structures for node children and append them in a vector
        auto _sons__ids = new std::vector<int>;// _sons__ids;
        std::vector<VVVdouble *> lk_sons(2);
        std::vector<VVVdouble *> lk_sons_empty(2);

        //if (!tree_->isLeaf(nodeID)) {
        if (!utree_->getNode(treemap_.left.at(nodeID))->isTerminalNode()) {
            (*_sons__ids) = utree_->getChildrenId(treemap_.left.at(nodeID));
            //(*_sons__ids) = tree_->getSonsId(nodeID);
            for (int l = 0; l < 2; l++) {

                lk_sons[l] = &likelihoodData_->getLikelihoodArray((*_sons__ids)[l]);
                lk_sons_empty[l] = &likelihoodEmptyData_->getLikelihoodArray((*_sons__ids)[l]);

            }
        }
        delete _sons__ids;


        // References
        VVVdouble *_likelihoods__node = &likelihoodData_->getLikelihoodArray(nodeID);
        VVVdouble *_likelihoods_empty__node = &likelihoodEmptyData_->getLikelihoodArray(nodeID);
        VVVdouble *pxy__node = &pxy_[nodeID];



        // Call kernel function
        _kernel_subtreelikelihood(nodeID, pxy__node, _likelihoods__node, _likelihoods_empty__node, &lk_sons, &lk_sons_empty);

    }

}

void RHomogeneousTreeLikelihood_PIP::computeSubtreeLikelihood(std::map<int, VVVdouble> *likelihoods,
                                                              std::map<int, VVVdouble> *likelihoods_empty,
                                                              const std::vector<int> &nodeList,
                                                              std::map<int, bool> *ts_node__data_origin,
                                                              tshlib::Utree &_utree__topology) {

    for (auto &nodeID:nodeList) {

        // Get references to data-structures for node children and append them in a vector
        auto _sons__ids = new std::vector<int>;// _sons__ids;
        std::vector<VVVdouble *> lk_sons(2);
        std::vector<VVVdouble *> lk_sons_empty(2);

        //if (!tree_->isLeaf(nodeID)) {
        if (!utree_->getNode(treemap_.left.at(nodeID))->isTerminalNode()) {
            (*_sons__ids) = _getMappedNodeChildren(nodeID, _utree__topology);

            for (int l = 0; l < 2; l++) {
                if ((*ts_node__data_origin)[(*_sons__ids)[l]]) {
                    lk_sons[l] = &(*likelihoods)[(*_sons__ids)[l]];
                    lk_sons_empty[l] = &(*likelihoods_empty)[(*_sons__ids)[l]];
                } else {
                    lk_sons[l] = &likelihoodData_->getLikelihoodArray((*_sons__ids)[l]);
                    lk_sons_empty[l] = &likelihoodEmptyData_->getLikelihoodArray((*_sons__ids)[l]);
                }
            }
        }
        delete _sons__ids;
        // References
        VVVdouble *_likelihoods__node = &(*likelihoods)[nodeID];
        VVVdouble *_likelihoods_empty__node = &(*likelihoods_empty)[nodeID];
        VVVdouble *pxy__node = &pxy_[nodeID];

        // Call kernel function
        _kernel_subtreelikelihood(nodeID, pxy__node, _likelihoods__node, _likelihoods_empty__node, &lk_sons, &lk_sons_empty);

    }

}


void RHomogeneousTreeLikelihood_PIP::_kernel_subtreelikelihood(int nodeID,
                                                               VVVdouble *pxy__node,
                                                               VVVdouble *_likelihoods__node,
                                                               VVVdouble *_likelihoods_empty__node,
                                                               std::vector<VVVdouble *> *lk_sons,
                                                               std::vector<VVVdouble *> *lk_sons_empty) {

    // Leaf node
    //if (tree_->isLeaf(nodeID)) {
    if (utree_->getNode(treemap_.left.at(nodeID))->isTerminalNode()) {
        VVdouble *indicator_node = &indicatorFun_[nodeID];
        _computePrTimesIndicator(pxy__node, indicator_node, _likelihoods__node);
        _computePrTimesIndicatorEmpty(pxy__node, indicator_node, _likelihoods_empty__node);
        return;
    }

    // Update FV children
    _computeHadamardFVSons((*lk_sons), _likelihoods__node);
    _computeHadamardFVSons((*lk_sons_empty), _likelihoods_empty__node);


    // Add transition probability for internal node (not root)

    //if (!tree_->isRoot(nodeID)) {
    if(!utree_->getNode(treemap_.left.at(nodeID))->isRootNode()){
        _computePrTimesFv(pxy__node, _likelihoods__node);
        _computePrTimesFv(pxy__node, _likelihoods_empty__node);
    }

    //--------------------------------------------------------------------------------------------
    // Debug ** array fetching + printing **
    //VVVdouble *_likelihoods_empty_node = &likelihoodEmptyData_->getLikelihoodArray(node->getId());
    //VVVdouble *_likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
    //_printFV(node, _likelihoods_node);
    //_printFV(node, _likelihoods_empty_node);
    //--------------------------------------------------------------------------------------------

}


void RHomogeneousTreeLikelihood_PIP::displayLikelihood(const Node *node) {
    DVLOG(2) << "Likelihoods at node " << node->getName() << ": ";
    displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId()));
    DVLOG(2) << "                                         ***";
}


void RHomogeneousTreeLikelihood_PIP::fireParameterChanged(const ParameterList &params) {

    applyParameters();

    if (rateDistribution_->getParameters().getCommonParametersWith(params).size() > 0
        || model_->getParameters().getCommonParametersWith(params).size() > 0) {
        //Rate parameter changed, need to recompute all probs:
        computeAllTransitionProbabilities();
    } else if (params.size() > 0) {
        //We may save some computations:
        for (size_t i = 0; i < params.size(); i++) {
            std::string s = params[i].getName();
            if (s.substr(0, 5) == "BrLen") {
                //Branch length parameter:
                computeTransitionProbabilitiesForNode(nodes_[TextTools::to<size_t>(s.substr(5))]);
            }
        }
        rootFreqs_ = model_->getFrequencies();
    }

    // Recompute the tree length
    tau_ = tree_->getTotalLength();                     //TODO: tau_ must be computed when: branch lengths
    computeNu();                                        //TODO: nu must be computed when: branch lengths + lambda + mu

    // Set all iotas
    setAllIotas();                                      //TODO: iotas must be computed when: branch lengths + mu

    // Set all betas
    setAllBetas();                                      //TODO: betas must be computed when: branch lengths + mu

    // Calls the routine to compute the FV values
    computeTreeLikelihood();

    minusLogLik_ = -getLogLikelihood();
}


void RHomogeneousTreeLikelihood_PIP::initialiseInsertionHistories() const {

    for (int i = 0; i < nbDistinctSites_; i++) {
        for (auto &nodeID:tree_->getNodesId()) {

            // Initialize vectors descCount_ and setA_ and indicatorFunctionVector
            std::vector<int> descCount_;
            std::vector<bool> setA_;
            descCount_.resize(nbDistinctSites_);

            setA_.resize(nbDistinctSites_);

            //if (tree_->isLeaf(nodeID))
            if (utree_->getNode(treemap_.left.at(nodeID))->isTerminalNode())
                indicatorFun_[nodeID].at(i).resize(nbStates_);

            descCountData_.insert(std::make_pair(nodeID, descCount_));
            setAData_.insert(std::make_pair(nodeID, setA_));
        }
    }
}


void RHomogeneousTreeLikelihood_PIP::setInsertionHistories(const SiteContainer &sites,
                                                           const std::vector<int> &nodeList,
                                                           std::map<int, std::vector<int>> *descCountData,
                                                           std::map<int, std::vector<bool>> *setAData) const {

    for (int i = 0; i < nbDistinctSites_; i++) {
        // Compute the total number of characters (exc. gap) for the current site
        int nonGaps_ = countNonGapCharacterInSite(i);
        for (auto &nodeID:nodeList) {
            // Computing descCount and descGapCount
            //if (tree_->isLeaf(nodeID)) {
            if (utree_->getNode(treemap_.left.at(nodeID))->isTerminalNode()) {
                (*descCountData)[nodeID][i] = (indicatorFun_[nodeID][i].back() == 1 ? 0 : 1);
            } else {
                auto _sons__ids = new std::vector<int>;// _sons__ids;
                //if (!tree_->isLeaf(nodeID))
                if (!utree_->getNode(treemap_.left.at(nodeID))->isTerminalNode())
                    //(*_sons__ids) = tree_->getSonsId(nodeID);
                    (*_sons__ids) = utree_->getChildrenId(treemap_.left.at(nodeID));

                // Reset value stored at the internal node
                (*descCountData)[nodeID][i] = 0;
                // Recompute it
                for (auto &l:(*_sons__ids))
                    (*descCountData)[nodeID][i] += (*descCountData)[l][i];

                delete _sons__ids;
            }
            // Activate or deactivate set A
            (*setAData)[nodeID][i] = ((*descCountData)[nodeID][i] == nonGaps_);
            DVLOG(3) << "setInsertionHistories [setA] (" << std::setfill('0')
                     << std::setw(3) << i << ") @node "
                     << tree_->getNode(nodeID)->getName()
                     << "\t" << (*setAData)[nodeID][i];
        }
    }
}

void RHomogeneousTreeLikelihood_PIP::setInsertionHistories(const SiteContainer &sites,
                                                           const std::vector<int> &nodeList,
                                                           std::map<int, std::vector<int>> *descCountData,
                                                           std::map<int, std::vector<bool>> *setAData,
                                                           std::map<int, bool> *ts_node__data_origin,
                                                           tshlib::Utree &_utree__topology) const {
    for (int i = 0; i < nbDistinctSites_; i++) {
        // Compute the total number of characters (exc. gap) for the current site
        int nonGaps_ = countNonGapCharacterInSite(i);

        for (auto &nodeID:nodeList) {
            // Reset value stored at the internal node
            (*descCountData)[nodeID][i] = 0;
            // Computing descCount and descGapCount
            //if (tree_->isLeaf(nodeID)) {
            if (utree_->getNode(treemap_.left.at(nodeID))->isTerminalNode()) {
                (*descCountData)[nodeID][i] = (indicatorFun_[nodeID][i].back() == 1 ? 0 : 1);
            } else {
                // Get the left and right son of the current node (mapped on the utree structure which can differ from tree_).
                std::vector<int> _sons__ids = _getMappedNodeChildren(nodeID, _utree__topology);
                for (auto &l:_sons__ids) {
                    if ((*ts_node__data_origin)[l]) {
                        (*descCountData)[nodeID][i] += (*descCountData)[l][i];
                    } else {
                        (*descCountData)[nodeID][i] += descCountData_[l][i];
                    }
                }
            }
            // Activate or deactivate set A
            (*setAData)[nodeID][i] = ((*descCountData)[nodeID][i] == nonGaps_);
            DVLOG(3) << "setInsertionHistories [setA] (" << std::setfill('0')
                     << std::setw(3) << i << ") @node "
                     << tree_->getNode(nodeID)->getName()
                     << "\t" << (*setAData)[nodeID][i];
        }
    }
}


int RHomogeneousTreeLikelihood_PIP::countNonGapCharacterInSite(int siteID) const {
    int nonGaps_ = 0;

    for (std::map<int, std::vector<std::vector<double>>>::iterator it = indicatorFun_.begin(); it != indicatorFun_.end(); ++it)
        nonGaps_ += (it->second.at(siteID).back() == 1 ? 0 : 1);


    return nonGaps_;
}


void RHomogeneousTreeLikelihood_PIP::setAllIotas() {

    Parameter mu_ = model_->getParameter("mu");

    //double _tree__totalLength = 0;
    //for(auto &brlen_name:getBranchLengthsParameters().getParameterNames())
    //    _tree__totalLength += getBranchLengthsParameters().getParameterValue(brlen_name);

    if (fabs(mu_.getValue()) < 1e-8)
        DLOG(WARNING) << "Constraint match at parameter mu, badValue = " << mu_.getValue() << " [ 1e-08; 10000]";


    // Compute tau value;
    double T_ = tree_->getTotalLength() + 1 / mu_.getValue();
    //double T_ = _tree__totalLength + 1 / mu_.getValue();

    if (fabs(T_) < 1e-8)
        DLOG(WARNING) << "Constraint match at parameter T, badValue = " << T_ << " [ 1e-08; 10000]";


    for (auto &nodeID:tree_->getNodesId()) {

        //if (!tree_->hasFather(nodeID)) {
        if (_utree__isRoot(nodeID)) {
            iotasData_[nodeID] = (1 / mu_.getValue()) / T_;
            continue;
        }
        //} else {
        iotasData_[nodeID] = tree_->getDistanceToFather(nodeID) / T_;

    }

}


void RHomogeneousTreeLikelihood_PIP::setAllBetas() {

    Parameter mu_ = model_->getParameter("mu");

    if (fabs(mu_.getValue()) < 1e-8)
        DLOG(WARNING) << "Constraint match at parameter mu, badValue = " << mu_.getValue() << " [ 1e-08; 10000]";


    for (auto &nodeID:tree_->getNodesId()) {

        //if (!tree_->hasFather(nodeID)) {
        if (_utree__isRoot(nodeID)) {
            betasData_[nodeID] = 1.0;
            continue;
        }

        double brLen__node = tree_->getDistanceToFather(nodeID);
        betasData_[nodeID] = (1.0 - exp(-mu_.getValue() * brLen__node)) / (mu_.getValue() * brLen__node);
            //node->vnode_beta = (1.0 - exp(-this->mu * vnode->vnode_branchlength)) / (this->mu * vnode->vnode_branchlength);

    }

}


double RHomogeneousTreeLikelihood_PIP::computePhi(double lkEmptyColumn) const {

    double p;
    double log_factorial_m;
    int m = (int) nbSites_;

    log_factorial_m = 0;
    for (int i = 1; i <= m; i++) {
        log_factorial_m += log(i);
    }

    p = -log_factorial_m + m * log(nu_) + (nu_ * (lkEmptyColumn - 1));

    return p;
}


void RHomogeneousTreeLikelihood_PIP::computeNu() {

    double lambda_ = model_->getParameter("lambda").getValue();
    double mu_ = model_->getParameter("mu").getValue();

    if (fabs(mu_) < 1e-8)
        DLOG(WARNING) << "Constraint match at parameter mu, badValue = " << mu_ << " [ 1e-08; 10000]";

    nu_ = lambda_ * (tau_ + 1 / mu_);

}


std::vector<int> RHomogeneousTreeLikelihood_PIP::getNodeListPostOrder(int startNodeID) const {
    std::vector<int> newList;
    getNodeListPostOrder_(newList, startNodeID);
    return newList;
}


void RHomogeneousTreeLikelihood_PIP::getNodeListPostOrder_(std::vector<int> &nodeList, int startNodeID) const {
    Node *node = tree_->getNode(startNodeID);
    if (!node->isLeaf()) {
        for (auto &son:node->getSons()) {
            getNodeListPostOrder_(nodeList, son->getId());
        }
        nodeList.push_back(startNodeID);
    } else {
        nodeList.push_back(startNodeID);
    }
}


void RHomogeneousTreeLikelihood_PIP::_printFV(Node *node, VVVdouble *likelihoodvector) const {

    // Print FV if required ------------------------------------------------------------------------------------------------
    std::ostringstream sout;
    std::ostringstream sout_siteline;

    unsigned long sites = (*likelihoodvector).size();
    unsigned long classes = (*likelihoodvector)[0].size();
    unsigned long states = (*likelihoodvector)[0][0].size();

    VVVdouble lknode_i_c = (*likelihoodvector);

    for (int c = 0; c < classes; c++) {

        DVLOG(3) << "Node: " << node->getName() << " - ASVR class [" << c << "]";

        for (int x = 0; x < states; x++) {

            for (int i = 0; i < sites; i++) {

                //unsigned long idxSite = likelihoodData_->getRootArrayPosition(i);

                sout_siteline << std::setfill('0') << std::setw(3) << i << "      | ";
                sout << std::fixed << std::setw(8) << std::setprecision(15) << (*likelihoodvector)[i][c][x] << " | ";
            }
            if (x == 0) {
                DVLOG(3) << sout_siteline.str();
                sout_siteline.str("");
                sout_siteline.clear();
            }
            DVLOG(3) << sout.str();
            sout.str("");
            sout.clear();
        }
    }
    DVLOG(3) << "-----------------------------------------------------";


}


void RHomogeneousTreeLikelihood_PIP::_printPrMatrix(Node *node, VVdouble *pr) {
    std::ostringstream sout;
    DVLOG(3) << "Node: " << node->getName() << " - distance [" << node->getDistanceToFather() << "]";

    for (size_t x = 0; x < nbStates_; x++) {
        for (size_t y = 0; y < nbStates_; y++) {

            sout << std::fixed << std::setw(8) << std::setprecision(6) << (*pr)[x][y] << "\t";

        }
        DVLOG(3) << sout.str();
        sout.str("");
        sout.clear();
        //DVLOG(3) << "\n";
    }
    DVLOG(3) << "---------------------------------";


}


void
RHomogeneousTreeLikelihood_PIP::_extendNodeListOnSetA(int qnodeID,
                                                      unsigned long site,
                                                      std::vector<int> &_node__list,
                                                      std::map<int, std::vector<bool>> *_ts__setadata,
                                                      tshlib::Utree &_utree__topology) const {

    _node__list.push_back(qnodeID);

    // Get corresponding tshlib::virtualnode
    //tshlib::VirtualNode *temp = treemap_.left.at(qnodeID);
    tshlib::VirtualNode *temp = _utree__topology.getNodeIdsMap().at(treemap_.left.at(qnodeID));

    do {

        if (temp->isTerminalNode()) {
            break;
        }

        if ((*_ts__setadata)[treemap_.right.at(temp->getNodeLeft()->getVnode_id())][site]) {

            temp = temp->getNodeLeft();

        } else if ((*_ts__setadata)[treemap_.right.at(temp->getNodeRight()->getVnode_id())][site]) {

            temp = temp->getNodeRight();

        } else {

            break;
        }

        //listNodes.push_back(tree_->getNode(treemap_.right.at(temp))->getId());
        _node__list.push_back(tree_->getNode(treemap_.right.at(temp->getVnode_id()))->getId());

    } while ((*_ts__setadata)[treemap_.right.at(temp->getVnode_id())][site]);

    std::reverse(_node__list.begin(), _node__list.end());

}


void RHomogeneousTreeLikelihood_PIP::_extendNodeListOnSetA(tshlib::VirtualNode *qnode,
                                                           std::vector<int> &listNodes,
                                                           unsigned long site,
                                                           tshlib::Utree &candUTree) const {

    tshlib::VirtualNode *temp = qnode;

    // Get the corresponding node on the BPP tree
    Node *node = tree_->getNode(treemap_.right.at(temp->getVnode_id()));
    listNodes.push_back(node->getId());

    do {

        if (temp->isTerminalNode()) {
            break;
        }

        if (setAData_[treemap_.right.at(temp->getNodeLeft()->getVnode_id())][site]) {

            temp = temp->getNodeLeft();

        } else if (setAData_[treemap_.right.at(temp->getNodeRight()->getVnode_id())][site]) {

            temp = temp->getNodeRight();

        } else {

            break;
        }

        listNodes.push_back(tree_->getNode(treemap_.right.at(temp->getVnode_id()))->getId());


    } while (setAData_[treemap_.right.at(temp->getVnode_id())][site]);

}


std::vector<int> RHomogeneousTreeLikelihood_PIP::remapVirtualNodeLists(std::vector<int> &inputList) const {

    std::vector<int> newList(inputList.size());

    for (int i = 0; i < inputList.size(); i++) {
        newList[i] = treemap_.right.at(inputList[i]);
    }

    return newList;
}


void RHomogeneousTreeLikelihood_PIP::setIndicatorFunction(const SiteContainer &sites) const {

    for (int i = 0; i < nbDistinctSites_; i++) {
        size_t indexRealSite = static_cast<size_t>(rootPatternLinksInverse_.at(i));

        for (auto &node:tree_->getNodes()) {
            if (node->isLeaf()) {

                size_t state = sites.getSequence(node->getName()).getValue(indexRealSite);
                std::vector<int> stateAliases = sites.getAlphabet()->getAlias(state);

                double stateValue = 1.0 / stateAliases.size();
                for (auto &salias : stateAliases)
                    indicatorFun_[node->getId()].at(i).at(salias) = stateValue;

            }
        }
    }

}


/********************************************************************************************************************************************/

void RHomogeneousTreeLikelihood_PIP::computeTreeLikelihood() {

    computeSubtreeLikelihood(likelihoodNodes_);
}

double RHomogeneousTreeLikelihood_PIP::computeLikelihoodWholeAlignmentEmptyColumn() const {

    // Add iota and beta quantities on nodes with actived SetA for empty column
    std::vector<int> _nodes__likelihood = utree_->getPostOrderNodeList(-1);
    std::vector<int> likelihoodNodes = remapVirtualNodeLists(_nodes__likelihood);
    //std::vector<int> likelihoodNodes = getNodeListPostOrder(tree_->getRootNode()->getId());
    std::vector<double> _node__partial_likelihood(likelihoodNodes.size());

    for (int idxNode = 0; idxNode < likelihoodNodes.size(); idxNode++) {

        auto _sons__ids = new std::vector<int>;
        std::vector<VVVdouble *> lk_sons_empty(2);
        //if (!tree_->isLeaf(likelihoodNodes[idxNode])) {
        //if (!utree_->getNode(treemap_.left.at(likelihoodNodes[idxNode]))->isTerminalNode()) {
        if (!_utree__isLeaf(likelihoodNodes[idxNode])) {
            //(*_sons__ids) = tree_->getSonsId(likelihoodNodes[idxNode]);
            //(*_sons__ids) = utree_->getChildrenId(treemap_.left.at(likelihoodNodes[idxNode]));
            (*_sons__ids) = utree_->getChildrenId(_utree__mapNode(likelihoodNodes[idxNode]));

            for (int l = 0; l < 2; l++) {
                lk_sons_empty[l] = &likelihoodEmptyData_->getLikelihoodArray((*_sons__ids)[l]);
            }
        }
        delete _sons__ids;
        _node__partial_likelihood[idxNode] = _kernel_likelihood_empty_forasite(likelihoodNodes[idxNode], &lk_sons_empty);
    }

    // Empty column likelihood
    double lk_site_empty = bpp::VectorTools::sum(_node__partial_likelihood);
    DVLOG(2) << "LK Empty [BPP] " << std::setprecision(18) << lk_site_empty;
    return lk_site_empty;
}


// this function is used during tree-search
double RHomogeneousTreeLikelihood_PIP::computeLikelihoodWholeAlignmentEmptyColumn(std::map<int, VVVdouble> *ts_lkemptydata,
                                                                                  std::map<int, bool> *ts_node__data_origin,
                                                                                  tshlib::Utree &_utree__topology) const {

    //std::vector<int> likelihoodNodes = getNodeListPostOrder(tree_->getRootNode()->getId());
    std::vector<int> _nodes__likelihood = _utree__topology.getPostOrderNodeList(-1);
    std::vector<int> likelihoodNodes = remapVirtualNodeLists(_nodes__likelihood);
    std::vector<double> _node__partial_likelihood(likelihoodNodes.size());

    for (int idxNode = 0; idxNode < likelihoodNodes.size(); idxNode++) {

        auto _sons__ids = new std::vector<int>;
        std::vector<VVVdouble *> lk_sons_empty(2);
        //if (!tree_->isLeaf(likelihoodNodes[idxNode])) {
        if (!utree_->getNode(treemap_.left.at(likelihoodNodes[idxNode]))->isTerminalNode()) {
            (*_sons__ids) = _getMappedNodeChildren(likelihoodNodes[idxNode], _utree__topology);

            for (int l = 0; l < 2; l++) {
                if ((*ts_node__data_origin)[(*_sons__ids)[l]]) {
                    lk_sons_empty[l] = &(*ts_lkemptydata)[(*_sons__ids)[l]];
                } else {
                    lk_sons_empty[l] = &likelihoodEmptyData_->getLikelihoodArray((*_sons__ids)[l]);
                }
            }
        }
        delete _sons__ids;
        _node__partial_likelihood[idxNode] = _kernel_likelihood_empty_forasite(likelihoodNodes[idxNode], &lk_sons_empty);

    }

    // Empty column likelihood
    double lk_site_empty = bpp::VectorTools::sum(_node__partial_likelihood);
    DVLOG(2) << "LK Empty [BPP] " << std::setprecision(18) << lk_site_empty;
    return lk_site_empty;
}

double RHomogeneousTreeLikelihood_PIP::_kernel_likelihood_empty_forasite(int nodeID, std::vector<VVVdouble *> *lk_sons_empty) const {

    double lk_site_empty = 0;
    std::vector<double> _cat__fvsite(nbClasses_, 0.);

    for (size_t c = 0; c < nbClasses_; c++) {

        //if (tree_->isLeaf(nodeID)) {
        if (utree_->getNode(treemap_.left.at(nodeID))->isTerminalNode()) {

            _cat__fvsite[c] += iotasData_[nodeID] * (1 - betasData_[nodeID]);

            DVLOG(3) << "lk empty at node[" << tree_->getNodeName(nodeID) << "] -> " << std::setprecision(18) << _cat__fvsite[c] << std::endl;
            continue;
        }

        std::vector<double> fvsons = _SingleRateCategoryHadamardMultFvSons(nodeID, 0, c, (*lk_sons_empty));

        //if (!tree_->isRoot(nodeID)) {
        if (treemap_.left.at(nodeID) != -1) {
            _cat__fvsite[c] += iotasData_[nodeID] * (1 - betasData_[nodeID] + betasData_[nodeID] * MatrixBppUtils::dotProd(&fvsons, &rootFreqs_));
        } else {
            _cat__fvsite[c] += iotasData_[nodeID] * betasData_[nodeID] * MatrixBppUtils::dotProd(&fvsons, &rootFreqs_);
        }

        DVLOG(3) << "lk empty at node[" << tree_->getNodeName(nodeID) << "] -> " << std::setprecision(18) << _cat__fvsite[c] << std::endl;

    }

    lk_site_empty = VectorTools::sum(VectorTools::kroneckerMult(_cat__fvsite, rateDistribution_->getProbabilities()));

    return lk_site_empty;

}

double RHomogeneousTreeLikelihood_PIP::computeLikelihoodForASite(std::vector<int> &_node__list, size_t i) const {
    std::vector<double> _node__partial_likelihood(_node__list.size());

    for (int idxNode = 0; idxNode < _node__list.size(); idxNode++) {

        auto _sons__ids = new std::vector<int>;
        std::vector<VVVdouble *> lk_sons(2);
        std::vector<VVVdouble *> lk_sons_empty(2);

        // If the node is internal, then fetch its children nodes and their data-structures
        //if (!tree_->isLeaf(_node__list[idxNode])) {
        //if (!utree_->getNode(treemap_.left.at(_node__list[idxNode]))->isTerminalNode()) {
        if (!_utree__isLeaf(_node__list[idxNode])) {
            //(*_sons__ids) = tree_->getSonsId(_node__list[idxNode]);
            //(*_sons__ids) = utree_->getChildrenId(treemap_.left.at(_node__list[idxNode]));
            (*_sons__ids) = utree_->getChildrenId(_utree__mapNode(_node__list[idxNode]));

            for (int l = 0; l < 2; l++) {
                lk_sons[l] = &likelihoodData_->getLikelihoodArray((*_sons__ids)[l]);
                lk_sons_empty[l] = &likelihoodEmptyData_->getLikelihoodArray((*_sons__ids)[l]);
            }
        }
        delete _sons__ids;
        // compute the likelihood at this node
        _node__partial_likelihood[idxNode] = _kernel_likelihood_forasite(i, _node__list[idxNode], &lk_sons, &lk_sons_empty, &setAData_);
    }

    return bpp::VectorTools::sum(_node__partial_likelihood);
}

// this function is used during tree-search
double RHomogeneousTreeLikelihood_PIP::computeLikelihoodForASite(size_t i,
                                                                 std::map<int, VVVdouble> *likelihoods,
                                                                 std::map<int, VVVdouble> *likelihoods_empty,
                                                                 std::map<int, std::vector<bool>> *ts_setadata,
                                                                 std::vector<int> &_node__list,
                                                                 std::map<int, bool> *ts_node__data_origin,
                                                                 tshlib::Utree &_utree__topology) const {

    std::vector<double> _node__partial_likelihood(_node__list.size(),0);

    for (int idxNode = 0; idxNode < _node__list.size(); idxNode++) {

        auto _sons__ids = new std::vector<int>;
        std::vector<VVVdouble *> lk_sons(2);
        std::vector<VVVdouble *> lk_sons_empty(2);
        
        //if (!tree_->isLeaf(_node__list[idxNode])) {
        if (!utree_->getNode(treemap_.left.at(_node__list[idxNode]))->isTerminalNode()) {
            (*_sons__ids) = _getMappedNodeChildren(_node__list[idxNode], _utree__topology);

            for (int l = 0; l < 2; l++) {
                if ((*ts_node__data_origin)[(*_sons__ids)[l]]) {
                    lk_sons[l] = &(*likelihoods)[(*_sons__ids)[l]];
                    lk_sons_empty[l] = &(*likelihoods_empty)[(*_sons__ids)[l]];
                } else {
                    lk_sons[l] = &likelihoodData_->getLikelihoodArray((*_sons__ids)[l]);
                    lk_sons_empty[l] = &likelihoodEmptyData_->getLikelihoodArray((*_sons__ids)[l]);
                }
            }
        }
        delete _sons__ids;
        _node__partial_likelihood[idxNode] = _kernel_likelihood_forasite(i, _node__list[idxNode], &lk_sons, &lk_sons_empty, ts_setadata);

    }

    return bpp::VectorTools::sum(_node__partial_likelihood);
}

double RHomogeneousTreeLikelihood_PIP::_kernel_likelihood_forasite(size_t i,
                                                                   int nodeID,
                                                                   std::vector<VVVdouble *> *lk_sons,
                                                                   std::vector<VVVdouble *> *lk_sons_empty,
                                                                   std::map<int, std::vector<bool>> *ts_setadata) const {

    double lk_site = 0;
    std::vector<double> _site__class_likelihood(nbClasses_,0);
    for (size_t c = 0; c < nbClasses_; c++) {
        double fv_site = 0;
        if ((*ts_setadata)[nodeID][i]) {
            DVLOG(3) << "[BPP] Likelihood for setA (" << i << ") @node " << tree_->getNodeName(nodeID);
            //if (!tree_->isLeaf(nodeID)) {
            if (!_utree__isLeaf(nodeID)) {
                std::vector<double> fvsons = _SingleRateCategoryHadamardMultFvSons(nodeID, i, c, (*lk_sons));
                fv_site += iotasData_[nodeID] * betasData_[nodeID] * MatrixBppUtils::dotProd(&fvsons, &rootFreqs_);
            } else {
                fv_site += (iotasData_[nodeID] * betasData_[nodeID]) * MatrixBppUtils::dotProd(&indicatorFun_[nodeID][i], &rootFreqs_);
            }
        }
        // Multiply the likelihood value of the site for the ASVR category probability
        //lk_site += (fv_site * rateDistribution_->getProbability(c));
        _site__class_likelihood[c] = (fv_site * rateDistribution_->getProbability(c));
    }
    return bpp::VectorTools::sum(_site__class_likelihood);
    //return lk_site;
}

/********************************************************************************************************************************************/

double RHomogeneousTreeLikelihood_PIP::getLogLikelihood() const {

    // Initialise lk
    double logLK;

    // Initialise vector of likelihood per each site
    std::vector<double> lk_sites(nbDistinctSites_);

    // compute likelihood empty column
    double lk_site_empty = computeLikelihoodWholeAlignmentEmptyColumn();

    const std::vector<unsigned int> *rootWeights = &likelihoodData_->getWeights();

    // For each site in the alignment
    for (int i = 0; i < nbDistinctSites_; i++) {

        lk_sites[i] = log(computeLikelihoodForASite(likelihoodNodes_, i)) * rootWeights->at(i);
        DVLOG(2) << "site log_lk[" << i << "]=" << std::setprecision(18) << lk_sites[i] << std::endl;

    }

    // Sum all the values stored in the lk vector
    logLK = MatrixBppUtils::sumVector(&lk_sites);
    DVLOG(2) << "LK Sites [BPP] " << std::setprecision(18) << logLK;

    // compute PHi
    double log_phi_value = computePhi(lk_site_empty);
    DVLOG(2) << "Log PHI [BPP] " << std::setprecision(18) << log_phi_value;

    // Add phi to site likelihood
    logLK += log_phi_value;

    return logLK;

}

/********************************************************************************************************************************************/

double RHomogeneousTreeLikelihood_PIP::getDLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {

    double dl = 0;
    Vdouble *dla = &likelihoodData_->getDLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
    for (size_t i = 0; i < nbStates_; i++) {
        dl += (*dla)[i] * rootFreqs_[i];
    }
    return dl;

}

double RHomogeneousTreeLikelihood_PIP::getDLikelihoodForASite(size_t site) const {
    // Derivative of the sum is the sum of derivatives:
    double dl = 0;

    for (size_t i = 0; i < nbClasses_; i++) {
        dl += getDLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
    }

    return dl;
}

double RHomogeneousTreeLikelihood_PIP::getDLogLikelihoodForASite(size_t site) const {
    // d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
    return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

double RHomogeneousTreeLikelihood_PIP::getDLogLikelihood() const {
    // Derivative of the sum is the sum of derivatives:
    double dl = 0;
    for (size_t i = 0; i < nbSites_; i++) {
        dl += getDLogLikelihoodForASite(i);
    }
    return dl;
}

void RHomogeneousTreeLikelihood_PIP::computeTreeDLikelihood(const std::string &variable) {
    // Get the node with the branch whose length must be derivated:
    size_t brI = TextTools::to<size_t>(variable.substr(5));
    const Node *branch = nodes_[brI];
    const Node *father = branch->getFather();
    VVVdouble *_dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());

    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    size_t nbSites = _dLikelihoods_father->size();
    for (size_t i = 0; i < nbSites; i++) {
        VVdouble *_dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++) {
            Vdouble *_dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
            for (size_t s = 0; s < nbStates_; s++) {
                (*_dLikelihoods_father_i_c)[s] = 1.;
            }
        }
    }

    size_t nbNodes = father->getNumberOfSons();
    for (size_t l = 0; l < nbNodes; l++) {
        const Node *son = father->getSon(l);

        std::vector<size_t> *_patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

        if (son == branch) {
            VVVdouble *dpxy__son = &dpxy_[son->getId()];
            for (size_t i = 0; i < nbSites; i++) {
                VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
                VVdouble *_dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
                for (size_t c = 0; c < nbClasses_; c++) {
                    Vdouble *_likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                    Vdouble *_dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
                    VVdouble *dpxy__son_c = &(*dpxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++) {
                        double dl = 0;
                        Vdouble *dpxy__son_c_x = &(*dpxy__son_c)[x];
                        for (size_t y = 0; y < nbStates_; y++) {
                            dl += (*dpxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
                        }
                        (*_dLikelihoods_father_i_c)[x] *= dl;
                    }
                }
            }
        } else {
            VVVdouble *pxy__son = &pxy_[son->getId()];
            for (size_t i = 0; i < nbSites; i++) {
                VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
                VVdouble *_dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
                for (size_t c = 0; c < nbClasses_; c++) {
                    Vdouble *_likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                    Vdouble *_dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
                    VVdouble *pxy__son_c = &(*pxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++) {
                        double dl = 0;
                        Vdouble *pxy__son_c_x = &(*pxy__son_c)[x];
                        for (size_t y = 0; y < nbStates_; y++) {
                            dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
                        }
                        (*_dLikelihoods_father_i_c)[x] *= dl;
                    }
                }
            }
        }
    }

    // Now we go down the tree toward the root node:
    computeDownSubtreeDLikelihood(father);
}

void RHomogeneousTreeLikelihood_PIP::computeDownSubtreeDLikelihood(const Node *node) {
    const Node *father = node->getFather();
    // We assume that the _dLikelihoods array has been filled for the current node 'node'.
    // We will evaluate the array for the father node.
    if (father == NULL) return; // We reached the root!

    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble *_dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());
    size_t nbSites = _dLikelihoods_father->size();
    for (size_t i = 0; i < nbSites; i++) {
        VVdouble *_dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++) {
            Vdouble *_dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
            for (size_t s = 0; s < nbStates_; s++) {
                (*_dLikelihoods_father_i_c)[s] = 1.;
            }
        }
    }

    size_t nbNodes = father->getNumberOfSons();
    for (size_t l = 0; l < nbNodes; l++) {
        const Node *son = father->getSon(l);

        VVVdouble *pxy__son = &pxy_[son->getId()];
        std::vector<size_t> *_patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

        if (son == node) {
            VVVdouble *_dLikelihoods_son = &likelihoodData_->getDLikelihoodArray(son->getId());
            for (size_t i = 0; i < nbSites; i++) {
                VVdouble *_dLikelihoods_son_i = &(*_dLikelihoods_son)[(*_patternLinks_father_son)[i]];
                VVdouble *_dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
                for (size_t c = 0; c < nbClasses_; c++) {
                    Vdouble *_dLikelihoods_son_i_c = &(*_dLikelihoods_son_i)[c];
                    Vdouble *_dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
                    VVdouble *pxy__son_c = &(*pxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++) {
                        double dl = 0;
                        Vdouble *pxy__son_c_x = &(*pxy__son_c)[x];
                        for (size_t y = 0; y < nbStates_; y++) {
                            dl += (*pxy__son_c_x)[y] * (*_dLikelihoods_son_i_c)[y];
                        }
                        (*_dLikelihoods_father_i_c)[x] *= dl;
                    }
                }
            }
        } else {
            VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
            for (size_t i = 0; i < nbSites; i++) {
                VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
                VVdouble *_dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
                for (size_t c = 0; c < nbClasses_; c++) {
                    Vdouble *_likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                    Vdouble *_dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
                    VVdouble *pxy__son_c = &(*pxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++) {
                        double dl = 0;
                        Vdouble *pxy__son_c_x = &(*pxy__son_c)[x];
                        for (size_t y = 0; y < nbStates_; y++) {
                            dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
                        }
                        (*_dLikelihoods_father_i_c)[x] *= dl;
                    }
                }
            }
        }
    }

    //Next step: move toward grand father...
    computeDownSubtreeDLikelihood(father);
}

double RHomogeneousTreeLikelihood_PIP::getFirstOrderDerivative(const std::string &variable) const throw(Exception) {
    if (!hasParameter(variable))
        throw ParameterNotFoundException("RHomogeneousTreeLikelihood::getFirstOrderDerivative().", variable);
    if (getRateDistributionParameters().hasParameter(variable)) {
        throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
    }
    if (getSubstitutionModelParameters().hasParameter(variable)) {
        throw Exception("Derivatives respective to substitution model parameters are not implemented.");
    }

    //const_cast<RHomogeneousTreeLikelihood_PIP*>(this)->computeTreeDLikelihood(variable);
    //double result = -getDLogLikelihood();
    double result = const_cast<RHomogeneousTreeLikelihood_PIP *>(this)->computeN1DerivativeLikelihood(variable);
    //auto der = new ThreePointsNumericalDerivative(const_cast<RHomogeneousTreeLikelihood_PIP *>(this));
    //double result2 = der->df(variable, this->getParameters());

    DVLOG(3) << "1D: " << result;
    return result;
}

double RHomogeneousTreeLikelihood_PIP::computeN1DerivativeLikelihood(const std::string &variable) {

    size_t brI = TextTools::to<size_t>(variable.substr(5));
    const Node *branch = nodes_[brI];

    double lk_1 = 0;
    double lk_2 = 0;
    double result = 0;
    double perturbation = 0.0001;

    double BranchLength_1 = branch->getDistanceToFather();
    double BranchLength_2 = (branch->getDistanceToFather() + perturbation);

    // Evaluate the likelihood function under the incoming parameter value
    lk_1 = -getValue();

    // Evaluate the likelihood function under the perturbed parameter value
    this->setParameterValue(variable, BranchLength_2);
    fireParameterChanged(getParameters());
    lk_2 = -getValue();

    // Restore previous branch length
    this->setParameterValue(variable, BranchLength_1);
    //fireParameterChanged(getParameters());

    // Compute the slope
    double diff_bl = BranchLength_2 - BranchLength_1;
    //double diff_lk = lk_2 - lk_1;

    //double op = lk_1 + log(1+exp(lk_2-lk_1));

    double op = lk_2 - lk_1;
    result = op / log(diff_bl);

    //VLOG(1) << "lk2 "<< std::setprecision(18) << lk_2 << " lk1 " << std::setprecision(18) << lk_1 << " op " << std::setprecision(18) << op << " ratio " << result ;

    // Save the first derivative
    d1bl_ = result;

    // Return the value to the first derivative
    return (d1bl_);


}

double RHomogeneousTreeLikelihood_PIP::evaluateLikelihoodPointForBranchDerivative(const std::string &variable, double new_branchlength) {

    double lk_2 = 0;

    // Evaluate the likelihood function under the perturbed parameter value
    this->setParameterValue(variable, new_branchlength);
    //fireParameterChanged(getParameters());
    lk_2 = -getValue();

    return lk_2;

}

/********************************************************************************************************************************************/

double RHomogeneousTreeLikelihood_PIP::getD2LikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {
    double d2l = 0;
    Vdouble *d2la = &likelihoodData_->getD2LikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
    for (size_t i = 0; i < nbStates_; i++) {
        d2l += (*d2la)[i] * rootFreqs_[i];
    }
    return d2l;
}

double RHomogeneousTreeLikelihood_PIP::getD2LikelihoodForASite(size_t site) const {
    // Derivative of the sum is the sum of derivatives:
    double d2l = 0;
    for (size_t i = 0; i < nbClasses_; i++) {
        d2l += getD2LikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
    }
    return d2l;
}

double RHomogeneousTreeLikelihood_PIP::getD2LogLikelihoodForASite(size_t site) const {
    return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
           - pow(getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

double RHomogeneousTreeLikelihood_PIP::getD2LogLikelihood() const {
    // Derivative of the sum is the sum of derivatives:
    double dl = 0;
    for (size_t i = 0; i < nbSites_; i++) {
        dl += getD2LogLikelihoodForASite(i);
    }
    return dl;
}

void RHomogeneousTreeLikelihood_PIP::computeTreeD2Likelihood(const std::string &variable) {
    // Get the node with the branch whose length must be derivated:
    size_t brI = TextTools::to<size_t>(variable.substr(5));
    const Node *branch = nodes_[brI];
    const Node *father = branch->getFather();

    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble *_d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
    size_t nbSites = _d2Likelihoods_father->size();
    for (size_t i = 0; i < nbSites; i++) {
        VVdouble *_d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++) {
            Vdouble *_d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
            for (size_t s = 0; s < nbStates_; s++) {
                (*_d2Likelihoods_father_i_c)[s] = 1.;
            }
        }
    }

    size_t nbNodes = father->getNumberOfSons();
    for (size_t l = 0; l < nbNodes; l++) {
        const Node *son = father->getSon(l);

        std::vector<size_t> *_patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

        if (son == branch) {
            VVVdouble *d2pxy__son = &d2pxy_[son->getId()];
            for (size_t i = 0; i < nbSites; i++) {
                VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
                VVdouble *_d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
                for (size_t c = 0; c < nbClasses_; c++) {
                    Vdouble *_likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                    Vdouble *_d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
                    VVdouble *d2pxy__son_c = &(*d2pxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++) {
                        double d2l = 0;
                        Vdouble *d2pxy__son_c_x = &(*d2pxy__son_c)[x];
                        for (size_t y = 0; y < nbStates_; y++) {
                            d2l += (*d2pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
                        }
                        (*_d2Likelihoods_father_i_c)[x] *= d2l;
                    }
                }
            }
        } else {
            VVVdouble *pxy__son = &pxy_[son->getId()];
            for (size_t i = 0; i < nbSites; i++) {
                VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
                VVdouble *_d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
                for (size_t c = 0; c < nbClasses_; c++) {
                    Vdouble *_likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                    Vdouble *_d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
                    VVdouble *pxy__son_c = &(*pxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++) {
                        double d2l = 0;
                        Vdouble *pxy__son_c_x = &(*pxy__son_c)[x];
                        for (size_t y = 0; y < nbStates_; y++) {
                            d2l += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
                        }
                        (*_d2Likelihoods_father_i_c)[x] *= d2l;
                    }
                }
            }
        }
    }

    // Now we go down the tree toward the root node:
    computeDownSubtreeD2Likelihood(father);
}

void RHomogeneousTreeLikelihood_PIP::computeDownSubtreeD2Likelihood(const Node *node) {
    const Node *father = node->getFather();
    // We assume that the _dLikelihoods array has been filled for the current node 'node'.
    // We will evaluate the array for the father node.
    if (father == NULL) return; // We reached the root!

    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble *_d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
    size_t nbSites = _d2Likelihoods_father->size();
    for (size_t i = 0; i < nbSites; i++) {
        VVdouble *_d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++) {
            Vdouble *_d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
            for (size_t s = 0; s < nbStates_; s++) {
                (*_d2Likelihoods_father_i_c)[s] = 1.;
            }
        }
    }

    size_t nbNodes = father->getNumberOfSons();
    for (size_t l = 0; l < nbNodes; l++) {
        const Node *son = father->getSon(l);

        VVVdouble *pxy__son = &pxy_[son->getId()];
        std::vector<size_t> *_patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

        if (son == node) {
            VVVdouble *_d2Likelihoods_son = &likelihoodData_->getD2LikelihoodArray(son->getId());
            for (size_t i = 0; i < nbSites; i++) {
                VVdouble *_d2Likelihoods_son_i = &(*_d2Likelihoods_son)[(*_patternLinks_father_son)[i]];
                VVdouble *_d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
                for (size_t c = 0; c < nbClasses_; c++) {
                    Vdouble *_d2Likelihoods_son_i_c = &(*_d2Likelihoods_son_i)[c];
                    Vdouble *_d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
                    VVdouble *pxy__son_c = &(*pxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++) {
                        double d2l = 0;
                        Vdouble *pxy__son_c_x = &(*pxy__son_c)[x];
                        for (size_t y = 0; y < nbStates_; y++) {
                            d2l += (*pxy__son_c_x)[y] * (*_d2Likelihoods_son_i_c)[y];
                        }
                        (*_d2Likelihoods_father_i_c)[x] *= d2l;
                    }
                }
            }
        } else {
            VVVdouble *_likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
            for (size_t i = 0; i < nbSites; i++) {
                VVdouble *_likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
                VVdouble *_d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
                for (size_t c = 0; c < nbClasses_; c++) {
                    Vdouble *_likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
                    Vdouble *_d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
                    VVdouble *pxy__son_c = &(*pxy__son)[c];
                    for (size_t x = 0; x < nbStates_; x++) {
                        double dl = 0;
                        Vdouble *pxy__son_c_x = &(*pxy__son_c)[x];
                        for (size_t y = 0; y < nbStates_; y++) {
                            dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
                        }
                        (*_d2Likelihoods_father_i_c)[x] *= dl;
                    }
                }
            }
        }
    }

    //Next step: move toward grand father...
    computeDownSubtreeD2Likelihood(father);
}

double RHomogeneousTreeLikelihood_PIP::computeN2DerivativeLikelihood(const std::string &variable) {

    size_t brI = TextTools::to<size_t>(variable.substr(5));
    const Node *branch = nodes_[brI];

    double lk_1 = 0;
    double lk_2 = 0;
    double dbl2 = 0;
    double perturbation = 0.0001;

    double BranchLength_1 = branch->getDistanceToFather();
    double BranchLength_2 = (branch->getDistanceToFather() + perturbation * 2);


    // Evaluate the likelihood function under the incoming parameter value
    lk_1 = -getValue();

    // Evaluate the second point of the function
    lk_2 = evaluateLikelihoodPointForBranchDerivative(variable, BranchLength_2);

    // Restore previous branch length
    this->setParameterValue(variable, BranchLength_1);
    //fireParameterChanged(getParameters());


    double op = lk_2 - lk_1;
    dbl2 = op / log(BranchLength_2 - BranchLength_1);


    d2bl_ = (dbl2 - d1bl_) / log(perturbation / 2);

    return d2bl_;
}

double RHomogeneousTreeLikelihood_PIP::getSecondOrderDerivative(const std::string &variable) const throw(Exception) {
    if (!hasParameter(variable))
        throw ParameterNotFoundException("RHomogeneousTreeLikelihood::getSecondOrderDerivative().", variable);
    if (getRateDistributionParameters().hasParameter(variable)) {
        throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
    }
    if (getSubstitutionModelParameters().hasParameter(variable)) {
        throw Exception("Derivatives respective to substitution model parameters are not implemented.");
    }

    //const_cast<RHomogeneousTreeLikelihood_PIP*>(this)->computeTreeD2Likelihood(variable);
    //double result = -getD2LogLikelihood();

    double result = const_cast<RHomogeneousTreeLikelihood_PIP *>(this)->computeN2DerivativeLikelihood(variable);

    DVLOG(3) << "2D: " << result;
    return result;
}

/********************************************************************************************************************************************/

void RHomogeneousTreeLikelihood_PIP::setLikelihoodNodes(std::vector<int> &nodeList) const {

    // Reset list containing the nodes on which computing the likelihood
    likelihoodNodes_.clear();

    // Overwrite the list with the one passed in the function arguments
    likelihoodNodes_ = nodeList;

}

void RHomogeneousTreeLikelihood_PIP::computeInDelDispersionOnTree(const SiteContainer &sites) {

    std::map<int, std::vector<std::string>> subsetAlignmentOnNode;
    std::map<int, double> localTreeLength;

    for (auto &nodeID:likelihoodNodes_) {
        // For each node in the tree (except the root node)
        //if (tree_->getNode(nodeID)->hasFather()) {
        if (tree_->getNode(nodeID)->isLeaf()) {
            // leaf node
            subsetAlignmentOnNode[nodeID].push_back(tree_->getNode(nodeID)->getName());
            localTreeLength[nodeID] = 0;
        } else {
            // internal node

            // Merge sons into internal node
            for (auto &sonID:tree_->getNode(nodeID)->getSonsId()) {
                for (auto &sonItem:subsetAlignmentOnNode[sonID]) {
                    subsetAlignmentOnNode[nodeID].push_back(sonItem);
                }
                localTreeLength[nodeID] += tree_->getNode(sonID)->getDistanceToFather() + localTreeLength[sonID];
            }

            tree_->getNode(nodeID)->setBranchProperty("local_tau",
                                                      *unique_ptr<Clonable>(new BppString(std::to_string((double) localTreeLength[nodeID]))));

            // get SubAlignment
            SiteContainer *subAlignment = getSubAlignment(sites, subsetAlignmentOnNode[nodeID]);

            // compute nh/ng
            setNhNgOnNode(*subAlignment, nodeID);

            delete subAlignment;

        }

    }
}

SiteContainer *RHomogeneousTreeLikelihood_PIP::getSubAlignment(const SiteContainer &sites, std::vector<std::string> subsetSequences) {

    auto subsetAlignment = new bpp::VectorSequenceContainer(sites.getAlphabet());
    for (auto &seqName:subsetSequences)
        subsetAlignment->addSequence(sites.getSequence(seqName));

    auto subsetSites = new bpp::VectorSiteContainer(*subsetAlignment);

    delete subsetAlignment;

    return subsetSites;
}

void RHomogeneousTreeLikelihood_PIP::setNhNgOnNode(SiteContainer &sites, int nodeID) {

    int gapCode = sites.getAlphabet()->getGapCharacterCode();

    int columnsWithGaps = 0;
    int columnsWithoutGaps = 0;
    // remove all the gap/unknown columns
    int exploredSites = 0;
    size_t originalSitesSize = sites.getNumberOfSites();
    size_t currSitePosition = 0;

    while (exploredSites < originalSitesSize) {

        int numGapsSeen = 0;
        int numCharSeen = 0;

        for (unsigned long s = 0; s < sites.getNumberOfSequences(); s++) {

            int currentChar = sites.getSite(currSitePosition).getValue(s);

            if (currentChar == gapCode) {
                numGapsSeen++;
            } else {
                numCharSeen++;
            }
        }

        if (numGapsSeen == sites.getSite(currSitePosition).size()) {
            sites.deleteSite(currSitePosition);

        } else {
            currSitePosition++;
        }
        exploredSites++;

    }

    size_t purgedAlignmentSize = sites.getNumberOfSites();

    for (unsigned long i = 0; i < sites.getNumberOfSites(); i++) {

        int numGapsSeen = 0;
        int numCharSeen = 0;

        for (unsigned long s = 0; s < sites.getNumberOfSequences(); s++) {

            int currentChar = sites.getSite(i).getValue(s);

            if (currentChar == gapCode) {
                numGapsSeen++;
            } else {
                numCharSeen++;
            }

        }
        // Count columns with gaps and columns without gaps
        if (numGapsSeen > 0) {
            columnsWithGaps++;

        } else {
            columnsWithoutGaps++;
        }

        //LOG_IF(FATAL, !nonGapSeen || !nonUnkownSeen) << "Column #" << i + 1 << " of the alignment contains only gaps. Please remove it and try again!";
    }

    tree_->getNode(nodeID)->setBranchProperty("node_age", *unique_ptr<Clonable>(new BppString(std::to_string((double) getNodeAge(nodeID)))));
    tree_->getNode(nodeID)->setBranchProperty("nh", *unique_ptr<Clonable>(new BppString(std::to_string(columnsWithoutGaps))));
    tree_->getNode(nodeID)->setBranchProperty("ng", *unique_ptr<Clonable>(new BppString(std::to_string(columnsWithGaps))));

    double ratio_nhng = 0;
    double weighted_nh = 0;
    double weighted_ng = 0;

    if (columnsWithoutGaps > 0 && columnsWithGaps > 0) {
        ratio_nhng = (double) columnsWithoutGaps / columnsWithGaps;
        weighted_nh = (double) columnsWithoutGaps / purgedAlignmentSize;
        weighted_ng = (double) columnsWithGaps / purgedAlignmentSize;
    }

    tree_->getNode(nodeID)->setBranchProperty("ratio_nhng", *unique_ptr<Clonable>(new BppString(std::to_string(ratio_nhng))));
    tree_->getNode(nodeID)->setBranchProperty("nh_weighted", *unique_ptr<Clonable>(new BppString(std::to_string(weighted_nh))));
    tree_->getNode(nodeID)->setBranchProperty("ng_weighted", *unique_ptr<Clonable>(new BppString(std::to_string(weighted_ng))));

    tree_->getNode(nodeID)->setBranchProperty("align_size", *unique_ptr<Clonable>(new BppString(std::to_string(purgedAlignmentSize))));


}

double RHomogeneousTreeLikelihood_PIP::getNodeAge(int nodeID) {
    double nodeAge = 0;

    int currNodeId = nodeID;

    while (tree_->getNode(currNodeId)->hasFather()) {

        nodeAge += tree_->getNode(currNodeId)->getDistanceToFather();
        currNodeId = tree_->getNode(currNodeId)->getFatherId();

    }

    return nodeAge;
}




