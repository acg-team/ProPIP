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
 * @file progressivePIP.cpp
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

#include <chrono>
#include <random>
#include <cfloat>

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <glog/logging.h>

#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"
#include "CompositePIPnode.hpp"
#include "CompositePIPnodeTBB.hpp"

using namespace bpp;

progressivePIP::progressivePIP(tshlib::Utree *utree,
                               bpp::Tree *tree,
                               bpp::SubstitutionModel *smodel,
                               UtreeBppUtils::treemap &inTreeMap,
                               bpp::SequenceContainer *sequences,
                               bpp::DiscreteDistribution *rDist,
                               long seed) {

    utree_ = utree; // tshlib tree
    tree_ = new TreeTemplate<Node>(*tree);  // bpp tree
    substModel_ = smodel;   // substitution model
    treemap_ = inTreeMap;   // tree map
    sequences_ = sequences; // sequences
    rDist_ = rDist; // rate-variation among site distribution
    alphabet_ = substModel_->getAlphabet(); // alphabet
    alphabetSize_ = alphabet_->getSize() - 1;   // original alphabet size
    extendedAlphabetSize_ = alphabet_->getSize();   // extended alphabet size
    seed_ = seed;   // seed for random number generation
};

void progressivePIP::_setLambda(double lambda) {

    // original lambda w/o rate variation
    lambda0_ = lambda;

    // insertion rate with rate variation among site r
    lambda_.resize(numCatg_);
    for (int i = 0; i < rDist_->getNumberOfCategories(); i++) {
        // lambda(r) = lambda * r
        lambda_.at(i) = lambda * rDist_->getCategories().at(i);
    }

}

void progressivePIP::_setMu(double mu) {

    // checks division by 0 or very small value
    if (fabs(mu) < SMALL_DOUBLE) {
        PLOG(FATAL) << "ERROR: mu is too small";
    }

    // original mu w/o rate variation
    mu0_ = mu;

    // deletion rate with rate variation among site r
    mu_.resize(numCatg_);
    for (int i = 0; i < rDist_->getCategories().size(); i++) {
        // mu(r) = mu *r
        mu_.at(i) = mu * rDist_->getCategories().at(i);
    }

}

void progressivePIP::_setPi(const Vdouble &pi) {

    // copy pi (steady state frequency distribution)
    // pi is a colMatrix (column array) to simplify the matrix multiplication
    pi_.resize(pi.size(), 1);
    for (int i = 0; i < pi.size(); i++) {
        pi_(i, 0) = pi.at(i);
    }

}

void progressivePIP::_setAllIotas() {

    // compute all the insertion probabilities (iota function)
    // N.B. the insertion probability at the root is different
    // from any other internal node

    double T;

    for(auto &node:compositePIPaligner_->pip_nodes_){

        node->iotasNode_.resize(numCatg_);

        if(node->_isRootNode()){

            // formula for the root PIPnode

            for (int catg = 0; catg < numCatg_; catg++) {

                // checks division by 0 or too small number
                if (fabs(mu_.at(catg)) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR in _setAllIotas: mu too small";
                }

                // T(r) = tau + 1/ (mu * r)
                T = tau_ + 1 / mu_.at(catg);

                // checks division by 0 or too small number
                if (fabs(T) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR in _setAllIotas: T too small";
                }

                // iota(root,r) = (lambda * r)/ (mu * r) / (lambda * r * (tau + 1/ (mu * r) ) )
                //              = 1 / (mu * r) / (tau + 1/ (mu *r) )
                node->iotasNode_.at(catg) =  (1 / mu_.at(catg)) / T;

            }

        }else{

            // formula for an internal PIPnode

            for (int catg = 0; catg < numCatg_; catg++) {


                // checks division by 0 or too small number
                if (fabs(mu_.at(catg)) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR in _setAllIotas: mu too small";
                }

                // T(r) = tau + 1/(mu * r)
                T = tau_ + 1 / mu_.at(catg);

                // checks division by 0 or too small number
                if (fabs(T) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR in _setAllIotas: T too small";
                }

                // checks division by 0 or too small number
                if (fabs(node->bnode_->getDistanceToFather()) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR in _setAllIotas: b(v) too small";
                }

                //iotasNode_[node->getId()][catg] = (node->getDistanceToFather() * rDist_->getCategory(catg) ) / T;
                // iota(v,r) = ( lambda * r * b(v) ) / (lambda * r * (tau + 1/ (mu *r) ) )
                //           = b(v) / (tau + 1/ (mu *r) )
                node->iotasNode_.at(catg) = node->bnode_->getDistanceToFather() / T;

            }

        }
    }

}

void progressivePIP::_setAllAlphas() {

    // alpha(v) = sum_from_v_to_root ( iota * beta * zeta )
    // zeta = exp(- mu *b ) is the "pure" survival probability

    for(auto &node:compositePIPaligner_->pip_nodes_){

        // resize for each gamma category
        node->alphaNode_.resize(numCatg_);

        for(int catg=0;catg<numCatg_;catg++){

            // checks division by 0 or too small number
            if (fabs(mu_.at(catg)) < SMALL_DOUBLE) {
                PLOG(WARNING) << "ERROR in _setAllAlphas: mu too small";
            }

            // initialize alpha with the probability at the starting node which is
            // alpha(v) = iota * beta
            node->alphaNode_.at(catg) = (1.0/mu_.at(catg)) / (tau_ + 1/mu_.at(catg));
        }

    }

}

void progressivePIP::_setAllEtas() {

    for(auto &node:compositePIPaligner_->pip_nodes_){

        // resize for each gamma category
        node->etaNode_.resize(numCatg_);

        for(int catg=0;catg<numCatg_;catg++) {

            if (node->_isRootNode()) {
                node->etaNode_.at(catg) = 0.0;
            } else {

                // checks division by 0 or too small number
                if (fabs(mu_.at(catg)) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR in _setAllEtas: mu too small";
                }

                // checks division by 0 or too small number
                if (fabs(node->bnode_->getDistanceToFather()) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR in _setAllEtas: b(v) too small";
                }

                node->etaNode_.at(catg) = node->alphaNode_.at(catg) * (1 - exp(-mu_.at(catg) * node->distanceToRoot_)) +
                                          (node->distanceToRoot_ / (tau_ + 1 / mu_.at(catg))) * \
                                          (1 - (1 - exp(-mu_.at(catg) * node->distanceToRoot_)) /
                                               (mu_.at(catg) * node->distanceToRoot_));
            }
        }

    }

}

void progressivePIP::_setAllBetas() {

    // compute all the survival probabilities (beta function)
    // N.B. the survival probability at the root is different
    // from any other internal node

    for(auto &node:compositePIPaligner_->pip_nodes_){

        node->betasNode_.resize(numCatg_);

        if(node->_isRootNode()){

            // formula for the root PIPnode

            for (int catg = 0; catg < rDist_->getNumberOfCategories(); catg++) {
                // by definition at the root (ev. local root) the survival probability is 1
                node->betasNode_.at(catg) = 1.0;
            }

        }else{

            // formula for an internal PIPnode

            for (int catg = 0; catg < rDist_->getCategories().size(); catg++) {

                // muT(r) = r * mu * b(v)
                double muT = rDist_->getCategory(catg) * mu_.at(catg) * node->bnode_->getDistanceToFather();

                // checks division by 0 or too small value
                if (fabs(muT) < SMALL_DOUBLE) {
                    PLOG(WARNING) << "ERROR mu * T is too small";
                }

                // survival probability on node v (different from (local)-root)
                // beta(v,r) = (1 - exp( -mu * r * b(v) )) / (mu * r * b(v))
                node->betasNode_.at(catg) = (1.0 - exp(-muT)) / muT;
            }

        }
    }

}

void progressivePIP::_computeLengthPathToRoot(){

    // get through the list of all the PIPnodes
    for (auto &node:compositePIPaligner_->pip_nodes_) {

        // get the bnode associated to the PIPnode
        bpp::Node *bnode = node->bnode_;

        // initialize the path length
        double T = 0.0;

        // climb the PIPnode tree
        // the root skips this loop
        while(bnode->hasFather()){

            // sum the branch length
            T += bnode->getDistanceToFather();

            // get the parent bnode
            bnode = bnode->getFather();
        }

        // save the path length at the given node
        // the path length from root to root is 0.0
        node->distanceToRoot_ = T;

    }

}

void progressivePIP::_buildPIPnodeTree() {

    // build a binary tree of PIPnode that dictates the alignmanet order

    // this method build a PIPnode binary tree with the same structure (same relations)
    // to the bpp tree

    for (auto &node:compositePIPaligner_->pip_nodes_) {

        bpp::Node *bnode = node->bnode_;

        if(bnode->hasFather()){
            // internal node
            int bnodeFatherId = bnode->getFatherId();
            node->parent_ = compositePIPaligner_->pip_nodes_.at(bnodeFatherId);
        }else{
            // root node
            node->parent_ = nullptr;
            PIPnodeRoot_ = node;
        }

        if(!bnode->isLeaf()){
            // internal node
            std::vector<int> bnodeSonsId = bnode->getSonsId();
            node->childL_ = compositePIPaligner_->pip_nodes_.at(bnodeSonsId.at(LEFT));
            node->childR_ = compositePIPaligner_->pip_nodes_.at(bnodeSonsId.at(RIGHT));
        }else{
            // leaf node
            node->childL_ = nullptr;
            node->childR_ = nullptr;
        }

    }

}

void progressivePIP::_computeTauRec_(PIPnode *pipnode) {

    // recursive computation of the total tree length and
    // of the total left/right subtree length

    if(pipnode->_isTerminalNode()){

        pipnode->subTreeLenL_ = 0.0;
        pipnode->subTreeLenR_ = 0.0;

    }else{

        _computeTauRec_(pipnode->childL_);
        _computeTauRec_(pipnode->childR_);

        pipnode->subTreeLenL_ = pipnode->childL_->subTreeLenL_ + \
                                pipnode->childL_->subTreeLenR_ + \
                                pipnode->childL_->bnode_->getDistanceToFather();

        pipnode->subTreeLenR_ = pipnode->childR_->subTreeLenL_ + \
                                pipnode->childR_->subTreeLenR_ + \
                                pipnode->childR_->bnode_->getDistanceToFather();

        tau_ += pipnode->childL_->bnode_->getDistanceToFather() +\
                pipnode->childR_->bnode_->getDistanceToFather();

    }

}

void progressivePIP::_computeTau_() {

    // compute the total tree length and the length of the left/right subtree
    // of the tree rooted at a given PIPnode

    // get the root Id
    int rootId = getBPProotNode()->getId();

    // get the root node
    PIPnode * rootPIPnode = compositePIPaligner_->pip_nodes_.at(rootId);

    // init the total tree length
    tau_ = 0.0;

    // compute the total tree length by means of a recursive method
    progressivePIP::_computeTauRec_(rootPIPnode);

}

void progressivePIP::_computeNu() {

    // compute the normalizing Poisson intensity (expected MSA length)

    // get the number of gamma categories
    nu_.resize(numCatg_);

    for (int catg = 0; catg < numCatg_; catg++) {


        // checks division by 0 or too small value
        if (fabs(mu_.at(catg)) < SMALL_DOUBLE) {
            PLOG(WARNING) << "ERROR mu is too small";
        }

        // computes the normalizing constant with discrete rate variation (gamma distribution)
        // nu(r) = lambda * r * (tau + 1/(mu *r))
        nu_.at(catg) = lambda_.at(catg) * (tau_ + 1 / mu_.at(catg));
    }

}

void progressivePIP::_initializePIP(std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                                    enumDP3Dversion DPversion,
                                    int num_sb,
                                    double temperature) {

    //***************************************************************************************
    // get dimensions
    numNodes_ = list_vnode_to_root.size(); // total number of nodes
    numCatg_ = rDist_->getNumberOfCategories(); // number of gamma categories
    num_sb_ = num_sb; // number of sub-optimal solutions
    temperature_ = temperature; // for SB algorithm
    //***************************************************************************************
    // computes lambda and mu with gamma
    // set lambdas with rate variation (gamma distribution)
    _setLambda(substModel_->getParameter("lambda").getValue());
    // set mus with rate variation (gamma distribution)
    _setMu(substModel_->getParameter("mu").getValue());
    //***************************************************************************************
    // set Pi
    // local copy of steady state frequency (Pi)
    _setPi(substModel_->getFrequencies());
    //***************************************************************************************

    nodeFactory *nodeFactory = new bpp::nodeFactory(); // Factory pattern for DP3D-CPU, DP3D-RAM,...

    switch (DPversion) {
        // TBB parallel_for cases
        case TBB_FOR_TASK:
            compositePIPaligner_ = new CompositePIPnodeTBB(numNodes_, true, true); // Composite pattern with array of PIPnodes
            break;
        case TBB_FOR_TASKBLOCK:
            compositePIPaligner_ = new CompositePIPnodeTBB(numNodes_, true, true); // Composite pattern with array of PIPnodes
            break;
        case TBB_FOR_BLOCK:
            compositePIPaligner_ = new CompositePIPnodeTBB(numNodes_, false, true); // Composite pattern with array of PIPnodes
            break;
            // TBB no parallel_for cases
        case TBB_TASK:
            compositePIPaligner_ = new CompositePIPnodeTBB(numNodes_, true, false); // Composite pattern with array of PIPnodes
            break;
        case TBB_TASKBLOCK:
            compositePIPaligner_ = new CompositePIPnodeTBB(numNodes_, true, false); // Composite pattern with array of PIPnodes
            break;
        case TBB_BLOCK:
            compositePIPaligner_ = new CompositePIPnodeTBB(numNodes_, false, false); // Composite pattern with array of PIPnodes
            break;
        default:
            compositePIPaligner_ = new CompositePIPnode(numNodes_); // Composite pattern with array of PIPnodes
            break;
    }

    for (auto &vnode:list_vnode_to_root) {

        // @ BOOST_ERROR
        //auto bnode = tree_->getNode(treemap_.right.at(vnode), false); // get bnode from vnode through the tree-map
        auto bnode = tree_->getNode(treemap_.right.at(vnode->getVnode_id()), false); // get bnode from vnode through the tree-map

        // create a PIPnode
        PIPnode * pip_node = nodeFactory->getPIPnode(DPversion, // PIPnode of type CPU, RAM, SB,... to access the correct DP version
                                                     this,      // PIPnode has access to progressivePIP fields through this pointer
                                                     vnode,     // PIPnode store the correponding vnode and
                                                     bnode);    // the bnode
        //***************************************************************************************
        // get Qs
        // set substitution/deletion probabilities with rate variation (gamma distribution)
        pip_node->_getPrFromSubstitutionModel();
        //***************************************************************************************

        compositePIPaligner_->addPIPnode(pip_node); // add PIPnode to composite array of PIPnodes

    }
    //***************************************************************************************

    _buildPIPnodeTree(); // build a tree of PIPnodes

    _computeTau_(); // compute the total tree length and the left/right subtree length
    // (length of the left/right subtree rooted at the given node) at each PIPnode

    _computeLengthPathToRoot();

    _computeNu();   // compute the Poisson normalizing intensity (corresponds to the expected MSA length)

    _setAllIotas(); // set iota (survival probability) on all nodes

    _setAllBetas(); // set beta (survival probability) on all nodes

    _setAllAlphas(); // alpha(v) = sum_from_v_to_root ( iota * beta * zeta )
    // zeta = exp(- mu *b ) is the "pure" survival probability

    _setAllEtas();
}

void progressivePIP::setSTFT_size(int size)
{
    for (auto &node : compositePIPaligner_->pip_nodes_)
    {
        if (node->_isTerminalNode() == false)
            node->_setSTFT_size(size);
    }
}


void progressivePIP::PIPnodeAlign(){

    // wrapper

    compositePIPaligner_->PIPnodeAlign();

}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
double progressivePIPutils::add_lns(double a_ln, double b_ln) {
    //ln(a + b) = ln{exp[ln(a) - ln(b)] + 1} + ln(b)

    double R;
    const double exp_precision =  log(pow(2,(double)DBL_MANT_DIG-1)-1);

    //ApplicationTools::displayResult("Mantissa precision", TextTools::toString(exp_precision, 50));

    if (std::isinf(a_ln) && std::isinf(b_ln)) {
        R = -std::numeric_limits<double>::infinity();
    } else if (std::isinf(a_ln)) {
        R = b_ln;
    } else if (std::isinf(b_ln)) {
        R = a_ln;
    } else if ((abs(a_ln - b_ln) >= exp_precision)) {
        //TODO:check this
        //2^52-1 = 4503599627370495.	log of that is 36.043653389117155867651465390794
        R = max(a_ln, b_ln);
    } else {
        R = log(exp(a_ln - b_ln) + 1) + b_ln;
    }

    return R;
}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************