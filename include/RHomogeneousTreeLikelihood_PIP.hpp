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
 * @file RHomogeneousTreeLikelihood_PIP.hpp
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
#ifndef CASTOR_RHOMOGENEOUSTREELIKELIHOOD_PIP_HPP
#define CASTOR_RHOMOGENEOUSTREELIKELIHOOD_PIP_HPP

#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRASRTreeLikelihoodData.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>

using namespace bpp;

#include "Utils.hpp"

namespace bpp {

    class RHomogeneousTreeLikelihood_PIP : public AbstractHomogeneousTreeLikelihood {
    protected:

        // Likelihood quantities (FV)
        mutable DRASRTreeLikelihoodData *likelihoodData_;           // Structure containing the 3D FV arrays at each internal node
        mutable DRASRTreeLikelihoodData *likelihoodEmptyData_;      // Structure containing the 3D FV arrays (empty column) at each internal node
        mutable std::vector<int> likelihoodNodes_;                  // The node is represented via its <int> ID

        // Insertion histories
        mutable std::map<int, std::vector<int>> descCountData_;     // Descendant count
        mutable std::map<int, std::vector<bool>> setAData_;         // SetA flags if a node should be included in the insertion histories

        // PIP quantities
        mutable std::map<int, double> iotasData_;                               //
        mutable std::map<int, double> betasData_;                               //
        mutable std::map<int, std::vector<std::vector<double>>> indicatorFun_;  //
        mutable double nu_;                                                     //
        mutable double tau_;                                                    //

        mutable UtreeBppUtils::treemap treemap_;                                // Bidirectional map Utree <-> BppTree
        mutable std::vector<unsigned long> rootPatternLinksInverse_;            //

        mutable std::map<int, std::map<int, std::vector<int>>> tsTemp_descCountData_;   // TS: collector for Descendant count vectors
        mutable std::map<int, std::map<int, std::vector<bool>>> tsTemp_setAData_;       // TS: collector for SetA flags vectors insertion hist.
        mutable std::map<int, std::map<int, bool>> tsTemp_node_data_origin;             // TS: origin of the data for the node

        mutable tshlib::Utree *utree_;

        double minusLogLik_;        // Log Likelihood values

        double d1bl_;               // value of the First-Derivative of the likelihood function
        double d2bl_;               // value of the Second-Derivative of the likelihood function

    public:
        /**
         * @brief Build a new RHomogeneousTreeLikelihood object without data.
         *
         * This constructor only initialize the parameters.
         * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
         *
         * @param tree The tree to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @param usePatterns Tell if recursive site compression should be performed.
         * @throw Exception in an error occured.
         */
        RHomogeneousTreeLikelihood_PIP(
                const Tree &tree,
                tshlib::Utree *utree,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                UtreeBppUtils::treemap *tm,
                bool checkRooted = true,
                bool verbose = true,
                bool usePatterns = true)
        throw(Exception);

        /**
         * @brief Build a new RHomogeneousTreeLikelihood object with data.
         *
         * This constructor initializes all parameters, data, and likelihood arrays.
         *
         * @param tree The tree to use.
         * @param data Sequences to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @param usePatterns Tell if recursive site compression should be performed.
         * @throw Exception in an error occured.
         */
        RHomogeneousTreeLikelihood_PIP(
                const Tree &tree,
                tshlib::Utree *utree,
                const SiteContainer &data,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                UtreeBppUtils::treemap *tm,
                bool checkRooted = true,
                bool verbose = true,
                bool usePatterns = true)
        throw(Exception);

        RHomogeneousTreeLikelihood_PIP(const RHomogeneousTreeLikelihood_PIP &lik);

        RHomogeneousTreeLikelihood_PIP &operator=(const RHomogeneousTreeLikelihood_PIP &lik);

        virtual ~RHomogeneousTreeLikelihood_PIP();

        RHomogeneousTreeLikelihood_PIP *clone() const { return new RHomogeneousTreeLikelihood_PIP(*this); }


    protected:

        /**
         * @brief Method called by constructors.
         */
        void init_(bool usePatterns) noexcept(false);

        void _computeHadamardFVSons(std::vector<VVVdouble *> inFVSons, VVVdouble *outFVParent) const;

        //std::vector<double> _SingleRateCategoryHadamardMultFvSons(int nodeID, unsigned long site, unsigned long rate) const;

        //std::vector<double> _SingleRateCategoryHadamardMultFvEmptySons(int nodeID, unsigned long rate) const;

        std::vector<double> _SingleRateCategoryHadamardMultFvSons(int nodeID,
                                                                  unsigned long site,
                                                                  unsigned long rate,
                                                                  std::vector<VVVdouble *> lk_sons) const;

        void _computePrTimesFv(VVVdouble *pxy__node, VVVdouble *_likelihoods_node) const;

        void _computePrTimesIndicator(VVVdouble *pxy__node, VVdouble *indicator_node, VVVdouble *_likelihoods_node) const;

        void _computePrTimesIndicatorEmpty(VVVdouble *pxy__node, VVdouble *indicator_node, VVVdouble *_likelihoods_node) const;

        std::vector<int> _getMappedNodeChildren(int nodeID, tshlib::Utree &_utree__topology) const;

        bool _utree__isLeaf(int nodeID) const { return utree_->getNode(_utree__mapNode(nodeID))->isTerminalNode(); };

        bool _utree__isRoot(int nodeID) const { return utree_->getNode(_utree__mapNode(nodeID))->isRootNode(); };

        int _utree__mapNode(int nodeID) const { return treemap_.left.at(nodeID); };

        void initialiseInsertionHistories() const;

    public:

        /**
         * @name The TreeLikelihood interface.
         *
         * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
         *
         * @{
         */
        void setData(const SiteContainer &sites) throw(Exception);

        double getNu() const { return nu_; }

        UtreeBppUtils::treemap &getTreemap() { return treemap_; }

        double getLikelihood() const {
            std::cerr << "getLikelihood()" << std::endl;
            return 0;
        };

        double getLogLikelihood() const;


        //double getLogLikelihood(std::vector<tshlib::VirtualNode *> &listNodes) const;
        //double getLogLikelihoodSubtree(const Node *node) const;
        //double getLogLikelihoodSubtreeForASite(size_t site) const;
        //double getLogLikelihoodSubtreeForASiteForARateClass(size_t site, size_t rateClass) const;


        /** @} */


        /**
         * @name The DiscreteRatesAcrossSites interface implementation:
         * @deprecated at the moment these methods are not used
         * @{
         */
        double getLikelihoodForASite(size_t site) const {
            std::cerr << "getLikelihoodForASite()" << std::endl;
            return 0;
        };

        double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {
            std::cerr << "getLikelihoodForASiteForARateClass()" << std::endl;
            return 0;
        };

        double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const {
            std::cerr << "getLikelihoodForASiteForARateClassForAState()" << std::endl;
            return 0;
        };

        double getLogLikelihoodForASite(size_t site) const {
            std::cerr << "getLogLikelihoodForASite()" << std::endl;
            return 0;
        };

        double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {
            std::cerr << "getLogLikelihoodForASiteForARateClass()" << std::endl;
            return 0;
        };

        double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const {
            std::cerr << "getLogLikelihoodForASiteForARateClassForAState()" << std::endl;
            return 0;
        };

        /** @} */

        /**
         * @brief Implements the Function interface.
         *
         * Update the parameter list and call the applyParameters() method.
         * Then compute the likelihoods at each node (computeLikelihood() method)
         * and call the getLogLikelihood() method.
         *
         * If a subset of the whole parameter list is passed to the function,
         * only these parameters are updated and the other remain constant (i.e.
         * equal to their last value).
         *
         * @param parameters The parameter list to pass to the function.
         */
        void setParameters(const ParameterList &parameters) throw(ParameterNotFoundException, ConstraintException);

        double getValue() const throw(Exception);

        size_t getSiteIndex(size_t site) const throw(IndexOutOfBoundsException) { return likelihoodData_->getRootArrayPosition(site); }

        /**
         * @name DerivableFirstOrder interface.
         *
         * @{
         */
        double getFirstOrderDerivative(const std::string &variable) const throw(Exception);

        double computeN1DerivativeLikelihood(const std::string &variable);

        double computeN2DerivativeLikelihood(const std::string &variable);

        double evaluateLikelihoodPointForBranchDerivative(const std::string &variable, double new_branchlength);
        /** @} */

        /**
         * @name DerivableSecondOrder interface.
         *
         * @{
         */
        double getSecondOrderDerivative(const std::string &variable) const throw(Exception);

        double getSecondOrderDerivative(const std::string &variable1,
                                        const std::string &variable2) const throw(Exception) { return 0; } // Not implemented for now.
        /** @} */

    public:    // Specific methods:


        DRASRTreeLikelihoodData *getLikelihoodData() { return likelihoodData_; }

        const DRASRTreeLikelihoodData *getLikelihoodData() const { return likelihoodData_; }

        std::vector<int> getNodeDescCounts(bpp::Node *node, int siteId) { return descCountData_[node->getId()]; }

        int getNodeDescCountForASite(int nodeID, int siteId) const { return descCountData_[nodeID][siteId]; }

        bool getSetAForANodeForASite(bpp::Node *node, int siteId) { return setAData_[node->getId()][siteId]; }


        /**
         * @name Interface to compute the likelihood components
         * @{
         */

        /**
         * @brief This method computes the likelihood of the tree  generating automatically postorder-traversal node list
         */
        void computeTreeLikelihood();

        /**
         * @brief This method computes the likelihood of the tree for a list of nodes computed using a postorder-traversal
         * @param nodeList
         */
        //void computeTreeLikelihood(std::vector<int> nodeList);

        /** @} */

        virtual double getDLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        virtual double getDLikelihoodForASite(size_t site) const;

        virtual double getDLogLikelihoodForASite(size_t site) const;

        virtual double getDLogLikelihood() const;

        virtual void computeTreeDLikelihood(const std::string &variable);

        virtual double getD2LikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        virtual double getD2LikelihoodForASite(size_t site) const;

        virtual double getD2LogLikelihoodForASite(size_t site) const;

        virtual double getD2LogLikelihood() const;

        virtual void computeTreeD2Likelihood(const std::string &variable);

        /**
         * @brief This method computes the likelihood after a tree rearrangment
         * @return The likelihood value using the intermediate partial values
         */
        //void fireTopologyChange(std::vector<int> nodeList);

        //double getLogLikelihoodOnTopologyChange() const;

        /**
         * @brief This method computes a list of nodes traversing the tree in postorder
         *
         */
        std::vector<int> getNodeListPostOrder(int startNodeID) const;

        void getNodeListPostOrder_(std::vector<int> &nodeList, int startNodeID) const;

        void setLikelihoodNodes(std::vector<int> &nodeList) const;

    protected:

        /**
         * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
         *        This function should be called only for filling the likelihood arrays (i.e. first traversal, parameter change -- but not topology changes).
         *
         * @param node The root of the subtree.
         */
        virtual void computeSubtreeLikelihood(const std::vector<int> &nodeList);

        virtual void computeSubtreeLikelihood() const;

        // this overloaded method is called during the tree-search
        virtual void computeSubtreeLikelihood(std::map<int, VVVdouble> *likelihoods,
                                              std::map<int, VVVdouble> *likelihoods_empty,
                                              const std::vector<int> &nodeList,
                                              std::map<int, bool> *ts_node__data_origin,
                                              tshlib::Utree &_utree__topology);


        virtual void _kernel_subtreelikelihood(int nodeID,
                                               VVVdouble *pxy__node,
                                               VVVdouble *_likelihoods__node,
                                               VVVdouble *_likelihoods_empty__node,
                                               std::vector<VVVdouble *> *lk_sons,
                                               std::vector<VVVdouble *> *lk_sons_empty);



//        virtual void _kernel_subtreelikelihood(int nodeID, VVVdouble *pxy__node,
//                                               VVVdouble *_likelihoods__node,
//                                               VVVdouble *_likelihoods_empty__node,
//                                               Vint *_sons__ids,
//                                               std::map<int, bool> *ts_node__data_origin);

        virtual void computeDownSubtreeDLikelihood(const Node *);

        virtual void computeDownSubtreeD2Likelihood(const Node *);

        void fireParameterChanged(const ParameterList &params);

        /**
         * @brief This method is mainly for debugging purpose.
         *
         * @param node The node at which likelihood values must be displayed.
         */
        virtual void displayLikelihood(const Node *node);

        /**
         * @brief This method sets DescCount (number of characters different from gap per column) value for all the nodes in the tree
         * @param sites SiteContainer of the aligned sites
         */
        virtual void setInsertionHistories(const SiteContainer &sites,
                                           const std::vector<int> &nodeList,
                                           std::map<int, std::vector<int>> *descCountData,
                                           std::map<int, std::vector<bool>> *setAData) const;

        /**
         * @brief This method method sets DescCount (number of characters different from gap per column) value for all the nodes in the tree during
         * the tree search
         * @param sites
         * @param nodeList
         * @param descCountData
         * @param setAData
         * @param _utree__topology
         */
        void setInsertionHistories(const SiteContainer &sites,
                                   const std::vector<int> &nodeList,
                                   std::map<int, std::vector<int>> *descCountData,
                                   std::map<int, std::vector<bool>> *setAData,
                                   std::map<int, bool> *ts_node__data_origin,
                                   tshlib::Utree &_utree__topology) const;
        /**
         * @brief This method sets the indicator for the number of evolutionary events (Insertions and Deletions) at each node of the topology
         * @param sites SiteContainer reference of the aligned sites
         */

        /**
         * @brief This method sets the setA (setA=1: possible insertion on that edge) value for all the nodes in the tree
         */
        /*virtual void setAllSetAData(const SiteContainer &sites) const;*/

        /**
         * @brief This method sets the iota value for all the nodes in the tree
         */
        virtual void setAllIotas();

        /**
         * @brief This method sets the beta value for all the nodes in the tree
         */
        virtual void setAllBetas();

        /**
         * @brief This method sets to 1 all the likelihood arrays recursively from a starting node
         * @param node The node at which the likelihood arrays must be reset
         */
        //virtual void resetNodeLikelihoodArrays(const Node *node);

        /**
         * @brief This method updates the likelihood arrays recursively from a starting node for a
         *        subtree
         * @param nodelist The postorder list of nodes at which the likelihood arrays must be updated
         */

        void setIndicatorFunction(const SiteContainer &sites) const;

        double computePhi(double lkEmptyColumn) const;

        void computeNu();

        void _printFV(Node *node, VVVdouble *likelihoodvector) const;

        void _printPrMatrix(Node *node, VVdouble *pr);

        std::vector<int> remapVirtualNodeLists(std::vector<int> &inputList) const;

        void _extendNodeListOnSetA(tshlib::VirtualNode *qnode, std::vector<int> &listNodes, unsigned long site, tshlib::Utree &candUTree) const;

        void _extendNodeListOnSetA(int qnodeID,
                                   unsigned long site,
                                   std::vector<int> &_node__list,
                                   std::map<int, std::vector<bool>> *_ts__setadata,
                                   tshlib::Utree &_utree__topology) const;


        double computeLikelihoodForASite(std::vector<int> &_node__list, size_t i) const;

        double computeLikelihoodForASite(size_t i,
                                         std::map<int, VVVdouble> *likelihoods,
                                         std::map<int, VVVdouble> *likelihoods_empty,
                                         std::map<int, std::vector<bool>> *ts_setadata,
                                         std::vector<int> &_node__list,
                                         std::map<int, bool> *ts_node__data_origin,
                                         tshlib::Utree &_utree__topology) const;

        double _kernel_likelihood_forasite(size_t i,
                                           int nodeID,
                                           std::vector<VVVdouble *> *lk_sons,
                                           std::vector<VVVdouble *> *lk_sons_empty,
                                           std::map<int, std::vector<bool>> *ts_setadata) const;


        double computeLikelihoodWholeAlignmentEmptyColumn() const;

        double computeLikelihoodWholeAlignmentEmptyColumn(std::map<int, VVVdouble> *ts_lkemptydata,
                                                          std::map<int, bool> *ts_node__data_origin,
                                                          tshlib::Utree &_utree__topology) const;

        double _kernel_likelihood_empty_forasite(int nodeID, std::vector<VVVdouble *> *lk_sons_empty) const;


        int countNonGapCharacterInSite(int siteID) const;

        void computeInDelDispersionOnTree(const SiteContainer &sites);

        SiteContainer *getSubAlignment(const SiteContainer &sites, std::vector<std::string> subsetSequences);

        void setNhNgOnNode(SiteContainer &sites, int nodeID);

        double getNodeAge(int nodeID);


    };
}
#endif //CASTOR_RHOMOGENEOUSTREELIKELIHOOD_PIP_HPP
