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
 * @file UnifiedTSHomogeneousTreeLikelihood.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 09 04 2018
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
#ifndef CASTOR_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_HPP
#define CASTOR_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_HPP

#include "UnifiedTSHSearchable.hpp"
#include "Utils.hpp"
#include "RHomogeneousTreeLikelihood_Generic.hpp"


namespace bpp {

    class UnifiedTSHomogeneousTreeLikelihood : public RHomogeneousTreeLikelihood_Generic, public UnifiedTSHSearchable {

    protected:
        //mutable DRASRTreeLikelihoodData *likelihoodData;
        mutable DRASRTreeLikelihoodData *likelihoodDataTest_;
        const SiteContainer* data_;
        mutable tshlib::Utree *utree_;
        mutable UtreeBppUtils::treemap treemap_;

        bool usePatterns_;


    public:
        UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                           TransitionModel *model,
                                           DiscreteDistribution *rDist,
                                           tshlib::Utree *utree_,
                                           UtreeBppUtils::treemap *treemap_,
                                           bool optNumericalDerivatives,
                                           std::map<std::string, std::string> &params,
                                           const std::string &suffix,
                                           bool checkRooted,
                                           bool verbose,
                                           bool usePatterns);


        UnifiedTSHomogeneousTreeLikelihood(const Tree &tree,
                                           const SiteContainer &data,
                                           TransitionModel *model,
                                           DiscreteDistribution *rDist,
                                           tshlib::Utree *utree_,
                                           UtreeBppUtils::treemap *treemap_,
                                           bool optNumericalDerivatives,
                                           std::map<std::string, std::string> &params,
                                           const std::string &suffix,
                                           bool checkRooted,
                                           bool verbose,
                                           bool usePatterns);

        ~UnifiedTSHomogeneousTreeLikelihood() ;

        UnifiedTSHomogeneousTreeLikelihood *clone() const override { return new UnifiedTSHomogeneousTreeLikelihood(*this); }

        const Tree &getTopology() const { return getTree(); }

        double getTopologyValue() const throw(Exception) { return getValue(); }

        tshlib::Utree *getUtreeTopology() {return utree_;}

        void init_(bool usePatterns);

        void fireTopologyChange(std::vector<int> nodeList, tshlib::Utree &_utree__topology);

        double updateLikelihoodOnTreeRearrangement(std::vector<int> &nodeList, tshlib::Utree &_utree__topology);

        double getLogLikelihoodOnTreeRearrangement() const;

        void topologyChangeSuccessful(std::vector<int> listNodes);

        void topologyCommitTree();

        std::vector<int> remapVirtualNodeLists(std::vector<int> &inputList) const;

        void updateLikelihoodArrays(Node *node, tshlib::Utree &_utree__topology);

        void resetLikelihoodsOnTopologyChangeSuccessful();

    };


}


#endif //CASTOR_UNIFIEDTSHOMOGENEOUSTREELIKELIHOOD_HPP
