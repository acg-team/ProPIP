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
 * @file RHomogeneousTreeLikelihood_Generic.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 23 04 2018
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
#ifndef CASTOR_RHOMOGENEOUSTREELIKELIHOOD_GENERIC_HPP
#define CASTOR_RHOMOGENEOUSTREELIKELIHOOD_GENERIC_HPP

#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Likelihood/DRASRTreeLikelihoodData.h>

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp {

    /**
     * @brief This class implement the 'traditional' way of computing likelihood for a tree.
     *
     * The substitution model is constant over the tree (homogeneous model).
     * A non uniform distribution of rates among the sites is allowed (ASRV models).</p>
     *
     * This class uses an instance of the DRASRTreeLikelihoodData for conditionnal likelihood storage.
     *
     * This class can also use a simple or recursive site compression.
     * In the simple case, computations for identical sites are not duplicated.
     * In the recursive case, computations for identical sub-sites (<i>site patterns </i>) are also not duplicated:
     * Following N. Galtier (personal communication ;-), we define a Pattern as a distinct site
     * in a sub-dataset corresponding to the dataset with sequences associated to a particular subtree.
     * The likelihood computation is the same for a given site, hence the idea is to save time from
     * performing many times the same coputation.
     * The network between all patterns is defined by the _patternLinks double map, initialized in the
     * initLikelihoodsWithPatterns() method. This initialisation takes more time than the classic
     * initTreeLikelihood one, where all likelihoods for a given site <i>i</i> are at the <i>i</i> coordinate
     * in the likelihood tensor, but is really faster when computing the likelihood (computeLikelihoods() method).
     * Hence, if you have to compute likelihood many times while holding the tree topology unchanged,
     * you should use patterns.
     * This decreases the likelihood computation time, but at a cost: some time is spent to establish the patterns
     * relationships. Whether to use or not patterns depends on what you actllay need:
     * - The more you compute likelihoods without changing the data or topology, the more patterns are interesting
     *   (this divides the cost of computing patterns by the number of computation performed).
     *   Patterns are hence usefull when you have a high number of computation to perform, while optimizing numerical
     *   parameters for instance).
     * - Patterns are more likely to occur whith small alphabet (nucleotides).
     */
    class RHomogeneousTreeLikelihood_Generic : public AbstractHomogeneousTreeLikelihood {
    protected:

        mutable DRASRTreeLikelihoodData *likelihoodData_;
        double minusLogLik_;


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
        RHomogeneousTreeLikelihood_Generic(
                const Tree &tree,
                TransitionModel *model,
                DiscreteDistribution *rDist,
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
        RHomogeneousTreeLikelihood_Generic(
                const Tree &tree,
                const SiteContainer &data,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                bool checkRooted = true,
                bool verbose = true,
                bool usePatterns = true)
        throw(Exception);

        RHomogeneousTreeLikelihood_Generic(const RHomogeneousTreeLikelihood_Generic &lik);

        RHomogeneousTreeLikelihood_Generic &operator=(const RHomogeneousTreeLikelihood_Generic &lik);

        virtual ~RHomogeneousTreeLikelihood_Generic();

        RHomogeneousTreeLikelihood_Generic *clone() const { return new RHomogeneousTreeLikelihood_Generic(*this); }

    protected:

        /**
         * @brief Method called by constructors.
         */
        void init_(bool usePatterns) throw(Exception);

    public:

        /**
         * @name The TreeLikelihood interface.
         *
         * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
         *
         * @{
         */
        void setData(const SiteContainer &sites) throw(Exception);

        double getLikelihood() const;

        double getLogLikelihood() const;

        double getLikelihoodForASite(size_t site) const;

        double getLogLikelihoodForASite(size_t site) const;
        /** @} */


        /**
         * @name The DiscreteRatesAcrossSites interface implementation:
         *
         * @{
         */
        double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

        double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;

        double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
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
        /** @} */

        /**
         * @name DerivableSecondOrder interface.
         *
         * @{
         */
        double getSecondOrderDerivative(const std::string &variable) const throw(Exception);

        double getSecondOrderDerivative(const std::string &variable1, const std::string &variable2) const throw(Exception) { return 0; } // Not implemented for now.
        /** @} */

    public:    // Specific methods:

        DRASRTreeLikelihoodData *getLikelihoodData() { return likelihoodData_; }

        const DRASRTreeLikelihoodData *getLikelihoodData() const { return likelihoodData_; }

        void computeTreeLikelihood();

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


    protected:

        /**
         * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
         *
         * @param node The root of the subtree.
         */
        virtual void computeSubtreeLikelihood(const Node *node); //Recursive method.
        virtual void computeDownSubtreeDLikelihood(const Node *);

        virtual void computeDownSubtreeD2Likelihood(const Node *);

        void fireParameterChanged(const ParameterList &params);

        /**
         * @brief This method is mainly for debugging purpose.
         *
         * @param node The node at which likelihood values must be displayed.
         */
        virtual void displayLikelihood(const Node *node);

        friend class RHomogeneousMixedTreeLikelihood;
    };


} //end of namespace bpp.


#endif //CASTOR_RHOMOGENEOUSTREELIKELIHOOD_GENERIC_HPP
