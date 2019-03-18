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
 * @file UnifiedDistanceEstimation.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 30 06 2018
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
#ifndef CASTOR_UNIFIEDDISTANCEESTIMATION_HPP
#define CASTOR_UNIFIEDDISTANCEESTIMATION_HPP

#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/DistanceMatrix.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Function/MetaOptimizer.h>
#include <Bpp/Numeric/Function/SimpleMultiDimensions.h>
#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>

namespace bpp {


/**
 * @brief This class is a simplified version of DRHomogeneousTreeLikelihood for 2-Trees.
 */
    class TwoTreeLikelihood_PIP :
            public AbstractDiscreteRatesAcrossSitesTreeLikelihood {

        enum class PairwiseSeqStates {
            cc, cg, gc, gg
        };
    private:
        SiteContainer *shrunkData_;
        std::vector<std::string> seqnames_;
        TransitionModel *model_;
        ParameterList brLenParameters_;

        mutable VVVdouble pxy_;

        mutable VVVdouble dpxy_;

        mutable VVVdouble d2pxy_;

        mutable std::vector<PairwiseSeqStates> setA_;

        /**
         * @brief As previous, but for the global container.
         *
         * The size of this vector is equal to the number of sites in the container,
         * each element corresponds to a site in the container and points to the
         * corresponding column in the likelihood array of the root node.
         * If the container contains no repeated site, there will be a strict
         * equivalence between each site and the likelihood array of the root node.
         * However, if this is not the case, some pointers may point toward the same
         * element in the likelihood array.
         */
        std::vector<size_t> rootPatternLinks_;

        /**
         * @brief The frequency of each site.
         */
        std::vector<unsigned int> rootWeights_;

        //some values we'll need:
        size_t nbSites_,         //the number of sites in the container
                nbClasses_,       //the number of rate classes
                nbStates_,        //the number of states in the alphabet
                nbDistinctSites_; //the number of distinct sites in the container

        mutable VVVdouble rootLikelihoods_;
        mutable VVdouble rootLikelihoodsS_;
        mutable Vdouble rootLikelihoodsSR_;
        mutable Vdouble dLikelihoods_;
        mutable Vdouble d2Likelihoods_;
        mutable VVdouble leafLikelihoods1_, leafLikelihoods2_;

        double minimumBrLen_;
        Constraint *brLenConstraint_;
        double brLen_;

    public:
        TwoTreeLikelihood_PIP(
                const std::string &seq1, const std::string &seq2,
                const SiteContainer &data,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                bool verbose) throw(Exception);

        TwoTreeLikelihood_PIP(const TwoTreeLikelihood_PIP &lik);

        TwoTreeLikelihood_PIP &operator=(const TwoTreeLikelihood_PIP &lik);

        TwoTreeLikelihood_PIP *clone() const { return new TwoTreeLikelihood_PIP(*this); }

        virtual ~TwoTreeLikelihood_PIP();

    public:

        /**
         * @name The TreeLikelihood interface.
         *
         * Other methods are implemented in the AbstractTreeLikelihood class.
         *
         * @{
         */
        size_t getNumberOfStates() const { return model_->getNumberOfStates(); }

        const std::vector<int> &getAlphabetStates() const { return model_->getAlphabetStates(); }

        int getAlphabetStateAsInt(size_t i) const { return model_->getAlphabetStateAsInt(i); }

        std::string getAlphabetStateAsChar(size_t i) const { return model_->getAlphabetStateAsChar(i); }

        TreeLikelihoodData *getLikelihoodData() throw(NotImplementedException) {
            throw NotImplementedException("TwoTreeLikelihood_PIP::getLikelihoodData.");
        }

        const TreeLikelihoodData *getLikelihoodData() const throw(NotImplementedException) {
            throw NotImplementedException("TwoTreeLikelihood_PIP::getLikelihoodData.");
        }

        double getLikelihood() const;

        double getLogLikelihood() const;

        double getLikelihoodForASite(size_t site) const;

        double getLogLikelihoodForASite(size_t site) const;

        ParameterList getBranchLengthsParameters() const;

        ParameterList getSubstitutionModelParameters() const;

        TransitionModel *getModelForSite(int nodeId, size_t siteIndex) { return model_; }

        const TransitionModel *getModelForSite(int nodeId, size_t siteIndex) const { return model_; }

        const std::vector<double> &getRootFrequencies(size_t siteIndex) const { return model_->getFrequencies(); }

        size_t getSiteIndex(size_t site) const throw(IndexOutOfBoundsException) { return rootPatternLinks_[site]; }

        /**
         * @brief This method is not applicable for this object.
         */
        VVVdouble getTransitionProbabilitiesPerRateClass(int nodeId, size_t siteIndex) const { return pxy_; }

        void setData(const SiteContainer &sites) throw(Exception) {}

        void initialize() throw(Exception);
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

        const TransitionModel *getModel() const { return model_; }

        /**
         * @brief Get the substitution model used for the computation.
         *
         * @return A pointer toward the substitution model of this instance.
         */
        TransitionModel *getModel() { return model_; }

        ConstBranchModelIterator *getNewBranchModelIterator(int nodeId) const throw(NotImplementedException) {
            throw NotImplementedException(
                    "TwoTreeLikelihood_PIP::getNewBranchSiteModelIterator. This class does not (yet) provide support for partition models.");
        }

        ConstSiteModelIterator *getNewSiteModelIterator(size_t siteIndex) const throw(NotImplementedException) {
            throw NotImplementedException(
                    "TwoTreeLikelihood_PIP::getNewSiteModelIterator. This class is for inner use only and does not provide site model iterators.");
        }


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

        /**
         * @name DerivableFirstOrder interface.
         *
         * @{
         */
        double getFirstOrderDerivative(const std::string &variable) const throw(Exception);
        /** @{ */

        /**
         * @name DerivableSecondOrder interface.
         *
         * @{
         */
        double getSecondOrderDerivative(const std::string &variable) const throw(Exception);

        double getSecondOrderDerivative(const std::string &variable1,
                                        const std::string &variable2) const throw(Exception) { return 0; } // Not implemented for now.
        /** @} */

        virtual void initBranchLengthsParameters();

        virtual void setMinimumBranchLength(double minimum) {
            minimumBrLen_ = minimum;
            if (brLenConstraint_) delete brLenConstraint_;
            brLenConstraint_ = new IntervalConstraint(1, minimumBrLen_, true);
            initBranchLengthsParameters();
        }

        virtual double getMinimumBranchLength() const { return minimumBrLen_; }

    protected:

        /**
         * @brief This method initializes the leaves according to a sequence container.
         *
         * Here the container shrunkData_ is used.
         * Likelihood is set to 1 for the state corresponding to the sequence site,
         * otherwise it is set to 0.
         *
         * The two likelihood arrays are initialized according to alphabet
         * size and sequences length, and filled with 1.
         *
         * NB: This method is recursive.
         *
         * @param sequences The sequence container to use.
         */
        virtual void initTreeLikelihoods(const SequenceContainer &sequences) throw(Exception);

        void fireParameterChanged(const ParameterList &params);

        virtual void computeTreeLikelihood();

        virtual void computeTreeDLikelihood();

        virtual void computeTreeD2Likelihood();

        virtual double computeEmptyColumnLikelihood() const;

        /**
         * @brief This builds the <i>parameters</i> list from all parametrizable objects,
         * <i>i.e.</i> substitution model, rate distribution and tree.
         */
        virtual void initParameters();

        /**
         * @brief All parameters are stores in a parameter list.
         *
         * This function apply these parameters to the substitution model,
         * to the rate distribution and to the branch lengths.
         */
        virtual void applyParameters() throw(Exception);

    };

/**
 * @brief Estimate a distance matrix from sequence data, according to a given model.
 *
 * By default, the parameters of the model are fixed to there given values.
 * It is possible to estimate one or several parameters by setting them with the
 * setAdditionalParameters() method.
 * Parameters will be estimated separately for each pair of sequence.
 *
 * For now it is not possible to retrieve estimated values.
 * You'll have to specify a 'profiler' to the optimizer and then look at the file
 * if you want to do so.
 */
    class UnifiedDistanceEstimation :
            public virtual Clonable {
    private:
        std::unique_ptr<TransitionModel> model_;
        std::unique_ptr<DiscreteDistribution> rateDist_;
        const SiteContainer *sites_;
        DistanceMatrix *dist_;
        Optimizer *optimizer_;
        MetaOptimizer *defaultOptimizer_;
        size_t verbose_;
        ParameterList parameters_;

    public:

        /**
         * @brief Create a new DistanceEstimation object according to a given substitution model and a rate distribution.
         *
         * This instance will own the model and distribution, and will take car of their recopy and destruction.
         *
         * @param model    The substitution model to use.
         * @param rateDist The discrete rate distribution to use.
         * @param verbose  The verbose level:
         *  - 0=Off,
         *  - 1=one * by row computation
         *  - 2=one * by row computation and one . by column computation
         *  - 3=2 + optimization verbose enabled
         *  - 4=3 + likelihood object verbose enabled
         */
        UnifiedDistanceEstimation(
                TransitionModel *model,
                DiscreteDistribution *rateDist,
                size_t verbose = 1) :
                model_(model),
                rateDist_(rateDist),
                sites_(0),
                dist_(0),
                optimizer_(0),
                defaultOptimizer_(0),
                verbose_(verbose),
                parameters_() {
            init_();
        }

        /**
         * @brief Create a new DistanceEstimation object and compute distances
         * according to a given substitution model and a rate distribution.
         *
         * This instance will own the model and distribution, and will take car of their recopy and destruction.
         *
         * @param model    The substitution model to use.
         * @param rateDist The discrete rate distribution to use.
         * @param sites    The sequence data.
         * @param verbose  The verbose level:
         *  - 0=Off,
         *  - 1=one * by row computation
         *  - 2=one * by row computation and one . by column computation
         *  - 3=2 + optimization verbose enabled
         *  - 4=3 + likelihood object verbose enabled
         *  @param computeMat if true the computeMatrix() method is called.
         */
        UnifiedDistanceEstimation(
                TransitionModel *model,
                DiscreteDistribution *rateDist,
                const SiteContainer *sites,
                size_t verbose = 1,
                bool computeMat = true) :
                model_(model),
                rateDist_(rateDist),
                sites_(sites),
                dist_(0),
                optimizer_(0),
                defaultOptimizer_(0),
                verbose_(verbose),
                parameters_() {
            init_();
            if (computeMat) computeMatrix();
        }

        /**
         * @brief Copy constructor.
         *
         * Only the distance matrix is hard-copied, if there is one.
         *
         * @param distanceEstimation The object to copy.
         */
        UnifiedDistanceEstimation(const UnifiedDistanceEstimation &distanceEstimation) :
                model_(distanceEstimation.model_->clone()),
                rateDist_(distanceEstimation.rateDist_->clone()),
                sites_(distanceEstimation.sites_),
                dist_(0),
                optimizer_(dynamic_cast<Optimizer *>(distanceEstimation.optimizer_->clone())),
                defaultOptimizer_(dynamic_cast<MetaOptimizer *>(distanceEstimation.defaultOptimizer_->clone())),
                verbose_(distanceEstimation.verbose_),
                parameters_(distanceEstimation.parameters_) {
            if (distanceEstimation.dist_ != 0)
                dist_ = new DistanceMatrix(*distanceEstimation.dist_);
            else
                dist_ = 0;
        }

        /**
         * @brief Assigment operator.
         *
         * Only the distance matrix is hard-copied, if there is one.
         *
         * @param distanceEstimation The object to copy.
         * @return A reference toward this object.
         */
        UnifiedDistanceEstimation &operator=(const UnifiedDistanceEstimation &distanceEstimation) {
            model_.reset(distanceEstimation.model_->clone());
            rateDist_.reset(distanceEstimation.rateDist_->clone());
            sites_ = distanceEstimation.sites_;
            if (distanceEstimation.dist_ != 0)
                dist_ = new DistanceMatrix(*distanceEstimation.dist_);
            else
                dist_ = 0;
            optimizer_ = dynamic_cast<Optimizer *>(distanceEstimation.optimizer_->clone());
            // _defaultOptimizer has already been initialized since the default constructor has been called.
            verbose_ = distanceEstimation.verbose_;
            parameters_ = distanceEstimation.parameters_;
            return *this;
        }

        virtual ~UnifiedDistanceEstimation() {
            if (dist_) delete dist_;
            delete defaultOptimizer_;
            delete optimizer_;
        }

        UnifiedDistanceEstimation *clone() const { return new UnifiedDistanceEstimation(*this); }

    private:
        void init_() {


            MetaOptimizerInfos *desc = new MetaOptimizerInfos();
            std::vector<std::string> name;
            name.push_back("BrLen");

            if (model_->getName().find("PIP") != string::npos) {
                desc->addOptimizer("Branch length", new SimpleMultiDimensions(0), name, 0, MetaOptimizerInfos::IT_TYPE_STEP);
            } else {
                desc->addOptimizer("Branch length", new PseudoNewtonOptimizer(0), name, 2, MetaOptimizerInfos::IT_TYPE_FULL);
            }

            ParameterList tmp = model_->getParameters();
            tmp.addParameters(rateDist_->getParameters());
            desc->addOptimizer("substitution model and rate distribution", new SimpleMultiDimensions(0), tmp.getParameterNames(), 0,
                               MetaOptimizerInfos::IT_TYPE_STEP);
            defaultOptimizer_ = new MetaOptimizer(0, desc);
            defaultOptimizer_->setMessageHandler(0);
            defaultOptimizer_->setProfiler(0);
            defaultOptimizer_->getStopCondition()->setTolerance(0.0001);
            optimizer_ = dynamic_cast<Optimizer *>(defaultOptimizer_->clone());
        }

    public:

        /**
         * @brief Perform the distance computation.
         *
         * Result can be called by the getMatrix() method.
         *
         * @throw NullPointerException if at least one of the model,
         * rate distribution or data are not initialized.
         */
        void computeMatrix() throw(NullPointerException);

        /**
         * @brief Get the distance matrix.
         *
         * @return A pointer toward the computed distance matrix.
         */
        DistanceMatrix *getMatrix() const { return dist_ == 0 ? 0 : new DistanceMatrix(*dist_); }

        bool hasModel() const { return model_.get(); }

        const TransitionModel &getModel() const throw(Exception) {
            if (hasModel())
                return *model_;
            else
                throw Exception("DistanceEstimation::getSubstitutionModel(). No model assciated to this instance.");
        }

        void resetSubstitutionModel(TransitionModel *model = 0) { model_.reset(model); }

        bool hasRateDistribution() const { return rateDist_.get(); }

        const DiscreteDistribution &getRateDistribution() const throw(Exception) {
            if (hasRateDistribution())
                return *rateDist_;
            else
                throw Exception("DistanceEstimation::getRateDistribution(). No rate distribution assciated to this instance.");
        }

        void resetRateDistribution(DiscreteDistribution *rateDist = 0) { rateDist_.reset(rateDist); }

        void setData(const SiteContainer *sites) { sites_ = sites; }

        const SiteContainer *getData() const { return sites_; }

        void resetData() { sites_ = 0; }

        void setOptimizer(const Optimizer *optimizer) {
            if (optimizer_) delete optimizer_;
            optimizer_ = dynamic_cast<Optimizer *>(optimizer->clone());
        }

        const Optimizer *getOptimizer() const { return optimizer_; }

        Optimizer *getOptimizer() { return optimizer_; }

        void resetOptimizer() { optimizer_ = dynamic_cast<Optimizer *>(defaultOptimizer_->clone()); }

        /**
         * @brief Specify a list of parameters to be estimated.
         *
         * Parameters will be estimated separately for each distance.
         *
         * @param parameters A list of parameters to estimate.
         */
        void setAdditionalParameters(const ParameterList &parameters) {
            parameters_ = parameters;
        }

        /**
         * @brief Reset all additional parameters.
         */
        void resetAdditionalParameters() {
            parameters_.reset();
        }

        /**
         * @param verbose Verbose level.
         */
        void setVerbose(size_t verbose) { verbose_ = verbose; }

        /**
         * @return Verbose level.
         */
        size_t getVerbose() const { return verbose_; }
    };
}

#endif //CASTOR_UNIFIEDDISTANCEESTIMATION_HPP
