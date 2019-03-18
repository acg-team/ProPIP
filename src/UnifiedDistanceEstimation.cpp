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
 * @file UnifiedDistanceEstimation.cpp
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
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Seq/SymbolListTools.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include "UnifiedDistanceEstimation.hpp"
#include "Utils.hpp"

using namespace bpp;


TwoTreeLikelihood_PIP::TwoTreeLikelihood_PIP(
        const std::string &seq1, const std::string &seq2,
        const SiteContainer &data,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        bool verbose) throw(Exception) :
        AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose),
        shrunkData_(0), seqnames_(2), model_(model), brLenParameters_(), pxy_(), dpxy_(), d2pxy_(),
        rootPatternLinks_(), rootWeights_(), nbSites_(0), nbClasses_(0), nbStates_(0), nbDistinctSites_(0),
        rootLikelihoods_(), rootLikelihoodsS_(), rootLikelihoodsSR_(), dLikelihoods_(), d2Likelihoods_(),
        leafLikelihoods1_(), leafLikelihoods2_(),
        minimumBrLen_(0.000001), brLenConstraint_(0), brLen_(0) {
    seqnames_[0] = seq1;
    seqnames_[1] = seq2;
    data_ = PatternTools::getSequenceSubset(data, seqnames_);
    if (data_->getAlphabet()->getAlphabetType()
        != model_->getAlphabet()->getAlphabetType())
        throw AlphabetMismatchException("TwoTreeTreeLikelihood::TwoTreeTreeLikelihood. Data and model must have the same alphabet type.",
                                        data_->getAlphabet(),
                                        model_->getAlphabet());

    nbSites_ = data_->getNumberOfSites();
    nbClasses_ = rateDistribution_->getNumberOfCategories();
    nbStates_ = model_->getNumberOfStates();
    if (verbose)
        ApplicationTools::displayMessage("Double-Recursive Homogeneous Tree Likelihood");

    // Initialize root patterns:
    SitePatterns pattern(data_);
    shrunkData_ = pattern.getSites();
    rootWeights_ = pattern.getWeights();
    rootPatternLinks_ = pattern.getIndices();
    nbDistinctSites_ = shrunkData_->getNumberOfSites();
    if (verbose)
        ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));

    // Init _likelihoods:
    if (verbose) ApplicationTools::displayTask("Init likelihoods arrays recursively");
    // Clone data for more efficiency on sequences access:
    const SiteContainer *sequences = new AlignedSequenceContainer(*shrunkData_);
    initTreeLikelihoods(*sequences);
    delete sequences;

    brLen_ = minimumBrLen_;
    brLenConstraint_ = new IntervalConstraint(1, minimumBrLen_, true);

    if (verbose) ApplicationTools::displayTaskDone();
}


TwoTreeLikelihood_PIP::TwoTreeLikelihood_PIP(const TwoTreeLikelihood_PIP &lik) :
        AbstractDiscreteRatesAcrossSitesTreeLikelihood(lik),
        shrunkData_(dynamic_cast<SiteContainer *>(lik.shrunkData_->clone())),
        seqnames_(lik.seqnames_),
        model_(lik.model_),
        brLenParameters_(lik.brLenParameters_),
        pxy_(lik.pxy_),
        dpxy_(lik.dpxy_),
        d2pxy_(lik.d2pxy_),
        rootPatternLinks_(lik.rootPatternLinks_),
        rootWeights_(lik.rootWeights_),
        nbSites_(lik.nbSites_),
        nbClasses_(lik.nbClasses_),
        nbStates_(lik.nbStates_),
        nbDistinctSites_(lik.nbDistinctSites_),
        rootLikelihoods_(lik.rootLikelihoods_),
        rootLikelihoodsS_(lik.rootLikelihoodsS_),
        rootLikelihoodsSR_(lik.rootLikelihoodsSR_),
        dLikelihoods_(lik.dLikelihoods_),
        d2Likelihoods_(lik.d2Likelihoods_),
        leafLikelihoods1_(lik.leafLikelihoods1_),
        leafLikelihoods2_(lik.leafLikelihoods2_),
        minimumBrLen_(lik.minimumBrLen_),
        brLenConstraint_(dynamic_cast<Constraint *>(lik.brLenConstraint_->clone())),
        brLen_(lik.brLen_) {}


TwoTreeLikelihood_PIP &TwoTreeLikelihood_PIP::operator=(const TwoTreeLikelihood_PIP &lik) {
    AbstractDiscreteRatesAcrossSitesTreeLikelihood::operator=(lik);
    shrunkData_ = dynamic_cast<SiteContainer *>(lik.shrunkData_->clone());
    seqnames_ = lik.seqnames_;
    model_ = lik.model_;
    brLenParameters_ = lik.brLenParameters_;
    pxy_ = lik.pxy_;
    dpxy_ = lik.dpxy_;
    d2pxy_ = lik.d2pxy_;
    rootPatternLinks_ = lik.rootPatternLinks_;
    rootWeights_ = lik.rootWeights_;
    nbSites_ = lik.nbSites_;
    nbClasses_ = lik.nbClasses_;
    nbStates_ = lik.nbStates_;
    nbDistinctSites_ = lik.nbDistinctSites_;
    rootLikelihoods_ = lik.rootLikelihoods_;
    rootLikelihoodsS_ = lik.rootLikelihoodsS_;
    rootLikelihoodsSR_ = lik.rootLikelihoodsSR_;
    dLikelihoods_ = lik.dLikelihoods_;
    d2Likelihoods_ = lik.d2Likelihoods_;
    leafLikelihoods1_ = lik.leafLikelihoods1_;
    leafLikelihoods2_ = lik.leafLikelihoods2_;
    minimumBrLen_ = lik.minimumBrLen_;
    brLenConstraint_ = dynamic_cast<Constraint *>(brLenConstraint_->clone());
    brLen_ = lik.brLen_;
    return *this;
}


TwoTreeLikelihood_PIP::~TwoTreeLikelihood_PIP() {
    delete shrunkData_;
    if (brLenConstraint_) delete brLenConstraint_;
}


void TwoTreeLikelihood_PIP::initialize() throw(Exception) {
    initParameters();
    initialized_ = true;
    fireParameterChanged(getParameters());
}


ParameterList TwoTreeLikelihood_PIP::getBranchLengthsParameters() const {
    if (!initialized_) throw Exception("TwoTreeLikelihood_PIP::getBranchLengthsParameters(). Object is not initialized.");
    return brLenParameters_.getCommonParametersWith(getParameters());
}


ParameterList TwoTreeLikelihood_PIP::getSubstitutionModelParameters() const {
    if (!initialized_) throw Exception("TwoTreeLikelihood_PIP::getSubstitutionModelParameters(). Object is not initialized.");
    return model_->getParameters().getCommonParametersWith(getParameters());
}


double TwoTreeLikelihood_PIP::getLikelihood() const {
    double l = 1.;
    for (size_t i = 0; i < nbDistinctSites_; i++) {
        l *= std::pow(rootLikelihoodsSR_[i], (int) rootWeights_[i]);
    }
    return l;
}


double TwoTreeLikelihood_PIP::getLogLikelihood() const {

    double logemptycol = computeEmptyColumnLikelihood();
    double ll = 0;
    for (size_t i = 0; i < nbDistinctSites_; i++) {
        ll += rootWeights_[i] * log(rootLikelihoodsSR_[i]);

    }
    return ll + logemptycol;
}

double TwoTreeLikelihood_PIP::computeEmptyColumnLikelihood() const {

    double lambda = model_->getParameter("lambda").getValue();
    double mu = model_->getParameter("mu").getValue();
    double nu = lambda * (brLen_ + 1 / mu);


    //compute empty column likelihood
    double rootBeta = 1;
    double rootIota = (1 / mu) / (brLen_ + 1 / mu);
    double branchBeta = (1 - exp(-mu * brLen_)) / (mu * brLen_);
    double branchIota = brLen_ / (brLen_ + 1 / mu);

    // compute root lk
    Vdouble fr = model_->getFrequencies();

    ColMatrix<double> emptycolumn(nbStates_, 1);
    ColMatrix<double> empty_l1(nbStates_, 1);
    ColMatrix<double> freqs(nbStates_, 1);
    ColMatrix<double> pxy(nbStates_, nbStates_);

    emptycolumn.operator()(model_->getAlphabet()->getGapCharacterCode(), 0) = 1;

    for (int i = 0; i < nbStates_; i++) {
        for (int o = 0; o < nbStates_; o++) {
            pxy.operator()(i, o) = pxy_[0][i][o];
        }
        freqs.operator()(i, 0) = fr.at(i);
    }


    MatrixTools::mult(pxy, emptycolumn, empty_l1);
    double llempty_root = MatrixBppUtils::dotProd(freqs, empty_l1);

    double lkemptyRoot = rootBeta * rootIota * llempty_root;
    double lkemptyBranch = branchIota * (1 - branchBeta + branchBeta * llempty_root);

    double llempty = lkemptyRoot + lkemptyBranch;

    //compute phi
    double logphy;
    double log_factorial_m;
    int m = (int) nbSites_;

    log_factorial_m = 0;
    for (int i = 1; i <= m; i++) {
        log_factorial_m += log(i);
    }

    logphy = -log_factorial_m + m * log(nu) + (nu * (llempty - 1));

    return logphy;
}

double TwoTreeLikelihood_PIP::getLikelihoodForASite(size_t site) const {
    return rootLikelihoodsSR_[rootPatternLinks_[site]];
}


double TwoTreeLikelihood_PIP::getLogLikelihoodForASite(size_t site) const {
    return log(rootLikelihoodsSR_[rootPatternLinks_[site]]);
}


double TwoTreeLikelihood_PIP::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {
    return rootLikelihoodsS_[rootPatternLinks_[site]][rateClass];
}


double TwoTreeLikelihood_PIP::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const {
    return log(rootLikelihoodsS_[rootPatternLinks_[site]][rateClass]);
}


double TwoTreeLikelihood_PIP::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const {
    return rootLikelihoods_[rootPatternLinks_[site]][rateClass][static_cast<size_t>(state)];
}


double TwoTreeLikelihood_PIP::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const {
    return log(rootLikelihoods_[rootPatternLinks_[site]][rateClass][static_cast<size_t>(state)]);
}


void TwoTreeLikelihood_PIP::initParameters() {
    // Reset parameters:
    resetParameters_();

    // Branch lengths:
    initBranchLengthsParameters();
    addParameters_(brLenParameters_);

    // Substitution model:
    addParameters_(model_->getIndependentParameters());

    // Rate distribution:
    addParameters_(rateDistribution_->getIndependentParameters());
}


void TwoTreeLikelihood_PIP::applyParameters() throw(Exception) {
    // Apply branch length:
    brLen_ = getParameterValue("BrLen");
    // Apply substitution model parameters:
    model_->matchParametersValues(getParameters());
    // Apply rate distribution parameters:
    rateDistribution_->matchParametersValues(getParameters());
}


void TwoTreeLikelihood_PIP::initBranchLengthsParameters() {
    if (brLen_ < minimumBrLen_) {
        ApplicationTools::displayWarning(
                "Branch length is too small: " + TextTools::toString(brLen_) + ". Value is set to " + TextTools::toString(minimumBrLen_));
        brLen_ = minimumBrLen_;
    }
    brLenParameters_.reset();
    brLenParameters_.addParameter(Parameter("BrLen", brLen_, brLenConstraint_));
}


void TwoTreeLikelihood_PIP::setParameters(const ParameterList &parameters)
throw(ParameterNotFoundException, ConstraintException) {
    setParametersValues(parameters);
}


void TwoTreeLikelihood_PIP::fireParameterChanged(const ParameterList &params) {
    applyParameters();

    // For now we ignore the parameter that changed and we recompute all arrays...

    // Computes all pxy and pyx once for all:
    pxy_.resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++) {
        VVdouble *pxy_c = &pxy_[c];
        pxy_c->resize(nbStates_);
        RowMatrix<double> Q = model_->getPij_t(brLen_ * rateDistribution_->getCategory(c));
        for (size_t x = 0; x < nbStates_; x++) {
            Vdouble *pxy_c_x = &(*pxy_c)[x];
            pxy_c_x->resize(nbStates_);
            for (size_t y = 0; y < nbStates_; y++) {
                (*pxy_c_x)[y] = Q(x, y);
            }
        }
    }
/*
    if (computeFirstOrderDerivatives_) {
        // Computes all dpxy/dt once for all:
        dpxy_.resize(nbClasses_);
        for (size_t c = 0; c < nbClasses_; c++) {
            VVdouble *dpxy_c = &dpxy_[c];
            dpxy_c->resize(nbStates_);
            double rc = rateDistribution_->getCategory(c);
            RowMatrix<double> dQ = model_->getdPij_dt(brLen_ * rc);
            for (size_t x = 0; x < nbStates_; x++) {
                Vdouble *dpxy_c_x = &(*dpxy_c)[x];
                dpxy_c_x->resize(nbStates_);
                for (size_t y = 0; y < nbStates_; y++) {
                    (*dpxy_c_x)[y] = rc * dQ(x, y);
                }
            }
        }
    }

    if (computeSecondOrderDerivatives_) {
        // Computes all d2pxy/dt2 once for all:
        d2pxy_.resize(nbClasses_);
        for (size_t c = 0; c < nbClasses_; c++) {
            VVdouble *d2pxy_c = &d2pxy_[c];
            d2pxy_c->resize(nbStates_);
            double rc = rateDistribution_->getCategory(c);
            RowMatrix<double> d2Q = model_->getd2Pij_dt2(brLen_ * rc);
            for (size_t x = 0; x < nbStates_; x++) {
                Vdouble *d2pxy_c_x = &(*d2pxy_c)[x];
                d2pxy_c_x->resize(nbStates_);
                for (size_t y = 0; y < nbStates_; y++) {
                    (*d2pxy_c_x)[y] = rc * rc * d2Q(x, y);
                }
            }
        }
    }*/

    computeTreeLikelihood();

    /*
    if (computeFirstOrderDerivatives_) {
        computeTreeDLikelihood();
    }
    if (computeSecondOrderDerivatives_) {
        computeTreeD2Likelihood();
    }
    */
}


double TwoTreeLikelihood_PIP::getValue() const
throw(Exception) {
    return -getLogLikelihood();
}


void TwoTreeLikelihood_PIP::initTreeLikelihoods(const SequenceContainer &sequences) throw(Exception) {
    const Sequence *seq1 = &sequences.getSequence(seqnames_[0]);
    const Sequence *seq2 = &sequences.getSequence(seqnames_[1]);

    const int gapChar = model_->getAlphabet()->getGapCharacterCode();

    leafLikelihoods1_.resize(nbDistinctSites_);
    leafLikelihoods2_.resize(nbDistinctSites_);
    setA_.resize(nbDistinctSites_);

    for (size_t i = 0; i < nbDistinctSites_; i++) {
        Vdouble *leafLikelihoods1_i = &leafLikelihoods1_[i];
        Vdouble *leafLikelihoods2_i = &leafLikelihoods2_[i];
        leafLikelihoods1_i->resize(nbStates_);
        leafLikelihoods2_i->resize(nbStates_);
        int state1 = seq1->getValue(i);
        int state2 = seq2->getValue(i);

        // Set state of the sequence
        if (state1 == gapChar && state2 != gapChar) {
            // GC state
            setA_[i] = PairwiseSeqStates::gc;
        } else if (state1 != gapChar && state2 == gapChar) {
            // CG state
            setA_[i] = PairwiseSeqStates::cg;
        } else if (state1 != gapChar && state2 != gapChar) {
            // CC state
            setA_[i] = PairwiseSeqStates::cc;
        } else {
            // GG state
            setA_[i] = PairwiseSeqStates::gg;
        }

        for (size_t s = 0; s < nbStates_; s++) {
            // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
            // otherwise value set to 0:
            try {
                (*leafLikelihoods1_i)[s] = model_->getInitValue(s, state1);
                (*leafLikelihoods2_i)[s] = model_->getInitValue(s, state2);
            }
            catch (SequenceNotFoundException &snfe) {
                throw SequenceNotFoundException("TwoTreeLikelihood_PIP::initTreelikelihoods. Leaf name in tree not found in site container: ",
                                                snfe.getSequenceId());
            }
        }
    }

    // Initialize likelihood vector:
    rootLikelihoods_.resize(nbDistinctSites_);
    rootLikelihoodsS_.resize(nbDistinctSites_);
    rootLikelihoodsSR_.resize(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++) {
        VVdouble *rootLikelihoods_i = &rootLikelihoods_[i];
        Vdouble *rootLikelihoodsS_i = &rootLikelihoodsS_[i];
        rootLikelihoods_i->resize(nbClasses_);
        rootLikelihoodsS_i->resize(nbClasses_);
        for (size_t c = 0; c < nbClasses_; c++) {
            Vdouble *rootLikelihoods_i_c = &(*rootLikelihoods_i)[c];
            rootLikelihoods_i_c->resize(nbStates_);
            for (size_t s = 0; s < nbStates_; s++) {
                (*rootLikelihoods_i_c)[s] = 1.; // All likelihoods are initialized to 1.
            }
        }
    }

    // Initialize d and d2 likelihoods:
    dLikelihoods_.resize(nbDistinctSites_);
    d2Likelihoods_.resize(nbDistinctSites_);
}


void TwoTreeLikelihood_PIP::computeTreeLikelihood() {
    // get parameters PIP
    double mu = model_->getParameter("mu").getValue();
    double tau = brLen_;

    // compute PIP composite parameters
    double betaRoot = 1;
    double iotaRoot = (1 / mu) / (tau + 1 / mu);
    double betaBranch = (1 - exp(-mu * brLen_)) / (mu * brLen_);
    double iotaBranch = brLen_ / (tau + 1 / mu);

    for (size_t i = 0; i < nbDistinctSites_; i++) {
        VVdouble *rootLikelihoods_i = &rootLikelihoods_[i];
        Vdouble *leafLikelihoods1_i = &leafLikelihoods1_[i];
        Vdouble *leafLikelihoods2_i = &leafLikelihoods2_[i];
        for (size_t c = 0; c < nbClasses_; c++) {
            Vdouble *rootLikelihoods_i_c = &(*rootLikelihoods_i)[c];
            VVdouble *pxy_c = &pxy_[c];
            for (size_t x = 0; x < nbStates_; x++) {
                Vdouble *pxy_c_x = &(*pxy_c)[x];
                double l = 0;
                double l1 = (*leafLikelihoods1_i)[x];
                for (size_t y = 0; y < nbStates_; y++) {
                    double l2 = (*leafLikelihoods2_i)[y];
                    l += l1 * l2 * (*pxy_c_x)[y];
                }
                (*rootLikelihoods_i_c)[x] = l;
            }
        }
    }

    Vdouble fr = model_->getFrequencies();
    Vdouble p = rateDistribution_->getProbabilities();
    for (size_t i = 0; i < nbDistinctSites_; i++) {
        // For each site in the sequence,
        VVdouble *rootLikelihoods_i = &rootLikelihoods_[i];
        Vdouble *rootLikelihoodsS_i = &rootLikelihoodsS_[i];
        rootLikelihoodsSR_[i] = 0;
        for (size_t c = 0; c < nbClasses_; c++) {
            (*rootLikelihoodsS_i)[c] = 0;
            // For each rate classe,
            Vdouble *rootLikelihoods_i_c = &(*rootLikelihoods_i)[c];
            for (size_t x = 0; x < nbStates_; x++) {
                // For each initial state,
                (*rootLikelihoodsS_i)[c] += fr[x] * (*rootLikelihoods_i_c)[x];
            }
            rootLikelihoodsSR_[i] += p[c] * (*rootLikelihoodsS_i)[c];
        }

        if (setA_[i] == PairwiseSeqStates::cc || setA_[i] == PairwiseSeqStates::cg) {

            rootLikelihoodsSR_[i] = betaRoot * iotaRoot * rootLikelihoodsSR_[i];

        } else if (setA_[i] == PairwiseSeqStates::gc) {

            Vdouble *leafLikelihoods2_i = &leafLikelihoods2_[i];

            double leaflk = 0;
            for (size_t x = 0; x < nbStates_; x++) {
                leaflk += fr[x] * (*leafLikelihoods2_i)[x];
            }
            //leaflk = betaBranch * iotaBranch * leaflk;

            //double rootlk = 0;
            //rootlk = betaRoot * iotaRoot * rootLikelihoodsSR_[i];
            //rootLikelihoodsSR_[i] = betaRoot * iotaRoot * rootLikelihoodsSR_[i];

            rootLikelihoodsSR_[i] = betaBranch * iotaBranch * leaflk;

        } else {
            // this is the case for -|-
            rootLikelihoodsSR_[i] = iotaBranch * (1 - betaBranch);
        }

        //else if (setA_[i] == PairwiseSeqStates::cg) {

        //    Vdouble *leafLikelihoods1_i = &leafLikelihoods1_[i];
        //    double leaflk = 0;
        //    double rootlk = 0;
        //    double tmp = rootLikelihoodsSR_[i];
        //for (size_t x = 0; x < nbStates_; x++) {
        //    leaflk += fr[x] * (*leafLikelihoods1_i)[x];
        //  }
        //leaflk = betaBranch * iotaBranch * leaflk;

        //leaflk = iotaBranch * (1-betaBranch);
        //    rootlk = betaRoot * iotaRoot * tmp;

        //    rootLikelihoodsSR_[i] = rootlk + leaflk;
        //}

    }
}


void TwoTreeLikelihood_PIP::computeTreeDLikelihood() {
    for (size_t i = 0; i < nbDistinctSites_; i++) {
        Vdouble *leafLikelihoods1_i = &leafLikelihoods1_[i];
        Vdouble *leafLikelihoods2_i = &leafLikelihoods2_[i];
        double dli = 0;
        for (size_t c = 0; c < nbClasses_; c++) {
            VVdouble *dpxy_c = &dpxy_[c];
            double dlic = 0;
            for (size_t x = 0; x < nbStates_; x++) {
                Vdouble *dpxy_c_x = &(*dpxy_c)[x];
                double l1 = (*leafLikelihoods1_i)[x];
                double dlicx = 0;
                for (size_t y = 0; y < nbStates_; y++) {
                    double l2 = (*leafLikelihoods2_i)[y];
                    dlicx += l1 * l2 * (*dpxy_c_x)[y];
                }
                dlic += dlicx * model_->freq(x);
            }
            dli += dlic * rateDistribution_->getProbability(c);
        }
        dLikelihoods_[i] = dli / rootLikelihoodsSR_[i];
    }
}


void TwoTreeLikelihood_PIP::computeTreeD2Likelihood() {
    for (size_t i = 0; i < nbDistinctSites_; i++) {
        Vdouble *leafLikelihoods1_i = &leafLikelihoods1_[i];
        Vdouble *leafLikelihoods2_i = &leafLikelihoods2_[i];
        double d2li = 0;
        for (size_t c = 0; c < nbClasses_; c++) {
            VVdouble *d2pxy_c = &d2pxy_[c];
            double d2lic = 0;
            for (size_t x = 0; x < nbStates_; x++) {
                Vdouble *d2pxy_c_x = &(*d2pxy_c)[x];
                double l1 = (*leafLikelihoods1_i)[x];
                double d2licx = 0;
                for (size_t y = 0; y < nbStates_; y++) {
                    double l2 = (*leafLikelihoods2_i)[y];
                    d2licx += l1 * l2 * (*d2pxy_c_x)[y];
                }
                d2lic += d2licx * model_->freq(x);
            }
            d2li += d2lic * rateDistribution_->getProbability(c);
        }
        d2Likelihoods_[i] = d2li / rootLikelihoodsSR_[i];
    }
}


double TwoTreeLikelihood_PIP::getFirstOrderDerivative(const string &variable) const
throw(Exception) {
    if (!hasParameter(variable))
        throw ParameterNotFoundException("TwoTreeLikelihood_PIP::getFirstOrderDerivative().", variable);
    if (getRateDistributionParameters().hasParameter(variable)) {
        cout << "DEBUB: WARNING!!! Derivatives respective to rate distribution parameter are not implemented." << endl;
        return log(-1.);
    }
    if (getSubstitutionModelParameters().hasParameter(variable)) {
        cout << "DEBUB: WARNING!!! Derivatives respective to substitution model parameters are not implemented." << endl;
        return log(-1.);
    }

    //
    // Computation for branch lengths:
    //

    // Get the node with the branch whose length must be derivated:
    double d = 0;
    for (size_t i = 0; i < nbDistinctSites_; i++) {
        d += rootWeights_[i] * dLikelihoods_[i];
    }
    return -d;
}


double TwoTreeLikelihood_PIP::getSecondOrderDerivative(const string &variable) const
throw(Exception) {
    if (!hasParameter(variable))
        throw ParameterNotFoundException("TwoTreeLikelihood_PIP::getSecondOrderDerivative().", variable);
    if (getRateDistributionParameters().hasParameter(variable)) {
        cout << "DEBUB: WARNING!!! Derivatives respective to rate distribution parameter are not implemented." << endl;
        return log(-1.);
    }
    if (getSubstitutionModelParameters().hasParameter(variable)) {
        cout << "DEBUB: WARNING!!! Derivatives respective to substitution model parameters are not implemented." << endl;
        return log(-1.);
    }

    //
    // Computation for branch lengths:
    //

    // Get the node with the branch whose length must be derivated:
    double d2 = 0;
    for (size_t i = 0; i < nbDistinctSites_; i++) {
        d2 += rootWeights_[i] * (d2Likelihoods_[i] - pow(dLikelihoods_[i], 2));
    }
    return -d2;
}


void UnifiedDistanceEstimation::computeMatrix() throw(NullPointerException) {
    size_t n = sites_->getNumberOfSequences();
    vector<string> names = sites_->getSequencesNames();
    if (dist_ != 0) delete dist_;
    dist_ = new DistanceMatrix(names);
    optimizer_->setVerbose(static_cast<unsigned int>(max(static_cast<int>(verbose_) - 2, 0)));
    for (size_t i = 0; i < n; ++i) {
        (*dist_)(i, i) = 0;
        if (verbose_ == 1) {
            ApplicationTools::displayGauge(i, n - 1, '=');
        }
        for (size_t j = i + 1; j < n; j++) {
            if (verbose_ > 1) {
                ApplicationTools::displayGauge(j - i - 1, n - i - 2, '=');
            }

            size_t d = SymbolListTools::getNumberOfDistinctPositions(sites_->getSequence(i), sites_->getSequence(j));
            size_t g = SymbolListTools::getNumberOfPositionsWithoutGap(sites_->getSequence(i), sites_->getSequence(j));

            AbstractDiscreteRatesAcrossSitesTreeLikelihood *lik;
            if (getModel().getName().find("PIP") != std::string::npos) {
                lik = new TwoTreeLikelihood_PIP(names[i], names[j], *sites_, model_.get(), rateDist_.get(), verbose_ > 3);
                lik->initialize();
                lik->enableDerivatives(false);
                lik->setParameterValue("BrLen", g == 0 ?
                                                dynamic_cast<TwoTreeLikelihood_PIP *>(lik)->getMinimumBranchLength() :
                                                std::max(dynamic_cast<TwoTreeLikelihood_PIP *>(lik)->getMinimumBranchLength(),
                                                         static_cast<double>(d) / static_cast<double>(g)));

            } else {
                lik = new TwoTreeLikelihood(names[i], names[j], *sites_, model_.get(), rateDist_.get(), verbose_ > 3);
                lik->initialize();
                lik->enableDerivatives(true);
                lik->setParameterValue("BrLen", g == 0 ?
                                                dynamic_cast<TwoTreeLikelihood *>(lik)->getMinimumBranchLength() :
                                                std::max(dynamic_cast<TwoTreeLikelihood *>(lik)->getMinimumBranchLength(),
                                                         static_cast<double>(d) / static_cast<double>(g)));
            }

            // Optimization:
            optimizer_->setFunction(lik);
            optimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
            ParameterList params = lik->getBranchLengthsParameters();
            params.addParameters(parameters_);
            optimizer_->init(params);
            optimizer_->optimize();
            // Store results:
            (*dist_)(i, j) = (*dist_)(j, i) = lik->getParameterValue("BrLen");


            delete lik;
        }
        if (verbose_ > 1 && ApplicationTools::message) ApplicationTools::message->endLine();
    }
}

