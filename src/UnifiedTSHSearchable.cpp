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
 * @file UnifiedTSHSearchable.cpp
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
#include <Bpp/Numeric/Function/SimpleMultiDimensions.h>
#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>
#include <glog/logging.h>
#include <Bpp/Text/KeyvalTools.h>
#include "UnifiedTSHSearchable.hpp"


using namespace bpp;

void UnifiedTSHSearchable::setOptimiser(AbstractHomogeneousTreeLikelihood *lk,
                                                      bool optNumericalDerivatives,
                                                      std::map<std::string, std::string> &params,
                                                      const string &suffix,
                                                      bool suffixIsOptional,
                                                      bool verbose,
                                                      int warn) {

    // Set likelihood function
    lk_ = lk;

    // -------------------------------------------------------------------------
    // Optimisation algorithm

    std::string PAR_optim_topology_algorithm = ApplicationTools::getStringParameter("optimization.topology.algorithm", params, "", suffix,
                                                                                    suffixIsOptional, warn + 1);
    std::string optTopology_MethodName;
    std::map<std::string, std::string> optTopology_MethodDetails;
    KeyvalTools::parseProcedure(PAR_optim_topology_algorithm, optTopology_MethodName, optTopology_MethodDetails);

    optMethodModel_ =ApplicationTools::getStringParameter("brlen_optimisation", optTopology_MethodDetails, "Brent", suffix, suffixIsOptional, warn + 1);


    // -------------------------------------------------------------------------
    // Message handler
    std::string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
    auto *messageHandler = static_cast<OutputStream *>((mhPath == "none") ? 0 : (mhPath == "std") ? ApplicationTools::message.get() : new StlOutputStream(
            new std::ofstream(mhPath.c_str(), std::ios::app)));

    // -------------------------------------------------------------------------
    // Profiler
    std::string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
    auto *profiler = static_cast<OutputStream *>((prPath == "none") ? nullptr : (prPath == "std") ? ApplicationTools::message.get() : new StlOutputStream(
            new std::ofstream(prPath.c_str(), std::ios::app)));
    if (profiler) profiler->setPrecision(20);


    auto nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, warn + 1);

    if (optNumericalDerivatives) {

        AbstractNumericalDerivative *dn_3points = new ThreePointsNumericalDerivative(lk);
        AbstractNumericalDerivative *dn_5points = new FivePointsNumericalDerivative(lk);

        // Initialise branch optimiser
        if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BRENT) {
            optimiser_ = new SimpleMultiDimensions(dn_5points);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BFGS) {
            optimiser_ = new BfgsMultiDimensions(dn_5points);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_NEWTON) {
            optimiser_ = new PseudoNewtonOptimizer(dn_3points);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_GRADIENT) {
            optimiser_ = new ConjugateGradientMultiDimensions(dn_5points);
        }else{
            LOG(FATAL) << "The branch optimisation method specified is unknown: " << optMethodModel_;
        }
    } else {
        // Initialise branch optimiser
        if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BRENT) {
            optimiser_ = new SimpleMultiDimensions(lk);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BFGS) {
            optimiser_ = new BfgsMultiDimensions(lk);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_NEWTON) {
            optimiser_ = new PseudoNewtonOptimizer(lk);
        } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_GRADIENT) {
            optimiser_ = new ConjugateGradientMultiDimensions(lk);
        }else{
            LOG(FATAL) << "The branch optimisation method specified is unknown: " << optMethodModel_;
        }
    }
    optimiser_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimiser_->setProfiler(profiler);
    optimiser_->setMessageHandler(messageHandler);
    optimiser_->setMaximumNumberOfEvaluations(nbEvalMax);
    optimiser_->setVerbose(0);

}


void UnifiedTSHSearchable::fireBranchOptimisation(std::vector<bpp::Node *> extractionNodes) {

    ParameterList parameters;
    // For each node involved in the move, get the corrisponding branch parameter (no root)
    for (auto &bnode:extractionNodes) {
        if (bnode->hasFather()) {
            Parameter brLen = lk_->getParameter("BrLen" + TextTools::toString(bnode->getId()));
            brLen.setName("BrLen" + TextTools::toString(bnode->getId()));
            parameters.addParameter(brLen);
        }
    }

    // set parameters on the likelihood function (inherited)
    lk_->setParameters(parameters);

    // Re-estimate branch length:
    if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BRENT) {
        auto optimiserInstance = dynamic_cast<SimpleMultiDimensions *>(optimiser_);
        optimiserInstance->init(parameters);
        optimiserInstance->doInit(parameters);
        optimiserInstance->getStopCondition()->setTolerance(0.001);
        optimiserInstance->optimize();

        // set parameters on the likelihood function (inherited)
        lk_->setParameters(optimiserInstance->getParameters());

    } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_BFGS) {
        auto optimiserInstance = dynamic_cast<BfgsMultiDimensions *>(optimiser_);
        optimiserInstance->init(parameters);
        optimiserInstance->doInit(parameters);
        optimiserInstance->getStopCondition()->setTolerance(0.001);
        optimiserInstance->optimize();

        // set parameters on the likelihood function (inherited)
        lk_->setParameters(optimiserInstance->getParameters());

    } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_NEWTON) {
        auto optimiserInstance = dynamic_cast<PseudoNewtonOptimizer *>(optimiser_);
        optimiserInstance->init(parameters);
        optimiserInstance->doInit(parameters);
        optimiserInstance->getStopCondition()->setTolerance(0.001);
        optimiserInstance->optimize();

        // set parameters on the likelihood function (inherited)
        lk_->setParameters(optimiserInstance->getParameters());

    } else if (optMethodModel_ == UnifiedTSHSearchable::OPTIMIZATION_GRADIENT) {
        auto optimiserInstance = dynamic_cast<ConjugateGradientMultiDimensions *>(optimiser_);
        optimiserInstance->init(parameters);
        optimiserInstance->doInit(parameters);
        optimiserInstance->getStopCondition()->setTolerance(0.001);
        optimiserInstance->optimize();

        // set parameters on the likelihood function (inherited)
        lk_->setParameters(optimiserInstance->getParameters());
    }


}

