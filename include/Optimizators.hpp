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
 * @file Optimizators.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 24 02 2018
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
#ifndef CASTOR_OPTIMIZATORS_HPP
#define CASTOR_OPTIMIZATORS_HPP
#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include "UnifiedDistanceEstimation.hpp"

namespace bpp {

    class Optimizators {

    public:

        Optimizators();

        virtual ~Optimizators();

        /**
       * @brief Optimize parameters according to options.
       *
       * @param tl               The TreeLikelihood function to optimize.
       * @param parameters       The initial list of parameters to optimize.
       *                         Use tl->getIndependentParameters() in order to estimate all parameters.
       * @param params           The attribute map where options may be found.
       * @param suffix           A suffix to be applied to each attribute name.
       * @param suffixIsOptional Tell if the suffix is absolutely required.
       * @param verbose          Print some info to the 'message' output stream.
       * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
       * @throw Exception        Any exception that may happen during the optimization process.
       * @return A pointer toward the final likelihood object.
       * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
       * clone this object. We may change this bahavior in the future...
       * You hence should write something like
       * @code
       * tl = PhylogeneticsApplicationTools::optimizeParameters(tl, ...);
       * @endcode
       */
        static TreeLikelihood *optimizeParameters(
                bpp::AbstractHomogeneousTreeLikelihood *inTL,
                const ParameterList &parameters,
                std::map<std::string, std::string> &params,
                const std::string &suffix = "",
                bool suffixIsOptional = true,
                bool verbose = true,
                int warn = 1) throw(Exception);


        /**
       * @brief Optimize numerical parameters (branch length, substitution model & rate distribution) of a TreeLikelihood function.
       *
       * Uses Newton's method for branch length and Brent or BFGS one dimensional method for other parameters.
       *
       * A condition over function values is used as a stop condition for the algorithm.
       *
       * @see BrentOneDimension, BFGSMultiDimensions
       *
       * @param tl             A pointer toward the TreeLikelihood object to optimize.
       * @param parameters     The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
       * @param listener       A pointer toward an optimization listener, if needed.
       * @param nstep          The number of progressive steps to perform (see NewtonBrentMetaOptimizer). 1 means full precision from start.
       * @param tolerance      The tolerance to use in the algorithm.
       * @param tlEvalMax      The maximum number of function evaluations.
       * @param messageHandler The massage handler.
       * @param profiler       The profiler.
       * @param reparametrization Tell if parameters should be transformed in order to remove constraints.
       *                          This can improve optimization, but is a bit slower.
       * @param verbose        The verbose level.
       * @param optMethodDeriv Optimization type for derivable parameters (first or second order derivatives).
       * @see OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
       * @param optMethodModel Optimization type for model parameters (Brent or BFGS).
       * @see OPTIMIZATION_BRENT, OPTIMIZATION_BFGS
       * @throw Exception any exception thrown by the Optimizer.
       */
        static unsigned int optimizeNumericalParametersUsingNumericalDerivatives(
                DiscreteRatesAcrossSitesTreeLikelihood *tl,
                const ParameterList &parameters,
                OptimizationListener *listener = 0,
                unsigned int nstep = 1,
                double tolerance = 0.000001,
                unsigned int tlEvalMax = 1000000,
                OutputStream *messageHandler = 0,
                OutputStream *profiler = 0,
                bool reparametrization = false,
                unsigned int verbose = 1,
                const std::string &optMethodDeriv = OPTIMIZATION_NEWTON,
                const std::string &optMethodModel = OPTIMIZATION_BRENT)
        throw(Exception);



        /**
       * @brief Build a tree using a distance method.
       *
       * This method estimate a distance matrix using a DistanceEstimation object, and then builds the phylogenetic tree using a AgglomerativeDistanceMethod object.
       * The main issue here is to estimate non-branch lengths parameters, as substitution model and rate distribution parameters.
       * Three options are provideed here:
       * - DISTANCEMETHOD_INIT (default) keep parameters to there initial value,
       * - DISTANCEMETHOD_PAIRWISE estimated parameters in a pairwise manner, which is standard but not that satisfying...
       * - DISTANCEMETHOD_ITERATIONS uses Ninio et al's iterative algorithm, which uses Maximum Likelihood to estimate these parameters, and then update the distance matrix.
       * Ninio M, Privman E, Pupko T, Friedman N.
       * Phylogeny reconstruction: increasing the accuracy of pairwise distance estimation using Bayesian inference of evolutionary rates.
       * Bioinformatics. 2007 Jan 15;23(2):e136-41.
       *
       * @param estimationMethod The distance estimation object to use.
       * @param reconstructionMethod The tree reconstruction object to use.
       * @param parametersToIgnore A list of parameters to ignore while optimizing parameters.
       * @param optimizeBrLen Tell if branch lengths should be optimized together with other parameters. This may lead to more accurate parameter estimation, but is slower.
       * @param param String describing the type of optimization to use.
       * @param tolerance Threshold on likelihood for stopping the iterative procedure. Used only with param=DISTANCEMETHOD_ITERATIONS.
       * @param tlEvalMax Maximum number of likelihood computations to perform when optimizing parameters. Used only with param=DISTANCEMETHOD_ITERATIONS.
       * @param profiler Output stream used by optimizer. Used only with param=DISTANCEMETHOD_ITERATIONS.
       * @param messenger Output stream used by optimizer. Used only with param=DISTANCEMETHOD_ITERATIONS.
       * @param verbose Verbose level.
       */
        static TreeTemplate <Node> *buildDistanceTreeGeneric(
                UnifiedDistanceEstimation &estimationMethod,
                AgglomerativeDistanceMethod &reconstructionMethod,
                const ParameterList &parametersToIgnore,
                bool optimizeBrLen = false,
                const std::string &param = DISTANCEMETHOD_INIT,
                double tolerance = 0.000001,
                unsigned int tlEvalMax = 1000000,
                OutputStream *profiler = 0,
                OutputStream *messenger = 0,
                unsigned int verbose = 0) throw(Exception);

        static TreeTemplate <Node> *buildDistanceTreeGenericFromDistanceMatrix(DistanceMatrix *dmatrix,
                                                                               AgglomerativeDistanceMethod &reconstructionMethod,
                                                                               unsigned int verbose);

        static std::string OPTIMIZATION_GRADIENT;
        static std::string OPTIMIZATION_NEWTON;
        static std::string OPTIMIZATION_BRENT;
        static std::string OPTIMIZATION_BFGS;

        static std::string DISTANCEMETHOD_INIT;
        static std::string DISTANCEMETHOD_PAIRWISE;
        static std::string DISTANCEMETHOD_ITERATIONS;
    };
} // end of namespace bpp.


#endif //CASTOR_OPTIMIZATORS_HPP
