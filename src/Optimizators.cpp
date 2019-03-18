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
 * @file Optimizators.cpp
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

#include <glog/logging.h>

// From bpp-core
#include <Bpp/Io/BppODiscreteDistributionFormat.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/FivePointsNumericalDerivative.h>


// From bpp-seq:
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/Likelihood/TreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>


#include <TreeRearrangment.hpp>

#include "Optimizators.hpp"
#include "Utils.hpp"
#include "UnifiedTSHTopologySearch.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_PIP.hpp"

using namespace bpp;

namespace bpp {

    std::string Optimizators::OPTIMIZATION_NEWTON = "newton";
    std::string Optimizators::OPTIMIZATION_GRADIENT = "gradient";
    std::string Optimizators::OPTIMIZATION_BRENT = "Brent";
    std::string Optimizators::OPTIMIZATION_BFGS = "BFGS";

    TreeLikelihood *Optimizators::optimizeParameters(
            bpp::AbstractHomogeneousTreeLikelihood *inTL,
            const ParameterList &parameters,
            std::map<std::string, std::string> &params,
            const std::string &suffix,
            bool suffixIsOptional,
            bool verbose,
            int warn)
    throw(Exception) {


        bpp::AbstractHomogeneousTreeLikelihood *tl = inTL;

        // -------------------------------------------------------------------------
        //  Entry point for optimization routines (both numerical and topology)
        std::string optimization = ApplicationTools::getStringParameter("optimization", params, "FullD(derivatives=Newton)", suffix, suffixIsOptional,
                                                                        warn);
        if (optimization == "None" || optimization == "none") {
            ApplicationTools::displayResult("Optimisations requested:", optimization);
            return tl;
        }

        // -------------------------------------------------------------------------
        // Parsing arguments
        std::string optName;
        std::map<std::string, std::string> optArgs;
        KeyvalTools::parseProcedure(optimization, optName, optArgs);

        // -------------------------------------------------------------------------
        // Verbosity of the optimization routines
        auto optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional, warn + 1);

        // -------------------------------------------------------------------------
        // Message handler
        std::string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional, "none",
                                                            warn + 1);
        auto *messageHandler = static_cast<OutputStream *>((mhPath == "none") ? 0 : (mhPath == "std") ? ApplicationTools::message.get()
                                                                                                      : new StlOutputStream(
                        new std::ofstream(mhPath.c_str(), std::ios::out)));

        if (verbose) ApplicationTools::displayResult("Numerical opt. | Message handler", mhPath);

        // -------------------------------------------------------------------------
        // Profiler
        std::string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional, "none",
                                                            warn + 1);
        auto *profiler = static_cast<OutputStream *>((prPath == "none") ? nullptr : (prPath == "std") ? ApplicationTools::message.get()
                                                                                                      : new StlOutputStream(
                        new std::ofstream(prPath.c_str(), std::ios::app)));
        if (profiler) profiler->setPrecision(20);

        if (verbose) ApplicationTools::displayResult("Numerical opt. | Optimizator profiler", prPath);

        // -------------------------------------------------------------------------
        // Scaling tree topology
        bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, warn + 1);
        if (scaleFirst) {
            // We scale the tree before optimizing each branch length separately:
            if (verbose) ApplicationTools::displayMessage("Scaling the tree before optimizing each branch length separately.");


            double tolerance = ApplicationTools::getDoubleParameter("optimization.scale_first.tolerance", params, .0001, suffix, suffixIsOptional,
                                                                    warn + 1);
            if (verbose) ApplicationTools::displayResult("Numerical opt. | Scaling tolerance:", tolerance);


            auto nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.scale_first.max_number_f_eval", params, 1000000, suffix,
                                                                          suffixIsOptional, warn + 1);
            if (verbose) ApplicationTools::displayResult("Numerical opt. | Scaling max # f eval:", nbEvalMax);


            OptimizationTools::optimizeTreeScale(tl, tolerance, nbEvalMax, messageHandler, profiler);
            if (verbose) ApplicationTools::displayResult("Numerical opt. | New tree likelihood:", -tl->getValue());

        }

        // -------------------------------------------------------------------------
        // Ignoring parameters: should I ignore some parameters?
        ParameterList parametersToEstimate = parameters;
        vector<string> parNames = parametersToEstimate.getParameterNames();

        if (params.find("optimization.ignore_parameter") != params.end())
            throw Exception("optimization.ignore_parameter is deprecated, use optimization.ignore_parameters instead!");
        string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", params, "", suffix, suffixIsOptional, warn + 1);
        StringTokenizer st(paramListDesc, ",");
        while (st.hasMoreToken()) {
            try {
                string param = st.nextToken();
                if (param == "BrLen") {
                    vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
                    parametersToEstimate.deleteParameters(vs);
                    if (verbose) ApplicationTools::displayResult("Numerical opt. | Parameter ignored", std::string("Branch lengths"));

                } else if (param == "Ancient") {
                    auto *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (!nhtl)
                        ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");

                    else {
                        vector<string> vs = nhtl->getRootFrequenciesParameters().getParameterNames();
                        parametersToEstimate.deleteParameters(vs);
                    }
                    ApplicationTools::displayResult("Numerical opt. | Parameter ignored", string("Root frequencies"));

                } else if (param == "Model") {
                    vector<string> vs;
                    vector<string> vs1 = tl->getSubstitutionModelParameters().getParameterNames();
                    auto *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (nhtl != nullptr) {
                        vector<string> vs2 = nhtl->getRootFrequenciesParameters().getParameterNames();
                        VectorTools::diff(vs1, vs2, vs);
                    } else
                        vs = vs1;

                    parametersToEstimate.deleteParameters(vs);
                    ApplicationTools::displayResult("Numerical opt. | Parameter ignored", string("Model"));

                } else if (param.find('*') != string::npos) {
                    vector<string> vs = ApplicationTools::matchingParameters(param, parNames);

                    for (auto it = vs.begin(); it != vs.end(); it++) {
                        parametersToEstimate.deleteParameter(*it);
                        ApplicationTools::displayResult("Numerical opt. | Parameter ignored", *it);


                    }
                } else {
                    parametersToEstimate.deleteParameter(param);
                    ApplicationTools::displayResult("Numerical opt. | Parameter ignored", param);

                }
            }
            catch (ParameterNotFoundException &pnfe) {
                ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");

            }
        }
        // -------------------------------------------------------------------------
        // Constrains: should I constrain some parameters?
        vector<string> parToEstNames = parametersToEstimate.getParameterNames();
        if (params.find("optimization.constrain_parameter") != params.end())
            throw Exception("optimization.constrain_parameter is deprecated, use optimization.constrain_parameters instead!");
        paramListDesc = ApplicationTools::getStringParameter("optimization.constrain_parameters", params, "", suffix, suffixIsOptional, warn + 1);

        string constraint;
        string pc, param;

        StringTokenizer st2(paramListDesc, ",");
        while (st2.hasMoreToken()) {
            try {
                pc = st2.nextToken();
                std::string::size_type index = pc.find('=');
                if (index == std::string::npos)
                    throw Exception("PhylogeneticsApplicationTools::optimizeParamaters. Bad constrain syntax, should contain `=' symbol: " + pc);
                param = pc.substr(0, index);
                constraint = pc.substr(index + 1);
                IntervalConstraint ic(constraint);

                vector<string> parNames2;

                if (param == "BrLen")
                    parNames2 = tl->getBranchLengthsParameters().getParameterNames();
                else if (param == "Ancient") {
                    auto *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (!nhtl)
                        ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");

                    else {
                        parNames2 = nhtl->getRootFrequenciesParameters().getParameterNames();
                        ApplicationTools::displayResult("Numerical opt. | Parameter ignored", string("Root frequencies"));
                    }
                } else if (param == "Model") {
                    vector<string> vs1 = tl->getSubstitutionModelParameters().getParameterNames();
                    auto *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (nhtl != nullptr) {
                        vector<string> vs2 = nhtl->getRootFrequenciesParameters().getParameterNames();
                        VectorTools::diff(vs1, vs2, parNames2);
                    } else
                        parNames2 = vs1;
                } else if (param.find('*') != std::string::npos)
                    parNames2 = ApplicationTools::matchingParameters(param, parToEstNames);
                else
                    parNames2.push_back(param);


                for (size_t i = 0; i < parNames2.size(); i++) {
                    Parameter &par = parametersToEstimate.getParameter(parNames2[i]);
                    if (par.hasConstraint()) {
                        par.setConstraint(ic & (*par.getConstraint()), true);
                        if (par.getConstraint()->isEmpty())
                            throw Exception("Empty interval for parameter " + parNames[i] + par.getConstraint()->getDescription());
                    } else
                        par.setConstraint(ic.clone(), true);

                    ApplicationTools::displayResult("Numerical opt. | Parameter constrained " + par.getName(), par.getConstraint()->getDescription());
                }
            }
            catch (ParameterNotFoundException &pnfe) {
                ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be constrained!");
            }
            catch (ConstraintException &pnfe) {

                string exception_desc = "Parameter '" + param + "' does not fit the constraint " + constraint;
                throw Exception(exception_desc);
            }
        }


        // -------------------------------------------------------------------------
        // Options for optimization routines
        // -------------------------------------------------------------------------

        // Number of max likelihood evaluations
        auto nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional,
                                                                      warn + 1);
        ApplicationTools::displayResult("Numerical opt. | Max # ML evaluations", TextTools::toString(nbEvalMax));

        // Tolerance
        double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional, warn + 1);
        ApplicationTools::displayResult("Numerical opt. | Tolerance", TextTools::toString(tolerance));

        // Backing up or restoring?
        unique_ptr<BackupListener> backupListener;
        string backupFile = ApplicationTools::getAFilePath("optimization.backup.file", params, false, false, suffix, suffixIsOptional, "none",
                                                           warn + 1);
        if (backupFile != "none") {
            ApplicationTools::displayResult("Parameter opt. | Parameters will be backup to", backupFile);

            backupListener.reset(new BackupListener(backupFile));
            if (FileTools::fileExists(backupFile)) {
                ApplicationTools::displayMessage("A backup file was found! Try to restore parameters from previous run...");

                ifstream bck(backupFile.c_str(), ios::in);
                vector<string> lines = FileTools::putStreamIntoVectorOfStrings(bck);
                double fval = TextTools::toDouble(lines[0].substr(5));
                ParameterList pl = tl->getParameters();
                for (size_t l = 1; l < lines.size(); ++l) {
                    if (!TextTools::isEmpty(lines[l])) {
                        StringTokenizer stp(lines[l], "=");
                        if (stp.numberOfRemainingTokens() != 2) {
                            cerr << "Corrupted backup file!!!" << endl;
                            cerr << "at line " << l << ": " << lines[l] << endl;
                        }
                        string pname = stp.nextToken();
                        string pvalue = stp.nextToken();
                        size_t p = pl.whichParameterHasName(pname);
                        pl.setParameter(p, AutoParameter(pl[p]));
                        pl[p].setValue(TextTools::toDouble(pvalue));
                    }
                }
                bck.close();
                tl->setParameters(pl);
                if (abs(tl->getValue() - fval) > 0.000001)
                    ApplicationTools::displayWarning("Warning, incorrect likelihood value after restoring from backup file.");

                ApplicationTools::displayResult("Restoring log-likelihood", -fval);

            }
        }

        // Topology optimisation
        bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, suffix, suffixIsOptional, warn + 1);
        ApplicationTools::displayResult("Parameter opt. | Optimize topology", optimizeTopo ? "yes" : "no");
        tshlib::TreeSearch *treesearch;

        // Derivatives
        string order = ApplicationTools::getStringParameter("derivatives", optArgs, "Newton", "", true, warn + 1);
        string optMethodDeriv;
        if (order == "Gradient") {
            optMethodDeriv = OptimizationTools::OPTIMIZATION_GRADIENT;
        } else if (order == "Newton") {
            optMethodDeriv = OptimizationTools::OPTIMIZATION_NEWTON;
        } else if (order == "BFGS") {
            optMethodDeriv = OptimizationTools::OPTIMIZATION_BFGS;
        } else if (order == "Brent") {
            optMethodDeriv = OptimizationTools::OPTIMIZATION_BRENT;
        } else
            throw Exception("Unknown derivatives algorithm: '" + order + "'.");

        if (verbose) ApplicationTools::displayResult("Numerical opt. | Optimization method", optName);
        if (verbose) ApplicationTools::displayResult("Numerical opt. | Algorithm Derivable parameters", order);

        // Reparametrization of the likelihood function
        bool reparam = ApplicationTools::getBooleanParameter("optimization.reparametrization", params, false, suffix, suffixIsOptional, warn + 1);
        if (verbose) ApplicationTools::displayResult("Numerical opt. | Reparametrization", (reparam ? "yes" : "no"));

        // See if we should use a molecular clock constraint:
        string clock = ApplicationTools::getStringParameter("optimization.clock", params, "None", suffix, suffixIsOptional, warn + 1);
        if (clock != "None" && clock != "Global")
            throw Exception("Molecular clock option not recognized, should be one of 'Global' or 'None'.");
        bool useClock = (clock == "Global");
        if (useClock && optimizeTopo)
            throw Exception("PhylogeneticsApplicationTools::optimizeParameters. Cannot optimize topology with a molecular clock.");
        if (verbose) ApplicationTools::displayResult("Numerical opt. | Molecular clock", clock);

        // Optimisation precision
        auto nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, warn + 1);
        if (verbose && nstep > 1) ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));

        unsigned int n = 0;
        if ((optName == "D-Brent") || (optName == "D-BFGS") || (optName == "ND-Brent") || (optName == "ND-BFGS")) {
            // Uses Newton-Brent method or Newton-BFGS method
            string optMethodModel;
            if ((optName == "D-Brent") || (optName == "ND-Brent"))
                optMethodModel = OptimizationTools::OPTIMIZATION_BRENT;
            else
                optMethodModel = OptimizationTools::OPTIMIZATION_BFGS;

            parametersToEstimate.matchParametersValues(tl->getParameters());

            if (optimizeTopo) {
                std::string PAR_optim_topology_algorithm = ApplicationTools::getStringParameter("optimization.topology.algorithm", params, "", suffix,
                                                                                                suffixIsOptional, warn + 1);

                // Parse string with tree search algorithm details
                std::string optTopology_MethodName;
                std::map<std::string, std::string> optTopology_MethodDetails;
                KeyvalTools::parseProcedure(PAR_optim_topology_algorithm, optTopology_MethodName, optTopology_MethodDetails);

                // Get method to define starting nodes
                std::string optTopology_StartingNodes_MethodName;
                std::map<std::string, std::string> optTopology_StartingNodes_MethodDetails;

                // Find parameters for tree search, if not found then set to default
                std::string PAR_optim_topology_coverage = ApplicationTools::getStringParameter("coverage", optTopology_MethodDetails, "best-search",
                                                                                               suffix, suffixIsOptional, warn + 1);
                std::string PAR_optim_topology_branchoptim = ApplicationTools::getStringParameter("brlen_optimisation", optTopology_MethodDetails,
                                                                                                  "Brent", suffix, suffixIsOptional, warn + 1);
                int PAR_optim_topology_maxcycles = ApplicationTools::getIntParameter("max_cycles", optTopology_MethodDetails, 100, suffix,
                                                                                     suffixIsOptional, warn + 1);
                double PAR_optim_topology_tolerance = ApplicationTools::getDoubleParameter("tolerance", optTopology_MethodDetails, 0.001, suffix,
                                                                                           suffixIsOptional, warn + 1);
                int PAR_optim_topology_threads = ApplicationTools::getIntParameter("threads", optTopology_MethodDetails, 1, suffix,
                                                                                     suffixIsOptional, warn + 1);
                std::string PAR_lkmove = ApplicationTools::getStringParameter("optimization.topology.likelihood", params, "bothways", "", true, true);

                // Prepare settings for the tree-search object (method + coverage)
                tshlib::TreeSearchHeuristics tsh = tshlib::TreeSearchHeuristics::nosearch;

                if (optTopology_MethodName == "Swap") {
                    tsh = tshlib::TreeSearchHeuristics::swap;
                } else if (optTopology_MethodName == "Phyml") {
                    tsh = tshlib::TreeSearchHeuristics::phyml;
                } else if (optTopology_MethodName == "Mixed") {
                    tsh = tshlib::TreeSearchHeuristics::mixed;
                }
                // Coverage setting
                tshlib::TreeRearrangmentOperations tro = tshlib::TreeRearrangmentOperations::classic_Mixed;

                if (PAR_optim_topology_coverage == "nni-search") {
                    tro = tshlib::TreeRearrangmentOperations::classic_NNI;
                } else if (PAR_optim_topology_coverage == "spr-search") {
                    tro = tshlib::TreeRearrangmentOperations::classic_SPR;
                } else if (PAR_optim_topology_coverage == "tbr-search") {
                    tro = tshlib::TreeRearrangmentOperations::classic_TBR;
                } else {
                    tro = tshlib::TreeRearrangmentOperations::classic_Mixed;
                }

                // if the user requests PhyML-like moves, then the search should cover only spr-like moves
                if (tsh == tshlib::TreeSearchHeuristics::phyml && tro != tshlib::TreeRearrangmentOperations::classic_SPR) {
                    tro = tshlib::TreeRearrangmentOperations::classic_SPR;
                    PAR_optim_topology_coverage = "spr-search";
                    if (verbose)
                        ApplicationTools::displayWarning("PhyML-like moves are supported only for SPR-coverage tree-search.\n           "
                                                         "The execution will continue with the following settings:");
                }

                // Print on screen if verbose level is sufficient
                if (verbose) ApplicationTools::displayResult("Topology optimization | Algorithm", optTopology_MethodName);
                if (verbose) ApplicationTools::displayResult("Topology optimization | Moves class", PAR_optim_topology_coverage);

                int PAR_optim_topology_startnodes = 1;
                if (optTopology_MethodDetails.find("starting_nodes") != optTopology_MethodDetails.end()) {
                    KeyvalTools::parseProcedure(optTopology_MethodDetails["starting_nodes"], optTopology_StartingNodes_MethodName,
                                                optTopology_StartingNodes_MethodDetails);
                    PAR_optim_topology_startnodes = ApplicationTools::getIntParameter("n", optTopology_StartingNodes_MethodDetails, 0, suffix,
                                                                                      suffixIsOptional, warn + 1);
                }

                if (verbose) ApplicationTools::displayResult("Topology optimization | Node seeding", optTopology_StartingNodes_MethodName);
                if (verbose) ApplicationTools::displayResult("Topology optimization | # start nodes", PAR_optim_topology_startnodes);
                if (verbose) ApplicationTools::displayResult("Topology optimization | BrLen optimisation", PAR_optim_topology_branchoptim);
                if (verbose) ApplicationTools::displayResult("Topology optimization | Max # cycles", PAR_optim_topology_maxcycles);
                if (verbose) ApplicationTools::displayResult("Topology optimization | Tolerance", PAR_optim_topology_tolerance);
                //if (verbose) ApplicationTools::displayResult("Topology optimization | Initial lk ", TextTools::toString(-tl->getValue(), 15));

                ApplicationTools::displayMessage("");

                // Instantiate tree-search object
                treesearch = new tshlib::TreeSearch;

                treesearch->setTreeSearchStrategy(tsh, tro);

                // Set initial score of the reference topology (up to this point in the optimisation process)
                treesearch->setScoringMethod(PAR_lkmove);

                // Pass the reference of the likelihood function to score the candidate topologies
                treesearch->setLikelihoodFunc(tl);

                // Set starting node method and number of nodes
                tshlib::StartingNodeHeuristics snh = tshlib::StartingNodeHeuristics::undef;
                if (optTopology_StartingNodes_MethodName == "Greedy") {
                    snh = tshlib::StartingNodeHeuristics::greedy;
                } else if (optTopology_StartingNodes_MethodName == "Hillclimbing") {
                    snh = tshlib::StartingNodeHeuristics::hillclimbing;
                }
                treesearch->setStartingNodeHeuristic(snh, PAR_optim_topology_startnodes);

                // Set stop condition and threshold to reach (either no. iterations or tolerance)
                treesearch->setTolerance(PAR_optim_topology_tolerance);
                treesearch->setMaxCycles(PAR_optim_topology_maxcycles);
                treesearch->setMaxNumThreads(PAR_optim_topology_threads);

                // Execute tree-search
                //treesearch->executeTreeSearch();

                //LOG(INFO) << "[TSH Cycle] Likelihood after tree-search lk=" << std::setprecision(18) << tl->getLogLikelihood();
            }

            // Execute numopt + treesearch iteratively until convergence is reached.
            double initScore;
            double cycleScore;
            double diffScore = std::abs(tl->getLogLikelihood());
            bool interrupt = false;

            while (tolerance < diffScore) {

                if (interrupt == true) { break; }

                initScore = tl->getLogLikelihood();
                ApplicationTools::displayResult("Numerical opt. cycle LK", TextTools::toString(initScore, 15));
                // Execute num-opt
                if ((optName == "D-Brent") || (optName == "D-BFGS")) {
                    n = OptimizationTools::optimizeNumericalParameters(
                            dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), parametersToEstimate,
                            backupListener.get(), nstep, tolerance, nbEvalMax, messageHandler, profiler, reparam, optVerbose, optMethodDeriv,
                            optMethodModel);
                } else {

                    n = Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(
                            dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), parametersToEstimate,
                            backupListener.get(), nstep, tolerance, nbEvalMax, messageHandler, profiler, reparam, optVerbose, optMethodDeriv,
                            optMethodModel);
                }

                if (optimizeTopo) {
                    // Execute tree-search
                    treesearch->executeTreeSearch();

                    if (treesearch->isTreeSearchSuccessful()) {
                        // Recompute the difference
                        cycleScore = tl->getLogLikelihood();
                        diffScore = std::abs(initScore) - std::abs(cycleScore);

                    } else {

                        interrupt = true;
                    }

                } else {

                    break;
                }

            }

        } else if (optName == "FullD") {
            // Uses Newton-raphson algorithm with numerical derivatives when required.
            if (optimizeTopo) {
                bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional,
                                                                         warn + 1);
                auto topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, suffix, suffixIsOptional,
                                                                               warn + 1);
                double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix,
                                                                        suffixIsOptional, warn + 1);
                double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix,
                                                                        suffixIsOptional, warn + 1);
                tl = OptimizationTools::optimizeTreeNNI2(
                        dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl), parametersToEstimate,
                        optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
                        reparam, optVerbose, optMethodDeriv, NNITopologySearch::PHYML);
            }

            parametersToEstimate.matchParametersValues(tl->getParameters());
            n = OptimizationTools::optimizeNumericalParameters2(
                    dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), parametersToEstimate,
                    backupListener.get(), tolerance, nbEvalMax, messageHandler, profiler, reparam, useClock, optVerbose, optMethodDeriv);
        } else
            throw Exception("Unknown optimization method: " + optName);

        string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, warn + 1);

        if (verbose) ApplicationTools::displayResult("\nFinal optimization step", finalMethod);

        Optimizer *finalOptimizer = nullptr;
        if (finalMethod == "none") {}
        else if (finalMethod == "simplex") {
            finalOptimizer = new DownhillSimplexMethod(tl);
        } else if (finalMethod == "powell") {
            finalOptimizer = new PowellMultiDimensions(tl);
        } else if (finalMethod == "bfgs") {

            parametersToEstimate.matchParametersValues(tl->getParameters());

            if ((optName == "D-Brent") || (optName == "D-BFGS")) {
                n = OptimizationTools::optimizeNumericalParameters(
                        dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),
                        parametersToEstimate,
                        backupListener.get(),
                        nstep,
                        tolerance,
                        nbEvalMax,
                        messageHandler,
                        profiler,
                        reparam,
                        optVerbose,
                        optMethodDeriv,
                        OptimizationTools::OPTIMIZATION_BFGS);
            } else {
                n = Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(
                        dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), parametersToEstimate,
                        backupListener.get(),
                        nstep,
                        tolerance,
                        nbEvalMax,
                        messageHandler,
                        profiler,
                        reparam,
                        optVerbose,
                        optMethodDeriv,
                        OptimizationTools::OPTIMIZATION_BFGS);


            }

        } else
            throw Exception("Unknown final optimization method: " + finalMethod);

        if (finalOptimizer) {
            parametersToEstimate.matchParametersValues(tl->getParameters());

            finalOptimizer->setProfiler(profiler);
            finalOptimizer->setMessageHandler(messageHandler);
            finalOptimizer->setMaximumNumberOfEvaluations(nbEvalMax);
            finalOptimizer->getStopCondition()->setTolerance(tolerance);
            finalOptimizer->setVerbose((unsigned int) verbose);
            finalOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
            finalOptimizer->init(parametersToEstimate);
            finalOptimizer->optimize();
            n += finalOptimizer->getNumberOfEvaluations();
            delete finalOptimizer;
        }

        if (verbose) ApplicationTools::displayResult("\nPerformed", TextTools::toString(n) + " function evaluations.");

        if (verbose) ApplicationTools::displayResult("Log likelihood after num/top optimisation", TextTools::toString(-tl->getValue(), 15));


        if (backupFile != "none") {
            string bf = backupFile + ".def";
            rename(backupFile.c_str(), bf.c_str());
        }
        return tl;
    }

    std::string Optimizators::DISTANCEMETHOD_INIT = "init";
    std::string Optimizators::DISTANCEMETHOD_PAIRWISE = "pairwise";
    std::string Optimizators::DISTANCEMETHOD_ITERATIONS = "iterations";


    TreeTemplate<Node> *Optimizators::buildDistanceTreeGeneric(
            UnifiedDistanceEstimation &estimationMethod,
            AgglomerativeDistanceMethod &reconstructionMethod,
            const ParameterList &parametersToIgnore,
            bool optimizeBrLen,
            const std::string &param,
            double tolerance,
            unsigned int tlEvalMax,
            OutputStream *profiler,
            OutputStream *messenger,
            unsigned int verbose) throw(Exception) {

        estimationMethod.resetAdditionalParameters();
        estimationMethod.setVerbose(verbose);

        if (param == DISTANCEMETHOD_PAIRWISE) {

            ParameterList tmp = estimationMethod.getModel().getIndependentParameters();
            tmp.addParameters(estimationMethod.getRateDistribution().getIndependentParameters());
            tmp.deleteParameters(parametersToIgnore.getParameterNames());
            estimationMethod.setAdditionalParameters(tmp);
        }

        TreeTemplate<Node> *tree = nullptr;
        TreeTemplate<Node> *previousTree = nullptr;

        bool test = true;
        while (test) {
            // Compute matrice:
            if (verbose > 0)
                ApplicationTools::displayTask("Estimating distance matrix", true);
            estimationMethod.computeMatrix();
            DistanceMatrix *matrix = estimationMethod.getMatrix();
            if (verbose > 0)
                ApplicationTools::displayTaskDone();

            // Compute tree:
            if (matrix->size() == 2) {
                //Special case, there is only one possible tree:
                auto *n1 = new Node(0);
                auto *n2 = new Node(1, matrix->getName(0));
                n2->setDistanceToFather((*matrix)(0, 0) / 2.);
                auto *n3 = new Node(2, matrix->getName(1));
                n3->setDistanceToFather((*matrix)(0, 0) / 2.);
                n1->addSon(n2);
                n1->addSon(n3);
                tree = new TreeTemplate<Node>(n1);
                break;
            }
            if (verbose > 0)
                ApplicationTools::displayTask("Building tree");
            reconstructionMethod.setDistanceMatrix(*matrix);
            reconstructionMethod.computeTree();
            previousTree = tree;
            delete matrix;
            tree = dynamic_cast<TreeTemplate<Node> *>(reconstructionMethod.getTree());
            if (verbose > 0)
                ApplicationTools::displayTaskDone();
            if (previousTree && verbose > 0) {
                int rf = TreeTools::robinsonFouldsDistance(*previousTree, *tree, false);
                ApplicationTools::displayResult("Topo. distance with previous iteration", TextTools::toString(rf));
                test = (rf == 0);
                delete previousTree;
            }
            if (param != DISTANCEMETHOD_ITERATIONS)
                break;  // Ends here.

            // Now, re-estimate parameters:
            unique_ptr<TransitionModel> model(estimationMethod.getModel().clone());
            unique_ptr<DiscreteDistribution> rdist(estimationMethod.getRateDistribution().clone());

            bpp::AbstractHomogeneousTreeLikelihood *tl = nullptr;
            tshlib::Utree *utree;
            UtreeBppUtils::treemap tm;

            if (estimationMethod.getModel().getName().find("PIP") == string::npos) {
                tl = new RHomogeneousTreeLikelihood_Generic(*tree, *estimationMethod.getData(), model.get(), rdist.get(), true, verbose > 1);
            } else {
                std::map<std::string, std::string> default_map;
                tree->setNodeName(tree->getRootId(), "root");
                UtreeBppUtils::renameInternalNodes(tree);
                utree = new tshlib::Utree();
                UtreeBppUtils::convertTree_b2u(tree, utree, tm);
                utree->addVirtualRootNode();
                // Once the tree has the root, then map it as well
                tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode->getVnode_id()));

                tl = new bpp::UnifiedTSHomogeneousTreeLikelihood_PIP(*tree, *estimationMethod.getData(), model.get(), rdist.get(), utree, &tm, true,
                                                                     default_map, "", false, false, false);
            }

            tl->initialize();
            ParameterList parameters = tl->getParameters();
            if (!optimizeBrLen) {
                vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
                parameters.deleteParameters(vs);
            }
            parameters.deleteParameters(parametersToIgnore.getParameterNames());
            if (estimationMethod.getModel().getName().find("PIP") == string::npos) {
                OptimizationTools::optimizeNumericalParameters(tl, parameters, NULL, 0, tolerance, tlEvalMax, messenger, profiler,
                                                               verbose > 0 ? verbose - 1 : 0);
            } else {

                Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(
                        dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), parameters, NULL, 0, tolerance, tlEvalMax, messenger, profiler,
                        false, verbose > 0 ? verbose - 1 : 0, OptimizationTools::OPTIMIZATION_BFGS, OptimizationTools::OPTIMIZATION_BFGS);

                //OptimizationTools::optimizeNumericalParameters(tl, parameters, NULL, 0, tolerance, tlEvalMax, messenger, profiler,
                //                                               verbose > 0 ? verbose - 1 : 0);
            }
            if (verbose > 0) {
                ParameterList tmp = tl->getSubstitutionModelParameters();
                for (unsigned int i = 0; i < tmp.size(); i++) {
                    ApplicationTools::displayResult(tmp[i].getName(), TextTools::toString(tmp[i].getValue()));
                }
                tmp = tl->getRateDistribution()->getParameters();
                for (unsigned int i = 0; i < tmp.size(); i++) {
                    ApplicationTools::displayResult(tmp[i].getName(), TextTools::toString(tmp[i].getValue()));
                }
            }

            delete utree;
        }
        return tree;
    }

    TreeTemplate<Node> *Optimizators::buildDistanceTreeGenericFromDistanceMatrix(DistanceMatrix *dmatrix,
                                                                                 AgglomerativeDistanceMethod &reconstructionMethod,
                                                                                 unsigned int verbose) {

        TreeTemplate<Node> *tree = nullptr;

        // Compute matrice:
        if (verbose > 0)
            ApplicationTools::displayTask("Importing distance matrix", true);
        DistanceMatrix *matrix = dmatrix;
        if (verbose > 0)
            ApplicationTools::displayTaskDone();

        // Compute tree:
        if (matrix->size() == 2) {
            //Special case, there is only one possible tree:
            Node *n1 = new Node(0);
            Node *n2 = new Node(1, matrix->getName(0));
            n2->setDistanceToFather((*matrix)(0, 0) / 2.);
            Node *n3 = new Node(2, matrix->getName(1));
            n3->setDistanceToFather((*matrix)(0, 0) / 2.);
            n1->addSon(n2);
            n1->addSon(n3);
            tree = new TreeTemplate<Node>(n1);
        }
        if (verbose > 0)
            ApplicationTools::displayTask("Building tree");
        reconstructionMethod.setDistanceMatrix(*matrix);
        reconstructionMethod.computeTree();
        delete matrix;
        tree = dynamic_cast<TreeTemplate<Node> *>(reconstructionMethod.getTree());
        if (verbose > 0)
            ApplicationTools::displayTaskDone();


        return tree;
    }


    unsigned int Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(
            DiscreteRatesAcrossSitesTreeLikelihood *tl,
            const ParameterList &parameters,
            OptimizationListener *listener,
            unsigned int nstep,
            double tolerance,
            unsigned int tlEvalMax,
            OutputStream *messageHandler,
            OutputStream *profiler,
            bool reparametrization,
            unsigned int verbose,
            const std::string &optMethodDeriv,
            const std::string &optMethodModel)
    throw(Exception) {
        DerivableSecondOrder *f = tl;
        ParameterList pl = parameters;

        // Shall we reparametrize the function to remove constraints?
        unique_ptr<DerivableSecondOrder> frep;
        if (reparametrization) {
            frep.reset(new ReparametrizationDerivableSecondOrderWrapper(f, parameters));
            f = frep.get();

            // Reset parameters to remove constraints:
            pl = f->getParameters().subList(parameters.getParameterNames());
        }

        // ///////////////
        // Build optimizer:

        // Branch lengths (via numerical derivatives)

        MetaOptimizerInfos *desc = new MetaOptimizerInfos();
        MetaOptimizer *poptimizer = 0;
        AbstractNumericalDerivative *fnum = new ThreePointsNumericalDerivative(f);
        AbstractNumericalDerivative *fnum5 = new FivePointsNumericalDerivative(f);

        if (optMethodDeriv == OPTIMIZATION_GRADIENT)
            desc->addOptimizer("Branch length parameters", new ConjugateGradientMultiDimensions(fnum5),
                               tl->getBranchLengthsParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
        else if (optMethodDeriv == OPTIMIZATION_NEWTON)
            desc->addOptimizer("Branch length parameters", new PseudoNewtonOptimizer(fnum), tl->getBranchLengthsParameters().getParameterNames(), 2,
                               MetaOptimizerInfos::IT_TYPE_FULL);
        else if (optMethodDeriv == OPTIMIZATION_BFGS)
            desc->addOptimizer("Branch length parameters", new BfgsMultiDimensions(fnum5), tl->getBranchLengthsParameters().getParameterNames(), 2,
                               MetaOptimizerInfos::IT_TYPE_FULL);
        else if (optMethodDeriv == OPTIMIZATION_BRENT)
            desc->addOptimizer("Branch length parameters", new SimpleMultiDimensions(fnum5), tl->getBranchLengthsParameters().getParameterNames(), 2,
                               MetaOptimizerInfos::IT_TYPE_FULL);
        else
            throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodDeriv);

        // Other parameters

        if (optMethodModel == OPTIMIZATION_BRENT) {
            ParameterList plsm = parameters.getCommonParametersWith(tl->getSubstitutionModelParameters());
            desc->addOptimizer("Substitution model parameter", new SimpleMultiDimensions(fnum5), plsm.getParameterNames(), 0,
                               MetaOptimizerInfos::IT_TYPE_STEP);


            ParameterList plrd = parameters.getCommonParametersWith(tl->getRateDistributionParameters());
            desc->addOptimizer("Rate distribution parameter", new SimpleMultiDimensions(fnum5), plrd.getParameterNames(), 0,
                               MetaOptimizerInfos::IT_TYPE_STEP);
            poptimizer = new MetaOptimizer(fnum5, desc, nstep);

        } else if (optMethodModel == OPTIMIZATION_BFGS) {
            vector<string> vNameDer;

            ParameterList plsm = parameters.getCommonParametersWith(tl->getSubstitutionModelParameters());
            vNameDer = plsm.getParameterNames();

            ParameterList plrd = parameters.getCommonParametersWith(tl->getRateDistributionParameters());

            vector<string> vNameDer2 = plrd.getParameterNames();

            vNameDer.insert(vNameDer.begin(), vNameDer2.begin(), vNameDer2.end());
            fnum->setParametersToDerivate(vNameDer);

            desc->addOptimizer("Rate & model distribution parameters", new BfgsMultiDimensions(fnum), vNameDer, 1, MetaOptimizerInfos::IT_TYPE_FULL);
            poptimizer = new MetaOptimizer(fnum, desc, nstep);
        } else
            throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodModel);

        poptimizer->setVerbose(verbose);
        poptimizer->setProfiler(profiler);
        poptimizer->setMessageHandler(messageHandler);
        poptimizer->setMaximumNumberOfEvaluations(tlEvalMax);
        poptimizer->getStopCondition()->setTolerance(tolerance);
        poptimizer->getDefaultStopCondition()->setTolerance(tolerance);

        // Optimize TreeLikelihood function:
        poptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        NaNListener *nanListener = new NaNListener(poptimizer, tl);
        poptimizer->addOptimizationListener(nanListener);
        if (listener)
            poptimizer->addOptimizationListener(listener);
        poptimizer->init(pl);
        poptimizer->optimize();

        if (verbose > 0)
            ApplicationTools::displayMessage("\n");

        // We're done.
        unsigned int nb = poptimizer->getNumberOfEvaluations();
        delete poptimizer;
        return nb;
    }


}
