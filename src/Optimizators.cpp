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

    void Optimizators::init(std::map<std::string, std::string> &params,std::string suffix,bool suffixIsOptional,bool verbose,int warn){

        /////////////////////////////////////////
        this->suffix = suffix;
        this->suffixIsOptional = suffixIsOptional;
        this->verbose = verbose;
        this->warn = warn;
        this->params = params;

        /////////////////////////////////////////
        // Verbosity of the optimization routines
        this->optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", this->params, 2, this->suffix, this->suffixIsOptional, this->warn + 1);

        /////////////////////////////////////////
        // Number of max likelihood evaluations
        this->nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", this->params, 1000000, this->suffix, this->suffixIsOptional,this->warn + 1);
        ApplicationTools::displayResult("Numerical opt. | Max # ML evaluations", TextTools::toString(this->nbEvalMax));

        /////////////////////////////////////////
        // Tolerance
        this->tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", this->params, .000001, this->suffix, this->suffixIsOptional, this->warn + 1);
        ApplicationTools::displayResult("Numerical opt. | Tolerance", TextTools::toString(this->tolerance));

        /////////////////////////////////////////
        // Message handler
        this->setMessageHandler();

        /////////////////////////////////////////
        // Profiler
        this->setProfiler();

        /////////////////////////////////////////
        // Set backup file
        this->backupFile  = ApplicationTools::getAFilePath("optimization.backup.file", this->params, false, false, this->suffix, this->suffixIsOptional, "none",this->warn + 1);

        /////////////////////////////////////////////////
        // Method

        this->optimization = ApplicationTools::getStringParameter("optimization", this->params, "FullD(derivatives=Newton)", this->suffix, this->suffixIsOptional, this->warn);

        this->finalMethod = ApplicationTools::getStringParameter("optimization.final", this->params, "none", this->suffix, this->suffixIsOptional, this->warn + 1);

        /////////////////////////////////////////////////
        // Reparametrization of the likelihood function

        this->reparam = ApplicationTools::getBooleanParameter("optimization.reparametrization", this->params, false, this->suffix, this->suffixIsOptional, this->warn + 1);

        /////////////////////////////////////////////////
        this->optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", this->params, true, this->suffix, this->suffixIsOptional,this->warn + 1);

        this->topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", this->params, 1, this->suffix, this->suffixIsOptional,this->warn + 1);

        this->tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", this->params, 100, this->suffix,this->suffixIsOptional, this->warn + 1);

        this->tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", this->params, 100, this->suffix,this->suffixIsOptional, this->warn + 1);

        /////////////////////////////////////////////////
        this->optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", this->params, false, this->suffix, this->suffixIsOptional, this->warn + 1);

        /////////////////////////////////////////////////
        // See if we should use a molecular clock constraint:
        this->clock = ApplicationTools::getStringParameter("optimization.clock", this->params, "None", this->suffix, this->suffixIsOptional, this->warn + 1);

    }

    void Optimizators::scaleTreeTopology(bpp::AbstractHomogeneousTreeLikelihood *tl){

        bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", this->params, false, this->suffix, this->suffixIsOptional, this->warn + 1);

        if (scaleFirst) {

            // We scale the tree before optimizing each branch length separately:
            if (this->verbose){
                ApplicationTools::displayMessage("Scaling the tree before optimizing each branch length separately.");
            }


            this->tolerance = ApplicationTools::getDoubleParameter("optimization.scale_first.tolerance", this->params, .0001, this->suffix, this->suffixIsOptional,this->warn + 1);

            if (this->verbose){
                ApplicationTools::displayResult("Numerical opt. | Scaling tolerance:", this->tolerance);
            }


            this->nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.scale_first.max_number_f_eval", this->params, 1000000, this->suffix,this->suffixIsOptional, this->warn + 1);

            if (this->verbose){
                ApplicationTools::displayResult("Numerical opt. | Scaling max # f eval:", this->nbEvalMax);
            }


            OptimizationTools::optimizeTreeScale(tl, this->tolerance, this->nbEvalMax, this->messageHandler, this->profiler);

            if (this->verbose){
                ApplicationTools::displayResult("Numerical opt. | New tree likelihood:", -tl->getValue());
            }

        }

    }

    void Optimizators::getIgnorParamsList(bpp::AbstractHomogeneousTreeLikelihood *tl,structParams *prms){

        prms->parNames = prms->parametersToEstimate.getParameterNames();

        if (this->params.find("optimization.ignore_parameter") != this->params.end()){
            throw Exception("optimization.ignore_parameter is deprecated, use optimization.ignore_parameters instead!");
        }

        prms->paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", this->params, "", this->suffix, this->suffixIsOptional, this->warn + 1);

        StringTokenizer st(prms->paramListDesc, ",");
        while (st.hasMoreToken()) {
            try {
                string param = st.nextToken();
                if (param == "BrLen") {
                    vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
                    prms->parametersToEstimate.deleteParameters(vs);
                    if (this->verbose){
                        ApplicationTools::displayResult("Numerical opt. | Parameter ignored", std::string("Branch lengths"));
                    }

                } else if (param == "Ancient") {
                    auto *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (!nhtl){
                        ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
                    }else {
                        vector<string> vs = nhtl->getRootFrequenciesParameters().getParameterNames();
                        prms->parametersToEstimate.deleteParameters(vs);
                    }
                    ApplicationTools::displayResult("Numerical opt. | Parameter ignored", string("Root frequencies"));

                } else if (param == "Model") {
                    vector<string> vs;
                    vector<string> vs1 = tl->getSubstitutionModelParameters().getParameterNames();
                    auto *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
                    if (nhtl != nullptr) {
                        vector<string> vs2 = nhtl->getRootFrequenciesParameters().getParameterNames();
                        VectorTools::diff(vs1, vs2, vs);
                    } else{
                        vs = vs1;
                    }

                    prms->parametersToEstimate.deleteParameters(vs);
                    ApplicationTools::displayResult("Numerical opt. | Parameter ignored", string("Model"));

                } else if (param.find('*') != string::npos) {

                    vector<string> vs = ApplicationTools::matchingParameters(param, prms->parNames);

                    for (auto it = vs.begin(); it != vs.end(); it++) {
                        prms->parametersToEstimate.deleteParameter(*it);
                        ApplicationTools::displayResult("Numerical opt. | Parameter ignored", *it);


                    }

                } else {
                    prms->parametersToEstimate.deleteParameter(param);
                    ApplicationTools::displayResult("Numerical opt. | Parameter ignored", param);

                }
            }catch (ParameterNotFoundException &pnfe) {
                ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");

            }
        }

    }

    void Optimizators::constrainParameters(bpp::AbstractHomogeneousTreeLikelihood *tl,structParams *prms){

        vector<string> parToEstNames = prms->parametersToEstimate.getParameterNames();

        if (params.find("optimization.constrain_parameter") != params.end()){
            throw Exception("optimization.constrain_parameter is deprecated, use optimization.constrain_parameters instead!");
        }

        prms->paramListDesc = ApplicationTools::getStringParameter("optimization.constrain_parameters", params, "", this->suffix, this->suffixIsOptional, this->warn + 1);

        std::string constraint;
        std::string pc;
        std::string param;

        StringTokenizer st2(prms->paramListDesc, ",");
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
                    Parameter &par = prms->parametersToEstimate.getParameter(parNames2[i]);
                    if (par.hasConstraint()) {
                        par.setConstraint(ic & (*par.getConstraint()), true);
                        if (par.getConstraint()->isEmpty())
                            throw Exception("Empty interval for parameter " + prms->parNames[i] + par.getConstraint()->getDescription());
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

    }

    void Optimizators::setBackUp(bpp::AbstractHomogeneousTreeLikelihood *tl){

        ApplicationTools::displayResult("Parameter opt. | Parameters will be backup to", this->backupFile);

        if(this->backupListener){
            delete this->backupListener;
            this->backupListener = nullptr;
        }

        this->backupListener = new BackupListener(this->backupFile);

        if (FileTools::fileExists(this->backupFile)) {

            ApplicationTools::displayMessage("A backup file was found! Try to restore parameters from previous run...");

            ifstream bck(this->backupFile.c_str(), ios::in);

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

            if (abs(tl->getValue() - fval) > 0.000001){
                ApplicationTools::displayWarning("Warning, incorrect likelihood value after restoring from backup file.");
            }

            ApplicationTools::displayResult("Restoring log-likelihood", -fval);

        }

    }

    structTreeSearchParams* Optimizators::getTreeSearchParameters(std::map<std::string, std::string> &optTopology_MethodDetails){

        structTreeSearchParams *stsp = nullptr;

        stsp = new structTreeSearchParams();

        // Get method to define starting nodes
        // Find parameters for tree search, if not found then set to default
        stsp->PAR_optim_topology_coverage = ApplicationTools::getStringParameter("coverage", optTopology_MethodDetails, "best-search",this->suffix, this->suffixIsOptional, this->warn + 1);

        stsp->PAR_optim_topology_branchoptim = ApplicationTools::getStringParameter("brlen_optimisation", optTopology_MethodDetails,"Brent", this->suffix, this->suffixIsOptional, this->warn + 1);

        stsp->PAR_optim_topology_maxcycles = ApplicationTools::getIntParameter("max_cycles", optTopology_MethodDetails, 100, this->suffix,this->suffixIsOptional, this->warn + 1);

        stsp->PAR_optim_topology_tolerance = ApplicationTools::getDoubleParameter("tolerance", optTopology_MethodDetails, 0.001, this->suffix,this->suffixIsOptional, this->warn + 1);

        stsp->PAR_optim_topology_threads = ApplicationTools::getIntParameter("threads", optTopology_MethodDetails, 1, this->suffix,this->suffixIsOptional, this->warn + 1);

        stsp->PAR_lkmove = ApplicationTools::getStringParameter("optimization.topology.likelihood", this->params, "bothways", "", true, true);

        return stsp;
    }

    void Optimizators::topologyOptimizer(bpp::AbstractHomogeneousTreeLikelihood *tl,tshlib::TreeSearch *treesearch){

        std::string optTopology_MethodName = "";
        std::string PAR_optim_topology_algorithm = "";
        std::string optTopology_StartingNodes_MethodName = "";

        std::map<std::string, std::string> optTopology_StartingNodes_MethodDetails;
        std::map<std::string, std::string> optTopology_MethodDetails;

        tshlib::TreeSearchHeuristics tsh;
        tshlib::TreeRearrangmentOperations tro;

        tshlib::StartingNodeHeuristics snh;

        structTreeSearchParams *stsp = nullptr;

        // Parse string with tree search algorithm details
        PAR_optim_topology_algorithm = ApplicationTools::getStringParameter("optimization.topology.algorithm", this->params, "", this->suffix, this->suffixIsOptional, this->warn + 1);

        KeyvalTools::parseProcedure(PAR_optim_topology_algorithm, optTopology_MethodName, optTopology_MethodDetails);

        // Get method to define starting nodes
        stsp = this->getTreeSearchParameters(optTopology_MethodDetails);

        // Prepare settings for the tree-search object (method + coverage)
        tsh = tshlib::TreeSearchHeuristics::nosearch;

        if (optTopology_MethodName == "Swap") {
            tsh = tshlib::TreeSearchHeuristics::swap;
        } else if (optTopology_MethodName == "Phyml") {
            tsh = tshlib::TreeSearchHeuristics::phyml;
        } else if (optTopology_MethodName == "Mixed") {
            tsh = tshlib::TreeSearchHeuristics::mixed;
        }

        // Coverage setting
        if (stsp->PAR_optim_topology_coverage == "nni-search") {
            tro = tshlib::TreeRearrangmentOperations::classic_NNI;
        } else if (stsp->PAR_optim_topology_coverage == "spr-search") {
            tro = tshlib::TreeRearrangmentOperations::classic_SPR;
        } else if (stsp->PAR_optim_topology_coverage == "tbr-search") {
            tro = tshlib::TreeRearrangmentOperations::classic_TBR;
        } else {
            tro = tshlib::TreeRearrangmentOperations::classic_Mixed;
        }

        // if the user requests PhyML-like moves, then the search should cover only spr-like moves
        if (tsh == tshlib::TreeSearchHeuristics::phyml && tro != tshlib::TreeRearrangmentOperations::classic_SPR) {

            tro = tshlib::TreeRearrangmentOperations::classic_SPR;

            stsp->PAR_optim_topology_coverage = "spr-search";

            if (verbose){
                ApplicationTools::displayWarning("PhyML-like moves are supported only for SPR-coverage tree-search. The execution will continue with the following settings:");
            }

        }

        // Print on screen if verbose level is sufficient
        if (verbose){
            ApplicationTools::displayResult("Topology optimization | Algorithm", optTopology_MethodName);
            ApplicationTools::displayResult("Topology optimization | Moves class", stsp->PAR_optim_topology_coverage);
        }

        int PAR_optim_topology_startnodes = 1;
        if (optTopology_MethodDetails.find("starting_nodes") != optTopology_MethodDetails.end()) {
            KeyvalTools::parseProcedure(optTopology_MethodDetails["starting_nodes"], optTopology_StartingNodes_MethodName,optTopology_StartingNodes_MethodDetails);
        }

        if (verbose){
            ApplicationTools::displayResult("Topology optimization | Node seeding", optTopology_StartingNodes_MethodName);
            ApplicationTools::displayResult("Topology optimization | # start nodes", PAR_optim_topology_startnodes);
            ApplicationTools::displayResult("Topology optimization | BrLen optimisation", stsp->PAR_optim_topology_branchoptim);
            ApplicationTools::displayResult("Topology optimization | Max # cycles", stsp->PAR_optim_topology_maxcycles);
            ApplicationTools::displayResult("Topology optimization | Tolerance", stsp->PAR_optim_topology_tolerance);
        }

        ApplicationTools::displayMessage("");

        // Instantiate tree-search object
        treesearch = new tshlib::TreeSearch;

        treesearch->setTreeSearchStrategy(tsh, tro);

        // Set initial score of the reference topology (up to this point in the optimisation process)
        treesearch->setScoringMethod(stsp->PAR_lkmove);

        // Pass the reference of the likelihood function to score the candidate topologies
        treesearch->setLikelihoodFunc(tl);

        // Set starting node method and number of nodes
        snh = tshlib::StartingNodeHeuristics::undef;
        if (optTopology_StartingNodes_MethodName == "Greedy") {
            snh = tshlib::StartingNodeHeuristics::greedy;
        } else if (optTopology_StartingNodes_MethodName == "Hillclimbing") {
            snh = tshlib::StartingNodeHeuristics::hillclimbing;
        }

        treesearch->setStartingNodeHeuristic(snh, PAR_optim_topology_startnodes);

        // Set stop condition and threshold to reach (either no. iterations or tolerance)
        treesearch->setTolerance(stsp->PAR_optim_topology_tolerance);
        treesearch->setMaxCycles(stsp->PAR_optim_topology_maxcycles);
        treesearch->setMaxNumThreads(stsp->PAR_optim_topology_threads);

        delete stsp;

    }

    void Optimizators::optimizeTopoBrentBFGS(bpp::AbstractHomogeneousTreeLikelihood *tl,tshlib::TreeSearch *treesearch,structParams *prms){

        // Uses Newton-Brent method or Newton-BFGS method
        if ((this->optName == "D-Brent") || (this->optName == "ND-Brent")){
            this->optMethodModel = OptimizationTools::OPTIMIZATION_BRENT;
        }else{
            this->optMethodModel = OptimizationTools::OPTIMIZATION_BFGS;
        }

        prms->parametersToEstimate.matchParametersValues(tl->getParameters());

        if (this->optimizeTopo) {
            this->topologyOptimizer(tl,treesearch);
        }

        // Execute numopt + treesearch iteratively until convergence is reached.
        this->iterativeOptimizeTreeAndParams(tl,treesearch,prms);

    }

    void Optimizators::iterativeOptimizeTreeAndParams(bpp::AbstractHomogeneousTreeLikelihood *tl,tshlib::TreeSearch *treesearch,structParams *prms){

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
            if ((this->optName == "D-Brent") || (this->optName == "D-BFGS")) {

                this->n = OptimizationTools::optimizeNumericalParameters(
                        dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), prms->parametersToEstimate,
                        this->backupListener, this->nstep, this->tolerance, this->nbEvalMax, this->messageHandler,
                        this->profiler, this->reparam, this->optVerbose, this->optMethodDeriv,this->optMethodModel);

            } else {

//                this->n = Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(
//                        dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),
//                        prms->parametersToEstimate,
//                        this->backupListener.get(),
//                        this->nstep,
//                        this->tolerance,
//                        this->nbEvalMax,
//                        this->messageHandler,
//                        this->profiler,
//                        this->reparam,
//                        this->optVerbose,
//                        this->optMethodDeriv,
//                        this->optMethodModel);
                this->n = Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(
                        dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),prms->parametersToEstimate);
            }

            if (this->optimizeTopo) {

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

    }

    void Optimizators::optimizeTopoFullD(bpp::AbstractHomogeneousTreeLikelihood *tl,structParams *prms){

        // Uses Newton-raphson algorithm with numerical derivatives when required.
        if (this->optimizeTopo) {

            tl = OptimizationTools::optimizeTreeNNI2(dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl), prms->parametersToEstimate,
                    this->optNumFirst, this->tolBefore, this->tolDuring, this->nbEvalMax, this->topoNbStep, this->messageHandler, this->profiler,
                    this->reparam, this->optVerbose, this->optMethodDeriv, NNITopologySearch::PHYML);

        }

        prms->parametersToEstimate.matchParametersValues(tl->getParameters());

        this->n = OptimizationTools::optimizeNumericalParameters2(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),
                prms->parametersToEstimate,this->backupListener, this->tolerance, this->nbEvalMax, this->messageHandler,
                this->profiler, this->reparam, this->useClock, this->optVerbose, this->optMethodDeriv);

    }

    tshlib::TreeSearch* Optimizators::optimizeTopology(bpp::AbstractHomogeneousTreeLikelihood *tl,structParams *prms){

        tshlib::TreeSearch *treesearch = nullptr;

        std::string order = "";

        this->nstep = 0;

        ApplicationTools::displayResult("Parameter opt. | Optimize topology", this->optimizeTopo ? "yes" : "no");

        // Derivatives
        order = ApplicationTools::getStringParameter("derivatives", prms->optArgs, "Newton", "", true, this->warn + 1);

        if (order == "Gradient") {
            this->optMethodDeriv = OptimizationTools::OPTIMIZATION_GRADIENT;
        } else if (order == "Newton") {
            this->optMethodDeriv = OptimizationTools::OPTIMIZATION_NEWTON;
        } else if (order == "BFGS") {
            this->optMethodDeriv = OptimizationTools::OPTIMIZATION_BFGS;
        } else if (order == "Brent") {
            this->optMethodDeriv = OptimizationTools::OPTIMIZATION_BRENT;
        } else{
            throw Exception("Unknown derivatives algorithm: '" + order + "'.");
        }

        if (verbose) {
            ApplicationTools::displayResult("Numerical opt. | Optimization method", this->optName);
            ApplicationTools::displayResult("Numerical opt. | Algorithm Derivable parameters", order);
        }

        if (this->verbose){
            ApplicationTools::displayResult("Numerical opt. | Reparametrization", (this->reparam ? "yes" : "no"));
        }

        if (this->clock != "None" && this->clock != "Global"){
            throw Exception("Molecular clock option not recognized, should be one of 'Global' or 'None'.");
        }

        this->useClock = (this->clock == "Global");
        if (this->useClock && this->optimizeTopo){
            throw Exception("PhylogeneticsApplicationTools::optimizeParameters. Cannot optimize topology with a molecular clock.");
        }

        if (this->verbose){
            ApplicationTools::displayResult("Numerical opt. | Molecular clock", this->clock);
        }

        // Optimisation precision
        this->nstep = ApplicationTools::getParameter<unsigned int>("nstep", prms->optArgs, 1, "", true, this->warn + 1);

        if (this->verbose && nstep > 1){
            ApplicationTools::displayResult("# of precision steps", TextTools::toString(this->nstep));
        }

        this->n = 0;

        if ((this->optName == "D-Brent") || (this->optName == "D-BFGS") || (this->optName == "ND-Brent") || (this->optName == "ND-BFGS")) {

            this->optimizeTopoBrentBFGS(tl,treesearch,prms);

        } else if (this->optName == "FullD") {

            this->optimizeTopoFullD(tl,prms);

        } else{

            throw Exception("Unknown optimization method: " + this->optName);

        }

        return treesearch;
    }

    void Optimizators::finalOptimization(bpp::AbstractHomogeneousTreeLikelihood *tl,structParams *prms){

        if (this->verbose){
            ApplicationTools::displayResult("\nFinal optimization step", this->finalMethod);
        }

        if (this->finalMethod == "none") {
            NULL;
        } else if (finalMethod == "simplex") {
            this->finOptimizer = new DownhillSimplexMethod(tl);
        } else if (finalMethod == "powell") {
            this->finOptimizer = new PowellMultiDimensions(tl);
        } else if (finalMethod == "bfgs") {

            prms->parametersToEstimate.matchParametersValues(tl->getParameters());

            this->optMethodModel = OptimizationTools::OPTIMIZATION_BFGS;

            if ((this->optName == "D-Brent") || (this->optName == "D-BFGS")) {

//                this->n = OptimizationTools::optimizeNumericalParameters(
//                        dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),
//                        prms->parametersToEstimate,
//                        this->backupListener.get(),
//                        this->nstep,
//                        this->tolerance,
//                        this->nbEvalMax,
//                        this->messageHandler,
//                        this->profiler,
//                        this->reparam,
//                        this->optVerbose,
//                        this->optMethodDeriv,
//                        OptimizationTools::OPTIMIZATION_BFGS);
                this->n = OptimizationTools::optimizeNumericalParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),prms->parametersToEstimate);

            } else {

//                this->n = Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(
//                        dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),
//                        prms->parametersToEstimate,
//                        this->backupListener.get(),
//                        this->nstep,
//                        this->tolerance,
//                        this->nbEvalMax,
//                        this->messageHandler,
//                        this->profiler,
//                        this->reparam,
//                        this->optVerbose,
//                        this->optMethodDeriv,
//                        OptimizationTools::OPTIMIZATION_BFGS);
                this->n = Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),prms->parametersToEstimate);

            }

        } else{
            throw Exception("Unknown final optimization method: " + this->finalMethod);
        }

        prms->parametersToEstimate.matchParametersValues(tl->getParameters());

        this->finOptimizer->setProfiler(this->profiler);
        this->finOptimizer->setMessageHandler(this->messageHandler);
        this->finOptimizer->setMaximumNumberOfEvaluations(this->nbEvalMax);
        this->finOptimizer->getStopCondition()->setTolerance(this->tolerance);
        this->finOptimizer->setVerbose((unsigned int) this->verbose);
        this->finOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        this->finOptimizer->init(prms->parametersToEstimate);
        this->finOptimizer->optimize();

        this-> n += this->finOptimizer->getNumberOfEvaluations();

    }

    void Optimizators::setMessageHandler(){

        OutputStream * messageHandler = nullptr;
        std::string mhPath = "";

        mhPath = ApplicationTools::getAFilePath("optimization.message_handler", this->params, false, false, this->suffix, this->suffixIsOptional, "none", this->warn + 1);

        messageHandler = static_cast<OutputStream *>((mhPath == "none") ? 0 : (mhPath == "std") ? ApplicationTools::message.get() : new StlOutputStream(new std::ofstream(mhPath.c_str(), std::ios::out)));

        if (this->verbose){
            ApplicationTools::displayResult("Numerical opt. | Message handler", mhPath);
        }

        this->messageHandler = messageHandler;
    }

    void Optimizators::setProfiler(){

        OutputStream *profiler = nullptr;
        std::string prPath = "";

        prPath = ApplicationTools::getAFilePath("optimization.profiler", this->params, false, false, this->suffix, this->suffixIsOptional, "none", this->warn + 1);

        profiler = static_cast<OutputStream *>((prPath == "none") ? nullptr : (prPath == "std") ? ApplicationTools::message.get(): new StlOutputStream(new std::ofstream(prPath.c_str(), std::ios::app)));

        if (profiler){
            profiler->setPrecision(20);
        }

        if (this->verbose){
            ApplicationTools::displayResult("Numerical opt. | Optimizator profiler", prPath);
        }

        this->profiler = profiler;
    }

    TreeLikelihood *Optimizators::optimizeParameters(
            bpp::AbstractHomogeneousTreeLikelihood *inTL,
            const ParameterList &parameters)
    throw(Exception) {

        bpp::AbstractHomogeneousTreeLikelihood *tl = inTL;
        tshlib::TreeSearch *treesearch;
        structParams *prms = nullptr;

        /////////////////////////////////////////////////
        //
        prms = new structParams();

        prms->parametersToEstimate = parameters;

        /////////////////////////////////////////////////
        // Check if optimization is requested

        if (this->optimization == "None" || this->optimization == "none") {
            ApplicationTools::displayResult("Optimisations requested:", this->optimization);
            return tl;
        }

        /////////////////////////////////////////////////
        // Parsing arguments

        KeyvalTools::parseProcedure(optimization, this->optName, prms->optArgs);

        /////////////////////////////////////////////////
        // Scaling tree topology

        this->scaleTreeTopology(tl);

        /////////////////////////////////////////////////
        // Ignoring parameters: should I ignore some parameters?

        this->getIgnorParamsList(tl,prms);

        /////////////////////////////////////////////////
        // Constrains: should I constrain some parameters?

        this->constrainParameters(tl,prms);

        /////////////////////////////////////////////////
        // Backing up or restoring?

        if (this->backupFile != "none") {
            this->setBackUp(tl);
        }

        /////////////////////////////////////////////////
        // Topology optimisation
        treesearch = this->optimizeTopology(tl,prms);

        if (this->finalOptimizer) {
            this->finalOptimization(tl,prms);
        }

        /////////////////////////////////////////////////
        if (this->verbose){
            ApplicationTools::displayResult("\nPerformed", TextTools::toString(this->n) + " function evaluations.");

            ApplicationTools::displayResult("Log likelihood after num/top optimisation", TextTools::toString(-tl->getValue(), 15));
        }

        /////////////////////////////////////////////////
        if (this->backupFile != "none") {
            string bf = this->backupFile + ".def";
            rename(this->backupFile.c_str(), bf.c_str());
        }

        /////////////////////////////////////////////////
        delete finalOptimizer;
        delete prms;


        return tl;
    }

    void Optimizators::reEstimateParameters(bpp::AbstractHomogeneousTreeLikelihood *tl,UnifiedDistanceEstimation &estimationMethod,
                                            TreeTemplate<Node> *tree,UtreeBppUtils::treemap &tm,
                                            bool optimizeBrLen,bpp::ParameterList &parametersToIgnore){

        ParameterList parameters;
        bool isPIP = false;

        unique_ptr<TransitionModel> model(estimationMethod.getModel().clone());
        unique_ptr<DiscreteDistribution> rdist(estimationMethod.getRateDistribution().clone());

        isPIP = (estimationMethod.getModel().getName().find("PIP") != string::npos);

        if (isPIP) {

            tshlib::Utree *utree = nullptr;
            std::map<std::string, std::string> default_map;

            tree->setNodeName(tree->getRootId(), "root");

            UtreeBppUtils::renameInternalNodes(tree);

            utree = new tshlib::Utree();

            UtreeBppUtils::convertTree_b2u(tree, utree, tm);

            utree->addVirtualRootNode();

            // Once the tree has the root, then map it as well
            tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode->getVnode_id()));

            tl = new bpp::UnifiedTSHomogeneousTreeLikelihood_PIP(*tree, *estimationMethod.getData(), model.get(), rdist.get(), utree, &tm, true,default_map, "", false, false, false);

            delete utree;

        } else {

            tl = new RHomogeneousTreeLikelihood_Generic(*tree, *estimationMethod.getData(), model.get(), rdist.get(), true, verbose > 1);

        }

        tl->initialize();

        parameters = tl->getParameters();

        if (!optimizeBrLen) {
            parameters.deleteParameters(tl->getBranchLengthsParameters().getParameterNames());
        }

        parameters.deleteParameters(parametersToIgnore.getParameterNames());


        if (isPIP) {


            this->backupListener = NULL;
            this->nstep = 0;
            this->reparam = false;
            this->optVerbose = false;
            this->optMethodDeriv = OptimizationTools::OPTIMIZATION_BFGS;
            this->optMethodModel = OptimizationTools::OPTIMIZATION_BFGS;

            Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(
                    dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),
                    parameters);

        } else {

            OptimizationTools::optimizeNumericalParameters(tl,
                    parameters, NULL, 0, this->tolerance, this->nbEvalMax,
                                                           this->messageHandler, this->profiler,this->verbose > 0 ? this->verbose - 1 : 0);

        }

    }

    TreeTemplate<Node> *Optimizators::buildDistanceTreeGeneric(UnifiedDistanceEstimation &estimationMethod,
                                                               AgglomerativeDistanceMethod &reconstructionMethod,
                                                               bpp::ParameterList &parametersToIgnore,
                                                               bool optimizeBrLen,
                                                               const std::string &param) throw(Exception) {

        TreeTemplate<Node> *tree = nullptr;
        TreeTemplate<Node> *previousTree = nullptr;
        DistanceMatrix *matrix = nullptr;
        bpp::AbstractHomogeneousTreeLikelihood *tl = nullptr;
        UtreeBppUtils::treemap tm;
        bool test = true;

        estimationMethod.resetAdditionalParameters();
        estimationMethod.setVerbose(verbose);

        if (param == DISTANCEMETHOD_PAIRWISE) {

            ParameterList tmp = estimationMethod.getModel().getIndependentParameters();

            tmp.addParameters(estimationMethod.getRateDistribution().getIndependentParameters());

            tmp.deleteParameters(parametersToIgnore.getParameterNames());

            estimationMethod.setAdditionalParameters(tmp);
        }

        while (test) {

            // Compute matrice:

            if (this->verbose > 0){
                ApplicationTools::displayTask("Estimating distance matrix", true);
            }

            estimationMethod.computeMatrix();

            matrix = estimationMethod.getMatrix();

            if (this->verbose > 0){
                ApplicationTools::displayTaskDone();
            }

            if (matrix->size() == 2) {
                tree = this->computeTree2leaves(matrix);
                break;
            }

            if (this->verbose > 0){
                ApplicationTools::displayTask("Building tree");
            }

            reconstructionMethod.setDistanceMatrix(*matrix);

            reconstructionMethod.computeTree();

            previousTree = tree;

            tree = dynamic_cast<TreeTemplate<Node> *>(reconstructionMethod.getTree());

            if (verbose > 0){
                ApplicationTools::displayTaskDone();
            }

            if (previousTree) {

                int rf = TreeTools::robinsonFouldsDistance(*previousTree, *tree, false);

                if(verbose > 0){
                    ApplicationTools::displayResult("Topo. distance with previous iteration", TextTools::toString(rf));
                }

                test = ( rf == 0);

                delete previousTree;
            }

            if (param != DISTANCEMETHOD_ITERATIONS){
                break;  // Ends here.
            }

            // Now, re-estimate parameters:
            this->reEstimateParameters(tl,estimationMethod,tree,tm,optimizeBrLen,parametersToIgnore);

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

            delete tl;
            delete matrix;
        }

        return tree;
    }

    TreeTemplate<Node>* Optimizators::computeTree2leaves(DistanceMatrix *dmatrix){

        TreeTemplate<Node> *tree = nullptr;

        Node *n1 = nullptr;
        Node *n2 = nullptr;
        Node *n3 = nullptr;

        //Special case, there is only one possible tree:
        n1 = new Node(0);
        n2 = new Node(1, dmatrix->getName(0));

        n2->setDistanceToFather((*dmatrix)(0, 0) / 2.);

        n3 = new Node(2, dmatrix->getName(1));

        n3->setDistanceToFather((*dmatrix)(0, 0) / 2.);
        n1->addSon(n2);
        n1->addSon(n3);

        tree = new TreeTemplate<Node>(n1);

        delete n1;
        delete n2;
        delete n3;

        return tree;
    }

    TreeTemplate<Node> *Optimizators::buildDistanceTreeGenericFromDistanceMatrix(DistanceMatrix *dmatrix,
                                                                                 AgglomerativeDistanceMethod &reconstructionMethod) {

        TreeTemplate<Node> *tree = nullptr;
        DistanceMatrix *matrix = nullptr;

        // Compute matrice:
        if (this->verbose > 0){
            ApplicationTools::displayTask("Importing distance matrix", true);
        }

        matrix = dmatrix;

        if (this->verbose > 0){
            ApplicationTools::displayTaskDone();
        }

        if (this->verbose > 0){
            ApplicationTools::displayTask("Building tree");
        }

        // Compute tree:
        if (matrix->size() == 2) {

            tree =  this->computeTree2leaves(dmatrix);

        }else{

            reconstructionMethod.setDistanceMatrix(*matrix);

            reconstructionMethod.computeTree();

            tree = dynamic_cast<TreeTemplate<Node> *>(reconstructionMethod.getTree());
        }

        if (this->verbose > 0){
            ApplicationTools::displayTaskDone();
        }

        delete matrix;

        return tree;
    }


    unsigned int Optimizators::optimizeNumericalParametersUsingNumericalDerivatives(DiscreteRatesAcrossSitesTreeLikelihood *tl,
            const ParameterList &parameters)throw(Exception) {

        MetaOptimizerInfos *desc = nullptr;
        DerivableSecondOrder *f = nullptr;
        MetaOptimizer *poptimizer = nullptr;
        NaNListener *nanListener = nullptr;
        AbstractNumericalDerivative *fnum  = nullptr;
        AbstractNumericalDerivative *fnum5 = nullptr;

        ParameterList plsm;
        ParameterList plrd;

        unsigned int nb = 0;

        unique_ptr<DerivableSecondOrder> frep; // Shall we reparametrize the function to remove constraints?
        ParameterList pl;

        f = tl;
        pl = parameters;

        if (this->reparam) {
            frep.reset(new ReparametrizationDerivableSecondOrderWrapper(f, parameters));
            f = frep.get();

            // Reset parameters to remove constraints:
            pl = f->getParameters().subList(parameters.getParameterNames());
        }

        // ///////////////
        // Build optimizer:

        // Branch lengths (via numerical derivatives)

        desc = new MetaOptimizerInfos();
        fnum = new ThreePointsNumericalDerivative(f);
        fnum5 = new FivePointsNumericalDerivative(f);

        if (optMethodDeriv == OPTIMIZATION_GRADIENT) {
            desc->addOptimizer("Branch length parameters", new ConjugateGradientMultiDimensions(fnum5),tl->getBranchLengthsParameters().getParameterNames(), 2,MetaOptimizerInfos::IT_TYPE_FULL);
        }else if (optMethodDeriv == OPTIMIZATION_NEWTON){
            desc->addOptimizer("Branch length parameters", new PseudoNewtonOptimizer(fnum), tl->getBranchLengthsParameters().getParameterNames(), 2,MetaOptimizerInfos::IT_TYPE_FULL);
        }else if (optMethodDeriv == OPTIMIZATION_BFGS){
            desc->addOptimizer("Branch length parameters", new BfgsMultiDimensions(fnum5), tl->getBranchLengthsParameters().getParameterNames(), 2,MetaOptimizerInfos::IT_TYPE_FULL);
        }else if (optMethodDeriv == OPTIMIZATION_BRENT){
            desc->addOptimizer("Branch length parameters", new SimpleMultiDimensions(fnum5), tl->getBranchLengthsParameters().getParameterNames(), 2,MetaOptimizerInfos::IT_TYPE_FULL);
        }else{
            throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodDeriv);
        }

        // Other parameters
        if (optMethodModel == OPTIMIZATION_BRENT) {

            plsm = parameters.getCommonParametersWith(tl->getSubstitutionModelParameters());

            desc->addOptimizer("Substitution model parameter", new SimpleMultiDimensions(fnum5), plsm.getParameterNames(), 0,MetaOptimizerInfos::IT_TYPE_STEP);

            plrd = parameters.getCommonParametersWith(tl->getRateDistributionParameters());

            desc->addOptimizer("Rate distribution parameter", new SimpleMultiDimensions(fnum5), plrd.getParameterNames(), 0,MetaOptimizerInfos::IT_TYPE_STEP);

            poptimizer = new MetaOptimizer(fnum5, desc, nstep);

        } else if (optMethodModel == OPTIMIZATION_BFGS) {

            vector<std::string> vNameDer;
            vector<std::string> vNameDer2;

            plsm = parameters.getCommonParametersWith(tl->getSubstitutionModelParameters());

            vNameDer = plsm.getParameterNames();

            plrd = parameters.getCommonParametersWith(tl->getRateDistributionParameters());

            vNameDer2 = plrd.getParameterNames();

            vNameDer.insert(vNameDer.begin(), vNameDer2.begin(), vNameDer2.end());

            fnum->setParametersToDerivate(vNameDer);

            desc->addOptimizer("Rate & model distribution parameters", new BfgsMultiDimensions(fnum), vNameDer, 1, MetaOptimizerInfos::IT_TYPE_FULL);

            poptimizer = new MetaOptimizer(fnum, desc, nstep);

        } else{
            throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodModel);
        }

        poptimizer->setVerbose(verbose);
        poptimizer->setProfiler(profiler);
        poptimizer->setMessageHandler(messageHandler);
        poptimizer->setMaximumNumberOfEvaluations(this->nbEvalMax);
        poptimizer->getStopCondition()->setTolerance(tolerance);
        poptimizer->getDefaultStopCondition()->setTolerance(tolerance);

        // Optimize TreeLikelihood function:
        poptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);

        nanListener = new NaNListener(poptimizer, tl);

        poptimizer->addOptimizationListener(nanListener);

        if (this->backupListener){
            poptimizer->addOptimizationListener(this->backupListener);
        }

        poptimizer->init(pl);

        poptimizer->optimize();

        if (this->verbose > 0){
            ApplicationTools::displayMessage("\n");
        }

        nb = poptimizer->getNumberOfEvaluations();

        delete poptimizer;
        delete desc;
        delete nanListener;

        return nb;
    }


}