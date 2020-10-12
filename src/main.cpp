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
 * @file main.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 20 12 2017
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
 * @see For more information visit: http://bit.ly/castor-tree
 *
 */

#include <iostream>
#include <fstream>
/*
* From Core:
*/
#include <Bpp/Io/OutputStream.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Text/KeyvalTools.h>

/*
* From SeqLib:
*/
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include <Bpp/Seq/AlphabetIndex/DefaultNucleotideScore.h>

/*
* From PhylLib:
*/
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Distance/PGMA.h>
#include <Bpp/Phyl/OptimizationTools.h>

/*
* From GLOG:
*/
#include <glog/logging.h>

/*
* From TSHLib:
*/
#include <TreeRearrangment.hpp>

/*
* From Boost:
*/
#include <boost/algorithm/string/split.hpp>

using namespace tshlib;

/*
* From miniJati:
*/
#include "Version.hpp"
#include "Utils.hpp"
#include "PIP.hpp"
#include "ExtendedAlphabet.hpp"
#include "RHomogeneousTreeLikelihood_PIP.hpp"
#include "RHomogeneousTreeLikelihood_Generic.hpp"
#include "CastorApplication.hpp"
#include "UnifiedTSHTopologySearch.hpp"
#include "Optimizators.hpp"
#include "SupportMeasures.hpp"
#include "UnifiedDistanceEstimation.hpp"

#include "pPIP.hpp"
#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"
#include "CompositePIPnode.hpp"

#include "inference_indel_rates.hpp"

/*
* From GSL:
*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "DistanceFactory.hpp"
#include "DistanceFactoryAngle.hpp"
#include "DistanceFactoryAlign.hpp"
#include "DistanceFactoryPrealigned.hpp"

int main(int argc, char *argv[]) {

    //FLAGS_logtostderr = 1;
    FLAGS_log_dir = ".";
    google::InitGoogleLogging(software::name.c_str());
    google::InstallFailureSignalHandler();

    try {

        bpp::CastorApplication castorapp(argc,
                                         argv,
                                         std::string(software::name + " " + software::version),
                                         std::string(software::releasegitbranch + " " + software::releasegitref),
                                         std::string(software::releasedate + ", " + software::releasetime));


        //////////////////////////////////////////////
        // START

        bpp::ApplicationTools::displayMessage("\n[Starting the application]");

        castorapp.start(argc);

        bpp::ApplicationTools::displayResult("Random seed set to", castorapp.getSeed());
        bpp::ApplicationTools::displayResult("Log files location", std::string("current execution path"));

        //////////////////////////////////////////////
        // CLI ARGUMENTS

        bpp::ApplicationTools::displayMessage("\n[Getting CLI arguments]");

        castorapp.getCLIarguments();

        //////////////////////////////////////////////
        // ALPHABET

        bpp::ApplicationTools::displayMessage("\n[Getting alphabet]");

        castorapp.getAlphabet();

        bpp::ApplicationTools::displayResult("Alphabet", TextTools::toString(castorapp.alphabetNoGaps->getAlphabetType()));
        bpp::ApplicationTools::displayBooleanResult("Allow gaps as extra character", castorapp.PAR_model_indels);
        bpp::ApplicationTools::displayResult("Genetic Code", castorapp.codeDesc);

        DLOG(INFO) << "alphabet:  " << castorapp.PAR_Alphabet << " | gap-extention " << (int) castorapp.PAR_model_indels;

        //////////////////////////////////////////////
        // DATA
        bpp::ApplicationTools::displayMessage("\n[Preparing input data]");

        castorapp.getData();

        bpp::ApplicationTools::displayBooleanResult("Aligned sequences", !castorapp.PAR_alignment);
        if (castorapp.PAR_alignment) {
            bpp::ApplicationTools::displayResult("Number of sequences",TextTools::toString(castorapp.sequences->getNumberOfSequences()));
        } else {
            bpp::ApplicationTools::displayResult("Number of sequences",TextTools::toString(castorapp.sites->getNumberOfSequences()));
            bpp::ApplicationTools::displayResult("Number of sites", TextTools::toString(castorapp.sites->getNumberOfSites()));
        }

        /////////////////////////////////////////
        // INITIAL TREE
        bpp::ApplicationTools::displayMessage("\n[Preparing initial tree]");

        /*
        bpp::Tree *tree = nullptr;
        string initTreeOpt = ApplicationTools::getStringParameter("init.tree", castorapp.getParams(), "user", "", false,
                                                                  1);
        ApplicationTools::displayResult("Initial tree", initTreeOpt);

        if (initTreeOpt == "user") {

            tree = PhylogeneticsApplicationTools::getTree(castorapp.getParams());
            DLOG(INFO) << "[Input tree parser] Number of leaves" << tree->getNumberOfLeaves();

        } else if (initTreeOpt == "random") {

            vector<std::string> names = castorapp.sites->getSequencesNames();
            tree = TreeTemplateTools::getRandomTree(names);
            tree->setBranchLengths(1.);

        } else if (initTreeOpt == "distance") {

            bpp::DistanceMatrix *distances;

            std::string PAR_distance_method = ApplicationTools::getStringParameter("init.distance.method",castorapp.getParams(), "nj");
            ApplicationTools::displayResult("Initial tree reconstruction method", PAR_distance_method);

            AgglomerativeDistanceMethod *distMethod = nullptr;
            std::string token = PAR_distance_method.substr(0, PAR_distance_method.find("-"));
            if (token == "wpgma") {
                auto *wpgma = new PGMA(true);
                distMethod = wpgma;
            } else if (token == "upgma") {
                auto *upgma = new PGMA(false);
                distMethod = upgma;
            } else if (token == "nj") {
                auto *nj = new NeighborJoining();
                nj->outputPositiveLengths(true);
                distMethod = nj;
            } else if (token == "bionj") {
                auto *bionj = new BioNJ();
                bionj->outputPositiveLengths(true);
                distMethod = bionj;
            } else if (token == "distmatrix") {
                // Use a distance matrix provided by the user
                ApplicationTools::displayResult("Initial tree method", std::string("LZ compression"));
                std::string PAR_distance_matrix;
                try {
                    PAR_distance_matrix = ApplicationTools::getAFilePath("init.distance.matrix.file",
                                                                         castorapp.getParams(), true, true, "", false,
                                                                         "",
                                                                         0);
                } catch (bpp::Exception &e) {
                    LOG(FATAL) << "Error when reading distance matrix file: " << e.message();
                }

                DLOG(INFO) << "initial tree method from LZ compression from matrix file" << PAR_distance_matrix;
                distances = InputUtils::parseDistanceMatrix(PAR_distance_matrix);
                bpp::BioNJ bionj(*distances, true, true, false);
                tree = bionj.getTree();
            } else if (token == "infere_distance_matrix") {

                int ALPHABET_DIM=0;
                int K=0;
                bool mldist_flag = true;
                bool mldist_gap_flag = false;
                double cutoff_dist = 1.0;
                double indel_rate = 1.0;

                if (castorapp.PAR_Alphabet.find("DNA") != std::string::npos) {
                    ALPHABET_DIM = 4;
                    K = 6;
                } else if (castorapp.PAR_Alphabet.find("Protein") != std::string::npos) {
                    ALPHABET_DIM = 20;
                    K = 2;
                } else if (castorapp.PAR_Alphabet.find("Codon") != std::string::npos) {
                    ALPHABET_DIM = 60;
                    K = 2;
                }

                if (!castorapp.sequences) {
                    bpp::Fasta seqReader;
                    castorapp.sequences = seqReader.readSequences(castorapp.PAR_input_sequences, castorapp.alphabet);
                }
                DistanceFactoryPrographMSA::DistanceFactory *dist_factory_angle = new DistanceFactoryPrographMSA::DistanceFactoryAngle(
                        ALPHABET_DIM, K);

                DistanceFactoryPrographMSA::DistanceMatrix dist_ml = dist_factory_angle->computePwDistances(castorapp.sequences,
                                                                                                            ALPHABET_DIM,
                                                                                                            K,
                                                                                                            mldist_flag,
                                                                                                            mldist_gap_flag,
                                                                                                            cutoff_dist,
                                                                                                            indel_rate);

                bpp::DistanceMatrix *dist_ = new DistanceMatrix(castorapp.sequences->getSequencesNames());

                for (int iii = 0; iii < castorapp.sequences->getNumberOfSequences(); iii++) {
                    for (int jjj = 0; jjj < castorapp.sequences->getNumberOfSequences(); jjj++) {
                        (*dist_)(iii, jjj) = (abs(dist_ml.distances(iii, jjj)) < DISTCUTOFF ? 0.0 : abs(
                                dist_ml.distances(iii, jjj)));
                    }
                }

                bpp::DistanceMethod *distMethod = nullptr;
                auto *bionj = new BioNJ(true, true, false);
                bionj->outputPositiveLengths(true);
                distMethod = bionj;
                distMethod->setDistanceMatrix(*dist_);
                distMethod->computeTree();
                tree = distMethod->getTree();

                delete dist_factory_angle;
                delete dist_;
                delete bionj;

            } else throw Exception("Unknown tree reconstruction method.");

            //----------------------------------------------------------------
            //m@x
            if (true) {

                // Compute bioNJ tree using the GTR model
                map<std::string, std::string> parmap;

                VectorSiteContainer *allSites;
                VectorSiteContainer *sitesDistMethod;
                bpp::Alphabet *alphabetDistMethod;

                if (castorapp.PAR_model_indels) {
                    alphabetDistMethod = castorapp.alphabet;
                    parmap["model"] = castorapp.modelMap["model"];
                } else {
                    alphabetDistMethod = castorapp.alphabetNoGaps;
                    parmap["model"] = castorapp.getParams()["model"];
                }



                allSites = SequenceApplicationTools::getSiteContainer(alphabetDistMethod, castorapp.getParams());
                sitesDistMethod = SequenceApplicationTools::getSitesToAnalyse(*allSites, castorapp.getParams());
                delete allSites;


                bpp::SubstitutionModel *smodel;
                if (castorapp.PAR_model_indels) {

                    // Instantiation of the canonical substitution model
                    if (castorapp.PAR_Alphabet.find("Codon") != std::string::npos ||
                            castorapp.PAR_Alphabet.find("Protein") != std::string::npos) {
                        smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(castorapp.alphabetNoGaps, castorapp.gCode.get(),
                                                                                          sitesDistMethod, castorapp.modelMap, "",
                                                                                          true,
                                                                                          false, 0);
                    } else {
                        smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabetDistMethod,
                                                                                          castorapp.gCode.get(), sitesDistMethod,
                                                                                          castorapp.modelMap,
                                                                                          "", true,
                                                                                          false, 0);
                    }


                    if(castorapp.modelMap.find("lambda") == castorapp.modelMap.end() || castorapp.modelMap.find("mu") == castorapp.modelMap.end()){


                        //==================================================================================
                        //==================================================================================
                        //==================================================================================
                        //==================================================================================
                        //==================================================================================
                        // m@x:: new code
                        if(!tree){

                            double lambda_tmp = 10.0;
                            double mu_tmp = 0.1;

                            bpp::SubstitutionModel *smodel_tmp;
                            bpp::SubstitutionModel *smodel_copy = smodel->clone();

                            VectorSiteContainer *sitesDistMethod_tmp = sitesDistMethod->clone();
                            bpp::Alphabet *alphabetDistMethod_tmp = alphabetDistMethod->clone();

                            // Instatiate the corrisponding PIP model given the alphabet
                            if (castorapp.PAR_Alphabet.find("DNA") != std::string::npos &&
                                    castorapp.PAR_Alphabet.find("Codon") == std::string::npos) {
                                smodel_tmp = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabetDistMethod_tmp), smodel_copy,
                                                     *sitesDistMethod_tmp, lambda_tmp, mu_tmp, false);
                            } else if (castorapp.PAR_Alphabet.find("Protein") != std::string::npos) {
                                smodel_tmp = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabetDistMethod_tmp), smodel_copy,
                                                    *sitesDistMethod_tmp, lambda_tmp, mu_tmp, false);
                            } else if (castorapp.PAR_Alphabet.find("Codon") != std::string::npos) {
                                smodel_tmp = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabetDistMethod_tmp), castorapp.gCode.get(),
                                                           smodel_copy, *sitesDistMethod_tmp,lambda_tmp,mu_tmp, false);
                                ApplicationTools::displayWarning(
                                        "Codon models are experimental in the current version... use with caution!");
                                DLOG(WARNING) << "CODONS activated but the program is not fully tested under these settings!";
                            }

                            //Initialize model to compute the distance tree
                            //TransitionModel *dmodel = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod, gCode.get(), sitesDistMethod, parmap);
                            TransitionModel *dmodel_tmp;
                            // Get transition model from substitution model
                            if (!castorapp.PAR_model_indels) {
                                dmodel_tmp = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod_tmp, castorapp.gCode.get(),
                                                                                               sitesDistMethod_tmp, parmap);
                            } else {
                                unique_ptr<TransitionModel> test;
                                test.reset(smodel_tmp);
                                dmodel_tmp = test.release();
                            }

                            // Add a ASRV distribution
                            DiscreteDistribution *rDist_tmp = nullptr;
                            if (dmodel_tmp->getNumberOfStates() > dmodel_tmp->getAlphabet()->getSize()) {
                                //Markov-modulated Markov model!
                                rDist_tmp = new ConstantRateDistribution();
                            } else {
                                rDist_tmp = PhylogeneticsApplicationTools::getRateDistribution(castorapp.getParams());
                            }

                            // Remove gap characters since we are roughly estimating the initial topology
                            if (!castorapp.PAR_model_indels) {
                                bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sitesDistMethod_tmp);
                            }

                            UnifiedDistanceEstimation distEstimation_tmp(dmodel_tmp, rDist_tmp, sitesDistMethod_tmp, 1, false);

                            // old LORENZO version



                            bpp::DistanceMatrix * dm_tmp = bpp::SiteContainerTools::computeSimilarityMatrix(*sitesDistMethod_tmp,true,"no full gap",true);



                            bpp::DistanceMethod *distMethod_tmp = nullptr;
                            auto *bionj_tmp = new BioNJ(true, true, false);
                            bionj_tmp->outputPositiveLengths(true);
                            distMethod_tmp = bionj_tmp;
                            distMethod_tmp->setDistanceMatrix((*dm_tmp));
                            distMethod_tmp->computeTree();
                            tree = distMethod_tmp->getTree();



                        }
                        //==================================================================================
                        //==================================================================================
                        //==================================================================================
                        //==================================================================================
                        //==================================================================================

                        inference_indel_rates::infere_indel_rates_from_sequences(castorapp.PAR_input_sequences,
                                                                                 castorapp.PAR_Alphabet,
                                                                                 castorapp.PAR_alignment,
                                                                                 castorapp.PAR_model_indels,
                                                                                 castorapp.getParams(),
                                                                                 tree,
                                                                                 castorapp.lambda,
                                                                                 castorapp.mu,
                                                                                 castorapp.gCode.get(),
                                                                                 castorapp.modelMap);


                    }else{
                        castorapp.lambda = std::stod(castorapp.modelMap["lambda"]);
                        castorapp.mu = std::stod(castorapp.modelMap["mu"]);
                    }

                    // Instatiate the corrisponding PIP model given the alphabet
                    if (castorapp.PAR_Alphabet.find("DNA") != std::string::npos &&
                            castorapp.PAR_Alphabet.find("Codon") == std::string::npos) {
                        smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabetDistMethod), smodel,
                                             *sitesDistMethod, castorapp.lambda, castorapp.mu, false);
                    } else if (castorapp.PAR_Alphabet.find("Protein") != std::string::npos) {
                        smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabetDistMethod), smodel,
                                            *sitesDistMethod, castorapp.lambda, castorapp.mu, false);
                    } else if (castorapp.PAR_Alphabet.find("Codon") != std::string::npos) {
                        smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabetDistMethod), castorapp.gCode.get(),
                                               smodel, *sitesDistMethod,
                                               castorapp.lambda,
                                               castorapp.mu, false);
                        ApplicationTools::displayWarning(
                                "Codon models are experimental in the current version... use with caution!");
                        DLOG(WARNING) << "CODONS activated but the program is not fully tested under these settings!";
                    }

                }

                //Initialize model to compute the distance tree
                //TransitionModel *dmodel = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod, gCode.get(), sitesDistMethod, parmap);
                TransitionModel *dmodel;
                // Get transition model from substitution model
                if (!castorapp.PAR_model_indels) {
                    dmodel = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod, castorapp.gCode.get(),
                                                                               sitesDistMethod, parmap);
                } else {
                    unique_ptr<TransitionModel> test;
                    test.reset(smodel);
                    dmodel = test.release();
                }

                // Add a ASRV distribution
                DiscreteDistribution *rDist = nullptr;
                if (dmodel->getNumberOfStates() > dmodel->getAlphabet()->getSize()) {
                    //Markov-modulated Markov model!
                    rDist = new ConstantRateDistribution();
                } else {
                    rDist = PhylogeneticsApplicationTools::getRateDistribution(castorapp.getParams());
                }

                // Remove gap characters since we are roughly estimating the initial topology
                if (!castorapp.PAR_model_indels) {
                    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sitesDistMethod);
                }

                UnifiedDistanceEstimation distEstimation(dmodel, rDist, sitesDistMethod, 1, false);

                if (PAR_distance_method.find("-ml") != std::string::npos) {

                    std::string PAR_optim_distance = ApplicationTools::getStringParameter(
                            "init.distance.optimization.method", castorapp.getParams(),
                            "init");
                    ApplicationTools::displayResult("Initial tree model parameters estimation method",
                                                    PAR_optim_distance);
                    if (PAR_optim_distance == "init") PAR_optim_distance = Optimizators::DISTANCEMETHOD_INIT;
                    else if (PAR_optim_distance == "pairwise")
                        PAR_optim_distance = Optimizators::DISTANCEMETHOD_PAIRWISE;
                    else if (PAR_optim_distance == "iterations")
                        PAR_optim_distance = Optimizators::DISTANCEMETHOD_ITERATIONS;
                    else throw Exception("Unknown parameter estimation procedure '" + PAR_optim_distance + "'.");

                    // Optimisation method verbosity
                    auto optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose",
                                                                                   castorapp.getParams(), 2);
                    string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", castorapp.getParams(),
                                                                   false, false);
                    auto *messenger = (mhPath == "none") ? nullptr : (mhPath == "std") ? ApplicationTools::message.get()
                                                                                       : new StlOutputStream(
                                    new ofstream(mhPath.c_str(), ios::out));
                    ApplicationTools::displayResult("Initial tree optimization handler", mhPath);

                    // Optimisation method profiler
                    string prPath = ApplicationTools::getAFilePath("optimization.profiler", castorapp.getParams(), false,
                                                                   false);
                    auto *profiler = (prPath == "none") ? nullptr : (prPath == "std") ? ApplicationTools::message.get()
                                                                                      : new StlOutputStream(
                                    new ofstream(prPath.c_str(), ios::out));
                    if (profiler) profiler->setPrecision(20);
                    ApplicationTools::displayResult("Initial tree optimization profiler", prPath);

                    // Should I ignore some parameters?
                    ParameterList allParameters = dmodel->getParameters();
                    allParameters.addParameters(rDist->getParameters());

                    ParameterList parametersToIgnore;
                    string paramListDesc = ApplicationTools::getStringParameter(
                            "init.distance.optimization.ignore_parameter", castorapp.getParams(),
                            "", "", true, false);
                    bool ignoreBrLen = false;
                    StringTokenizer st(paramListDesc, ",");

                    while (st.hasMoreToken()) {
                        try {
                            string param = st.nextToken();
                            if (param == "BrLen")
                                ignoreBrLen = true;
                            else {
                                if (allParameters.hasParameter(param)) {
                                    Parameter *p = &allParameters.getParameter(param);
                                    parametersToIgnore.addParameter(*p);
                                } else ApplicationTools::displayWarning("Parameter '" + param + "' not found.");
                            }
                        } catch (ParameterNotFoundException &pnfe) {
                            ApplicationTools::displayError(
                                    "Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
                        }
                    }

                    auto nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval",
                                                                                  castorapp.getParams(), 1000000);
                    ApplicationTools::displayResult("Initial tree optimization | max # ML evaluations",
                                                    TextTools::toString(nbEvalMax));

                    double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance",
                                                                            castorapp.getParams(), .000001);
                    ApplicationTools::displayResult("Initial tree optimization | Tolerance",
                                                    TextTools::toString(tolerance));

                    //Here it is:
                    tree = Optimizators::buildDistanceTreeGeneric(distEstimation, *distMethod, parametersToIgnore,
                                                                  !ignoreBrLen, PAR_optim_distance,
                                                                  tolerance, nbEvalMax, profiler, messenger,
                                                                  optVerbose);

                } else {
                    // Fast but rough estimate of the initial tree topology (distance based without optimisation -ML)

                    if(!tree) {
                        // LORENZO version
                        distEstimation.computeMatrix();
                        DistanceMatrix *dm = distEstimation.getMatrix();
                        distMethod->setDistanceMatrix((*dm));
                        distMethod->computeTree();
                        tree = distMethod->getTree();
                    }

                }


                delete sitesDistMethod;
                delete distMethod;

            } //m@x
            //----------------------------------------------------------------

        } else throw Exception("Unknown init tree method.");


        // If the tree has multifurcation, then resolve it with midpoint rooting
        auto ttree_ = new TreeTemplate<Node>(*tree);
        if (ttree_->getRootNode()->getNumberOfSons() > 2) {
            TreeTemplateTools::midRoot(*(ttree_), TreeTemplateTools::MIDROOT_VARIANCE, false);
            tree = ttree_;
        }

        // Rename internal nodes with standard Vxx * where xx is a progressive number
        tree->setNodeName(tree->getRootId(), "root");
        UtreeBppUtils::renameInternalNodes(tree);

        // Try to write the current tree to file. This will be overwritten by the optimized tree,
        // but allow to check file existence before running optimization!
        PhylogeneticsApplicationTools::writeTree(*tree, castorapp.getParams());

        // Setting branch lengths?
        string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", castorapp.getParams(), "Input",
                                                                      "", true, 1);
        string cmdName;
        map<string, string> cmdArgs;
        KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);
        if (cmdName == "Input") {
            // Is the root has to be moved to the midpoint position along the branch that contains it ? If no, do nothing!
            bool midPointRootBrLengths = ApplicationTools::getBooleanParameter("midpoint_root_branch", cmdArgs, false,
                                                                               "", true, 2);
            if (midPointRootBrLengths)
                TreeTools::constrainedMidPointRooting(*tree);
        } else if (cmdName == "Equal") {
            double value = ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, 2);
            if (value <= 0)
                throw Exception("Value for branch length must be superior to 0");
            ApplicationTools::displayResult("Branch lengths set to", value);
            tree->setBranchLengths(value);
        } else if (cmdName == "Clock") {
            TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
        } else if (cmdName == "Grafen") {
            string grafenHeight = ApplicationTools::getStringParameter("height", cmdArgs, "input", "", true, 2);
            double h;
            if (grafenHeight == "input") {
                h = TreeTools::getHeight(*tree, tree->getRootId());
            } else {
                h = TextTools::toDouble(grafenHeight);
                if (h <= 0) throw Exception("Height must be positive in Grafen's method.");
            }
            ApplicationTools::displayResult("Total height", TextTools::toString(h));

            double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, 2);
            ApplicationTools::displayResult("Grafen's rho", rho);
            TreeTools::computeBranchLengthsGrafen(*tree, rho);
            double nh = TreeTools::getHeight(*tree, tree->getRootId());
            tree->scaleTree(h / nh);
        } else throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");
        ApplicationTools::displayResult("Branch lengths", cmdName);

        DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(tree);

        // Convert the bpp into utree for tree-search engine
        auto utree = new Utree();
        UtreeBppUtils::treemap tm;
        UtreeBppUtils::convertTree_b2u(tree, utree, tm);
        if (castorapp.PAR_alignment) {
            UtreeBppUtils::associateNode2Alignment(castorapp.sequences, utree);
        } else {
            UtreeBppUtils::associateNode2Alignment(castorapp.sites, utree);
        }

        DLOG(INFO) << "Bidirectional map size: " << tm.size();
        DLOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

        utree->addVirtualRootNode();
        // Once the tree has the root, then map it as well
        //tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode));

        tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode->getVnode_id()));

        */

        castorapp.getTree();


        ApplicationTools::displayResult("Initial tree total length", TextTools::toString(castorapp.tree->getTotalLength(), 6));
        /////////////////////////
        // SUBSTITUTION MODEL

        bpp::ApplicationTools::displayMessage("\n[Setting up substitution model]");
        
        castorapp.getSubstitutionModel(castorapp.tree);

        /////////////////////////
        // AMONG-SITE-RATE-VARIATION

        bpp::ApplicationTools::displayMessage("\n[Getting ASRV distribution]");

        castorapp.getASRV();

        /////////////////////////
        // ALIGN SEQUENCES

        if (castorapp.PAR_alignment) {

            bpp::ApplicationTools::displayMessage("\n[Computing the alignment]");

            castorapp.computeMSA(castorapp.tree,castorapp.utree,castorapp.tm);

            bpp::ApplicationTools::displayResult("Aligner optimised for:", castorapp.PAR_alignment_version);

            bpp::ApplicationTools::displayBooleanResult("Stochastic backtracking active", castorapp.PAR_alignment_sbsolutions > 1);
            if (castorapp.PAR_alignment_sbsolutions > 1) {
                bpp::ApplicationTools::displayResult("Number of stochastic solutions:",TextTools::toString(castorapp.PAR_alignment_sbsolutions));
            }

            bpp::ApplicationTools::displayResult("Alignment log likelihood", TextTools::toString(castorapp.proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getScore(), 15));
        }

        /////////////////////////
        // HOMOGENEOUS MODELING - initialization likelihood functions

        castorapp.initLkFun(castorapp.tree,castorapp.utree,castorapp.tm);

        /////////////////////////
        // PARAMETER SANITY CHECK

        ApplicationTools::displayMessage("\n[Parameter sanity check]");

        castorapp.parameterSanityCheck();

        /////////////////////////
        // OPTIMISE PARAMETERS

        bpp::ApplicationTools::displayMessage("\n[Executing numerical parameters and topology optimization]");

        castorapp.optimizeParameters();

        /////////////////////////
        // OUTPUT

        bpp::ApplicationTools::displayMessage("\n[Printing parameters]");

        castorapp.output();

        /////////////////////////
        // BOOTSTRAPING

        castorapp.bootstrapping(castorapp.utree,castorapp.tm);

        /////////////////////////
        // CLEAR OBJECTS

        // Delete objects and free memory
        delete castorapp.alphabet;
        delete castorapp.alphabetNoGaps;
        delete castorapp.sites;
        delete castorapp.rDist;
        delete castorapp.tl;
        delete castorapp.sequences;
        delete castorapp.tree;

        /////////////////////////
        // CLOSE APPLICATION
        castorapp.done();
        google::ShutdownGoogleLogging();
        exit(EXIT_SUCCESS);

    } catch (exception &e) {
        cout << e.what() << endl;
        google::ShutdownGoogleLogging();
        exit(EXIT_SUCCESS);
    }
}
