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

        if (argc < 2) {
            castorapp.help();
            exit(0);
        } else {
            castorapp.banner();
            castorapp.startTimer();
        };

        bpp::ApplicationTools::displayResult("Random seed set to", castorapp.getSeed());
        ApplicationTools::displayResult("Log files location", std::string("current execution path"));

        //////////////////////////////////////////////
        // CLI ARGUMENTS
        std::string PAR_model_substitution = ApplicationTools::getStringParameter("model", castorapp.getParams(), "JC69",
                                                                                  "", true, true);


        bool PAR_alignment = ApplicationTools::getBooleanParameter("alignment", castorapp.getParams(), false);

        // Split model string description and test if PIP is required
        std::string modelStringName;
        std::map<std::string, std::string> modelMap;
        KeyvalTools::parseProcedure(PAR_model_substitution, modelStringName, modelMap);
        bool PAR_model_indels = modelStringName == "PIP";


        /* ***************************************************
         * Standard workflow
         * [INPUT]
         * 1. tree + alignment => (1.1) just parse everything
         * 2. alignment  => (2.1) parse alignment => (2.2) generate tree using bioNJ
         * 3. sequences  => (3.1) parse sequences => (3.2) generate tree using bioNJ => (3.3) perform alignment
         * (it should not be supported in production)
         * 4. sequences + tree => (4.1) parse sequence => (4.2) parse tree => (4.3) perform alignment
         */

        //////////////////////////////////////////////
        // ALPHABET
        // The alphabet object contains the not-extended alphabet as requested by the user,
        // while alpha contains the extended version of the same alphabet.
        std::string PAR_Alphabet = ApplicationTools::getStringParameter("alphabet", castorapp.getParams(), "DNA", "",
                                                                        true, true);

        // Alphabet without gaps
        bpp::Alphabet *alphabetNoGaps = bpp::SequenceApplicationTools::getAlphabet(castorapp.getParams(), "", false,
                                                                                   false);

        // Genetic code
        unique_ptr<GeneticCode> gCode;

        // Codon alphabet ?
        bool codonAlphabet = false;

        // Alphabet used for all the computational steps (it can allows for gap extension)
        bpp::Alphabet *alphabet;
        if (PAR_model_indels) {

            if (PAR_Alphabet.find("DNA") != std::string::npos || PAR_Alphabet.find("Codon") != std::string::npos) {
                alphabet = new bpp::DNA_EXTENDED();
            } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
                alphabet = new bpp::ProteicAlphabet_Extended();
            }

            // This is additional to the alphabet instance
            if (PAR_Alphabet.find("Codon") != std::string::npos) {
                alphabet = new CodonAlphabet_Extended(dynamic_cast<bpp::NucleicAlphabet *>(alphabet));
                codonAlphabet = true;
            }

        } else {

            if (PAR_Alphabet.find("DNA") != std::string::npos || PAR_Alphabet.find("Codon") != std::string::npos) {
                alphabet = new bpp::DNA();
            } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
                alphabet = new bpp::ProteicAlphabet();
            } else {}

            // This is additional to the alphabet instance
            if (PAR_Alphabet.find("Codon") != std::string::npos) {
                alphabet = new CodonAlphabet(dynamic_cast<bpp::NucleicAlphabet *>(alphabet));
                codonAlphabet = true;
            }
        }

        // Alphabet used for codon models
        //CodonAlphabet *codonAlphabet = dynamic_cast<bpp::CodonAlphabet *>(alphabet);
        if (codonAlphabet) {
            std::string codeDesc = ApplicationTools::getStringParameter("genetic_code", castorapp.getParams(), "Standard",
                                                                        "", true, true);
            ApplicationTools::displayResult("Genetic Code", codeDesc);
            //if (PAR_model_indels) {
            //    gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(dynamic_cast<bpp::CodonAlphabet_Extended *>(alphabet)->getNucleicAlphabet(), codeDesc));
            //} else {
            gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(
                    dynamic_cast<bpp::CodonAlphabet *>(alphabetNoGaps)->getNucleicAlphabet(),
                    codeDesc));
            //}
        }

        ApplicationTools::displayResult("Alphabet", TextTools::toString(alphabetNoGaps->getAlphabetType()));
        ApplicationTools::displayBooleanResult("Allow gaps as extra character", PAR_model_indels);
        DLOG(INFO) << "alphabet:  " << PAR_Alphabet << " | gap-extention " << (int) PAR_model_indels;


        //////////////////////////////////////////////
        // DATA
        ApplicationTools::displayMessage("\n[Preparing input data]");
        std::string PAR_input_sequences = ApplicationTools::getAFilePath("input.sequence.file", castorapp.getParams(),
                                                                         true, true, "", false, "", 1);

        bpp::SequenceContainer *sequences = nullptr;
        bpp::SiteContainer *sites = nullptr;

        try {
            ApplicationTools::displayBooleanResult("Aligned sequences", !PAR_alignment);

            if (PAR_alignment) {

                // If the user requires the computation of an alignment, then the input file is made of unaligned sequences
                bpp::Fasta seqReader;
                sequences = seqReader.readSequences(PAR_input_sequences, alphabet);
                ApplicationTools::displayResult("Number of sequences",
                                                TextTools::toString(sequences->getNumberOfSequences()));

            } else {

                VectorSiteContainer *allSites = SequenceApplicationTools::getSiteContainer(alphabet,
                                                                                           castorapp.getParams());

                sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, castorapp.getParams(), "", true,
                                                                    !PAR_model_indels, true, 1);
                delete allSites;


                AlignmentUtils::checkAlignmentConsistency(*sites);
                ApplicationTools::displayResult("Number of sequences",
                                                TextTools::toString(sites->getNumberOfSequences()));
                ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

            }


        } catch (bpp::Exception &e) {
            LOG(FATAL) << "Error when reading sequence file due to: " << e.message();
        }


        /////////////////////////////////////////
        // INITIAL TREE
        ApplicationTools::displayMessage("\n[Preparing initial tree]");

        bpp::Tree *tree = nullptr;
        string initTreeOpt = ApplicationTools::getStringParameter("init.tree", castorapp.getParams(), "user", "", false,
                                                                  1);
        ApplicationTools::displayResult("Initial tree", initTreeOpt);

        if (initTreeOpt == "user") {

            tree = PhylogeneticsApplicationTools::getTree(castorapp.getParams());
            DLOG(INFO) << "[Input tree parser] Number of leaves" << tree->getNumberOfLeaves();

        } else if (initTreeOpt == "random") {

            vector<string> names = sites->getSequencesNames();
            tree = TreeTemplateTools::getRandomTree(names);
            tree->setBranchLengths(1.);

        } else if (initTreeOpt == "distance") {

            bpp::DistanceMatrix *distances;

            std::string PAR_distance_method = ApplicationTools::getStringParameter("init.distance.method",
                                                                                   castorapp.getParams(), "nj");
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
            } else throw Exception("Unknown tree reconstruction method.");

            //if (!PAR_alignment) {

            // Compute bioNJ tree using the GTR model
            map<std::string, std::string> parmap;

            VectorSiteContainer *allSites;
            VectorSiteContainer *sitesDistMethod;
            bpp::Alphabet *alphabetDistMethod;

            if (PAR_model_indels) {
                alphabetDistMethod = alphabet;
                parmap["model"] = modelMap["model"];
            } else {
                alphabetDistMethod = alphabetNoGaps;
                parmap["model"] = castorapp.getParams()["model"];
            }

            allSites = SequenceApplicationTools::getSiteContainer(alphabetDistMethod, castorapp.getParams());
            sitesDistMethod = SequenceApplicationTools::getSitesToAnalyse(*allSites, castorapp.getParams());
            delete allSites;

            bpp::SubstitutionModel *smodel;
            if (PAR_model_indels) {

                // Instantiation of the canonical substitution model
                if (PAR_Alphabet.find("Codon") != std::string::npos ||
                    PAR_Alphabet.find("Protein") != std::string::npos) {
                    smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabetNoGaps, gCode.get(),
                                                                                      sitesDistMethod, modelMap, "",
                                                                                      true,
                                                                                      false, 0);
                } else {
                    smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabetDistMethod,
                                                                                      gCode.get(), sitesDistMethod,
                                                                                      modelMap,
                                                                                      "", true,
                                                                                      false, 0);
                }

                double lambda = (modelMap.find("lambda") == modelMap.end()) ? 0.1 : std::stod(modelMap["lambda"]);
                double mu = (modelMap.find("mu") == modelMap.end()) ? 0.2 : std::stod(modelMap["mu"]);

                // Instatiate the corrisponding PIP model given the alphabet
                if (PAR_Alphabet.find("DNA") != std::string::npos &&
                    PAR_Alphabet.find("Codon") == std::string::npos) {
                    smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabetDistMethod), smodel,
                                         *sitesDistMethod, lambda, mu, false);
                } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
                    smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabetDistMethod), smodel,
                                        *sitesDistMethod, lambda, mu, false);
                } else if (PAR_Alphabet.find("Codon") != std::string::npos) {
                    smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabetDistMethod), gCode.get(),
                                           smodel, *sitesDistMethod,
                                           lambda,
                                           mu, false);
                    ApplicationTools::displayWarning(
                            "Codon models are experimental in the current version... use with caution!");
                    DLOG(WARNING) << "CODONS activated byt the program is not fully tested under these settings!";
                }

            }

            //Initialize model to compute the distance tree
            //TransitionModel *dmodel = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod, gCode.get(), sitesDistMethod, parmap);
            TransitionModel *dmodel;
            // Get transition model from substitution model
            if (!PAR_model_indels) {
                dmodel = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod, gCode.get(),
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
            if (!PAR_model_indels) {
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


                distEstimation.computeMatrix();
                DistanceMatrix *dm = distEstimation.getMatrix();
                distMethod->setDistanceMatrix((*dm));
                distMethod->computeTree();
                tree = distMethod->getTree();
                std::cout << std::endl;

            }

            delete sitesDistMethod;
            delete distMethod;

//            } else {
//
//                // Use a distance matrix provided by the user
//
//                ApplicationTools::displayResult("Initial tree method", std::string("LZ compression"));
//                std::string PAR_distance_matrix;
//                try {
//                    PAR_distance_matrix = ApplicationTools::getAFilePath("init.distance.matrix.file",
//                                                                         castorapp.getParams(), true, true, "", false, "",
//                                                                         0);
//                } catch (bpp::Exception &e) {
//                    LOG(FATAL) << "Error when reading distance matrix file: " << e.message();
//                }
//
//                DLOG(INFO) << "initial tree method from LZ compression from matrix file" << PAR_distance_matrix;
//                distances = InputUtils::parseDistanceMatrix(PAR_distance_matrix);
//                bpp::BioNJ bionj(*distances, true, true, false);
//                tree = bionj.getTree();
//            }

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
        if (PAR_alignment) {
            UtreeBppUtils::associateNode2Alignment(sequences, utree);
        } else {
            UtreeBppUtils::associateNode2Alignment(sites, utree);
        }

        DLOG(INFO) << "Bidirectional map size: " << tm.size();
        DLOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

        utree->addVirtualRootNode();
        // Once the tree has the root, then map it as well
        //tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode));

        tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode->getVnode_id()));

        ApplicationTools::displayResult("Initial tree total length", TextTools::toString(tree->getTotalLength(), 6));
        /////////////////////////
        // SUBSTITUTION MODEL

        ApplicationTools::displayMessage("\n[Setting up substitution model]");

        bpp::SubstitutionModel *smodel = nullptr;
        bpp::TransitionModel *model = nullptr;

        double lambda;
        double mu;
        bool estimatePIPparameters = false;

        // Instantiate a substitution model and extend it with PIP
        if (PAR_model_indels) {

            bool computeFrequenciesFromData = false;

            // If frequencies are estimated from the data, but there is no alignment, then flag it.
            //if (!PAR_alignment) {
            std::string baseModel;

            std::map<std::string, std::string> basemodelMap;
            KeyvalTools::parseProcedure(modelMap["model"], baseModel, basemodelMap);

            std::vector<std::string> keys;
            for (auto it = basemodelMap.begin(); it != basemodelMap.end(); ++it) keys.push_back(it->first);

            if (!keys.empty()) {
                baseModel += "(";
                for (auto &key:keys) {
                    if (key != "initFreqs") {
                        baseModel += key + "=" + basemodelMap[key];
                    } else {
                        if (basemodelMap[key] == "observed") {
                            computeFrequenciesFromData = true;
                        }
                    }
                    baseModel += ",";
                }
                baseModel.pop_back();
                baseModel += ")";
                modelMap["model"] = baseModel;
            }

            //}

            // Instantiation of the canonical substitution model
            if (PAR_Alphabet.find("Codon") != std::string::npos || PAR_Alphabet.find("Protein") != std::string::npos) {
                smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabetNoGaps, gCode.get(), sites,
                                                                                  modelMap, "", true, false, 0);
            } else {
                smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites,
                                                                                  modelMap, "", true, false, 0);
            }

            // If PIP, then check if lambda/mu initial values are estimated from the data
            if (modelMap.find("estimated") != modelMap.end()) {
                ApplicationTools::displayError(
                        "The use of the tag [observed] is obsolete. Use the tag [initFreqs] instead");
                exit(1);
            } else if (modelMap.find("initFreqs") != modelMap.end()) {
                if (modelMap["initFreqs"] == "observed") {
                    estimatePIPparameters = true;
                }
            }

            if (estimatePIPparameters) {

//                if (PAR_alignment) {
//                    lambda = bpp::estimateLambdaFromData(tree, sequences, PAR_proportion);
//                    mu = bpp::estimateMuFromData(tree, PAR_proportion);
//
//                } else {
                lambda = bpp::estimateLambdaFromData(tree, sites);
                mu = bpp::estimateMuFromData(tree, sites);
                //}

                DLOG(INFO) << "[PIP model] Estimated PIP parameters from data using input sequences (lambda=" <<
                           lambda << ",mu=" << mu << "," "I=" << lambda * mu << ")";

            } else {
                lambda = (modelMap.find("lambda") == modelMap.end()) ? 0.1 : std::stod(modelMap["lambda"]);
                mu = (modelMap.find("mu") == modelMap.end()) ? 0.2 : std::stod(modelMap["mu"]);
            }

            DLOG(INFO) << "[PIP model] Fixed PIP parameters to (lambda=" << lambda << ",mu=" << mu << "," "I="
                       << lambda * mu << ")";

            // Instantiate the corrisponding PIP model given the alphabet
            if (PAR_alignment) {
                if (PAR_Alphabet.find("DNA") != std::string::npos && PAR_Alphabet.find("Codon") == std::string::npos) {

                    smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabet), smodel, *sequences, lambda, mu,
                                         computeFrequenciesFromData);

                } else if (PAR_Alphabet.find("Protein") != std::string::npos) {

                    smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabet), smodel, *sequences, lambda, mu,
                                        computeFrequenciesFromData);

                } else if (PAR_Alphabet.find("Codon") != std::string::npos) {

                    smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabet), gCode.get(), smodel,
                                           *sequences, lambda, mu,
                                           computeFrequenciesFromData);

                    ApplicationTools::displayWarning(
                            "Codon models are experimental in the current version... use with caution!");
                    DLOG(WARNING) << "CODONS activated byt the program is not fully tested under these settings!";
                }
            } else {
            if (PAR_Alphabet.find("DNA") != std::string::npos && PAR_Alphabet.find("Codon") == std::string::npos) {

                smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabet), smodel, *sites, lambda, mu,
                                     computeFrequenciesFromData);

            } else if (PAR_Alphabet.find("Protein") != std::string::npos) {

                smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabet), smodel, *sites, lambda, mu,
                                    computeFrequenciesFromData);

            } else if (PAR_Alphabet.find("Codon") != std::string::npos) {

                smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabet), gCode.get(), smodel,
                                       *sites, lambda, mu,
                                       computeFrequenciesFromData);

                ApplicationTools::displayWarning(
                        "Codon models are experimental in the current version... use with caution!");
                DLOG(WARNING) << "CODONS activated byt the program is not fully tested under these settings!";
            }
            }

        } else {
            // if the alphabet is not extended, then the gap character is not supported
            //if (!PAR_alignment) bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);
            bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);
            smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites,
                                                                              castorapp.getParams(), "", true, false, 0);
        }

        DLOG(INFO) << "[Substitution model] Number of states: " << (int) smodel->getNumberOfStates();

        ApplicationTools::displayResult("Substitution model", smodel->getName());
        if (PAR_model_indels)
            ApplicationTools::displayResult("Indel parameter initial value",
                                            (estimatePIPparameters) ? "estimated" : "fixed");

        ParameterList parameters = smodel->getParameters();
        for (size_t i = 0; i < parameters.size(); i++) {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }

        for (size_t i = 0; i < smodel->getFrequencies().size(); i++) {

            ApplicationTools::displayResult("eq.freq(" + smodel->getAlphabet()->getName(i) + ")",
                                            TextTools::toString(smodel->getFrequencies()[i], 4));
        }


        if (PAR_model_indels) {
            double pip_intensity =
                    parameters.getParameter("PIP.lambda").getValue() * parameters.getParameter("PIP.mu").getValue();
            ApplicationTools::displayResult("PIP.intensity", TextTools::toString(pip_intensity));
        }

        bpp::StdStr s1;
        bpp::PhylogeneticsApplicationTools::printParameters(smodel, s1, 1, true);
        DLOG(INFO) << s1.str();

        /////////////////////////
        // AMONG-SITE-RATE-VARIATION

        bpp::DiscreteDistribution *rDist = nullptr;

        // Among site rate variation (ASVR)
        if (smodel->getNumberOfStates() >= 2 * smodel->getAlphabet()->getSize()) {
            // Markov-modulated Markov model!
            rDist = new ConstantRateDistribution();
        } else {
            rDist = PhylogeneticsApplicationTools::getRateDistribution(castorapp.getParams());
        }

        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************
        // ********  3D DP PIP  ************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************


        /*
        FILE *fid = fopen("/Users/max/castor/data/WAG_Q","w");
        for(int i=0;i<21;i++){
            for(int j=0;j<21;j++){
                fprintf(fid,"%18.16lf ",smodel->getGenerator()(i,j));
            }
            fprintf(fid,"\n");
        }
        fclose(fid);

        fid = fopen("/Users/max/castor/data/WAG_Pi","w");
        for(int i=0;i<21;i++) {
            fprintf(fid,"%18.16lf\n",smodel->getFrequencies()[i]);
        }
        fclose(fid);

        exit(0);
        */

        /*
        FILE *fid=fopen("Q.data","w");
        for(int i=0;i<20;i++){
            for(int j=0;j<20;j++) {
                if(i==j){
                    double val=smodel->getGenerator().operator()(i,j)+smodel->getParameter("mu").getValue();

                    fprintf(fid,"%18.16lf ",val);
                }else{
                    fprintf(fid,"%18.16lf ",smodel->getGenerator().operator()(i,j));
                }
            }
            fprintf(fid,"\n");
        }
        fclose(fid);
    */

        // COMPUTE ALIGNMENT USING PROGRESSIVE-PIP
        pPIP *alignment = nullptr;
        progressivePIP *proPIP = nullptr;



        if (PAR_alignment) {

            //std::string PAR_output_file_lk = ApplicationTools::getAFilePath("output.lk.file", castorapp.getParams(), false,
            //                                                                false, "", true, "", 1);

            std::string PAR_output_file_msa = ApplicationTools::getAFilePath("output.msa.file", castorapp.getParams(), false,
                                                                             false, "", true, "", 1);

            std::string PAR_alignment_version = ApplicationTools::getStringParameter("alignment.version",
                                                                                     castorapp.getParams(), "cpu", "",
                                                                                     true, 0);

            int PAR_alignment_sbsolutions = ApplicationTools::getIntParameter("alignment.sb_solutions",
                                                                              castorapp.getParams(), 1, "", true, 0);

            double PAR_alignment_sbtemperature = ApplicationTools::getDoubleParameter("alignment.sb_temperature",
                                                                              castorapp.getParams(), 1.0, "", true, 0);


            ApplicationTools::displayMessage("\n[Computing the multi-sequence alignment]");
            //ApplicationTools::displayResult("Proportion gappy sites", TextTools::toString(PAR_proportion, 4));
            ApplicationTools::displayResult("Aligner optimised for:", PAR_alignment_version);

            ApplicationTools::displayBooleanResult("Stochastic backtracking active", PAR_alignment_sbsolutions > 1);
            if (PAR_alignment_sbsolutions > 1) {
                ApplicationTools::displayResult("Number of stochastic solutions:",
                                                TextTools::toString(PAR_alignment_sbsolutions));
            }

            DLOG(INFO) << "[Alignment sequences] Starting MSA_t inference using Pro-PIP...";

            // Execute alignment on post-order node list
            std::vector<tshlib::VirtualNode *> ftn = utree->getPostOrderNodeList();


            int num_sb;// = 0; // number of sub-optimal MSAs
            double temperature;// = 1; // temperature for SB version

            enumDP3Dversion DPversion = CPU; // DP3D version

            if (PAR_alignment_version.find("cpu") != std::string::npos) {
                DPversion = CPU; // slower but uses less memory
                num_sb = 1;
            } else if (PAR_alignment_version.find("ram") != std::string::npos) {
                DPversion = RAM; // faster but uses more memory
                num_sb = 1;
            } else if (PAR_alignment_version.find("sb") != std::string::npos) {
                DPversion = SB;  // stochastic backtracking version
                num_sb = PAR_alignment_sbsolutions;
                temperature = PAR_alignment_sbtemperature;
            } else {
                ApplicationTools::displayError("The user specified an unknown alignment.version. "
                                               "The execution will not continue.");
            }

            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            proPIP = new bpp::progressivePIP(utree,               // tshlib tree
                                             tree,                // bpp tree
                                             smodel,              // substitution model
                                             tm,                  // tree-map
                                             sequences,           // un-aligned input sequences
                                             rDist,               // rate-variation among site distribution
                                             castorapp.getSeed());  // seed for random number generation

            proPIP->_initializePIP(ftn,          // list of tshlib nodes in on post-order (correct order of execution)
                                   DPversion,    // version of the alignment algorithm
                                   num_sb,       // number of suboptimal MSAs
                                   temperature); // to tune the greedyness of the sub-optimal solution

            proPIP->PIPnodeAlign(); // align input sequences with a DP algorithm under PIP

            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
            std::cout << "\nAlignment elapsed time (msec): " << duration << std::endl;
            ApplicationTools::displayResult("\nAlignment elapsed time (msec):",
                                            TextTools::toString((double) duration, 4));


            // convert PIPmsa into a sites objects
            sites = PIPmsaUtils::PIPmsa2Sites(proPIP->alphabet_,
                                              *(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getseqNames()),
                                              *(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getMSA()));


            // Export alignment to file
            if (PAR_output_file_msa.find("none") == std::string::npos) {
                DLOG(INFO) << "[Alignment sequences]\t The final alignment can be found in " << PAR_output_file_msa;
                bpp::Fasta seqWriter;
                seqWriter.writeAlignment(TextUtils::appendToFilePath(PAR_output_file_msa, "initial"), *sites, true);
            }

            double MSAscore = proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getScore();

            ApplicationTools::displayResult("Alignment log likelihood", TextTools::toString(MSAscore, 15));


        }
        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************
        //**********************************************************************************************************

        /////////////////////////
        // Homogeneous modeling - initialization likelihood functions

        ApplicationTools::displayMessage("\n[Setting up likelihood functions]");

        // Initialization likelihood functions
        bpp::AbstractHomogeneousTreeLikelihood *tl;

        // Get transition model from substitution model
        if (!PAR_model_indels) {
            model = bpp::PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), sites,
                                                                           castorapp.getParams(), "", true, false, 0);
        } else {
            unique_ptr<TransitionModel> test;
            test.reset(smodel);
            model = test.release();
        }

        // Initialise likelihood functions
        if (!PAR_model_indels) {
            //tl = new bpp::RHomogeneousTreeLikelihood_Generic(*tree, *sites, model, rDist, false, false, false);
            tl = new bpp::UnifiedTSHomogeneousTreeLikelihood(*tree, *sites, model, rDist, utree, &tm, true,
                                                             castorapp.getParams(), "", false, false,
                                                             false);

        } else {
            //tl = new bpp::RHomogeneousTreeLikelihood_PIP(*tree, *sites, model, rDist, &tm, false, false, false);
            tl = new bpp::UnifiedTSHomogeneousTreeLikelihood_PIP(*tree, *sites, model, rDist, utree, &tm, true,
                                                                 castorapp.getParams(), "", false,
                                                                 false, false);
        }
        ApplicationTools::displayResult("Tree likelihood model", std::string("Homogeneous"));


        tl->initialize();

        /////////////////////////
        // Parameter sanity check

        ApplicationTools::displayMessage("\n[Parameter sanity check]");

        //Listing parameters
        string paramNameFile = ApplicationTools::getAFilePath("output.parameter_names.file", castorapp.getParams(), false,
                                                              false, "", true, "none", 1);
        if (paramNameFile != "none") {
            ApplicationTools::displayResult("List parameters to", paramNameFile);
            ofstream pnfile(paramNameFile.c_str(), ios::out);
            ParameterList pl = tl->getParameters();
            for (size_t i = 0; i < pl.size(); ++i) {
                pnfile << pl[i].getName() << endl;
            }
            pnfile.close();
        }

        //Check initial likelihood:
        double logL = tl->getValue();
        if (std::isinf(logL)) {
            // This may be due to null branch lengths, leading to null likelihood!
            ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
            ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
            ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
            ParameterList pl = tl->getBranchLengthsParameters();
            for (unsigned int i = 0; i < pl.size(); i++) {
                if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
            }
            tl->matchParametersValues(pl);
            logL = tl->getLogLikelihood();
        }
        ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
        if (std::isinf(logL)) {
            ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
            if (codonAlphabet) {
                bool f = false;
                size_t s;
                for (size_t i = 0; i < sites->getNumberOfSites(); i++) {
                    if (std::isinf(tl->getLogLikelihoodForASite(i))) {
                        const Site &site = sites->getSite(i);
                        s = site.size();
                        for (size_t j = 0; j < s; j++) {
                            if (gCode->isStop(site.getValue(j))) {
                                (*ApplicationTools::error << "Stop Codon at site " << site.getPosition()
                                                          << " in sequence "
                                                          << sites->getSequence(j).getName()).endLine();
                                f = true;
                            }
                        }
                    }
                }
                if (f)
                    exit(-1);
            }
            bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites",
                                                                         castorapp.getParams(), false, "",
                                                                         true, 1);
            if (!removeSaturated) {
                ofstream debug("DEBUG_likelihoods.txt", ios::out);
                for (size_t i = 0; i < sites->getNumberOfSites(); i++) {
                    debug << "Position " << sites->getSite(i).getPosition() << " = " << tl->getLogLikelihoodForASite(i)
                          << endl;
                }
                debug.close();
                ApplicationTools::displayError(
                        "!!! Site-specific likelihood have been written in file DEBUG_likelihoods.txt .");
                ApplicationTools::displayError(
                        "!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
                ApplicationTools::displayError(
                        "!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
                exit(1);
            } else {
                ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
                for (size_t i = sites->getNumberOfSites(); i > 0; --i) {
                    if (std::isinf(tl->getLogLikelihoodForASite(i - 1))) {
                        ApplicationTools::displayResult("Ignore saturated site", sites->getSite(i - 1).getPosition());
                        sites->deleteSite(i - 1);
                    }
                }
                ApplicationTools::displayResult("Number of sites retained", sites->getNumberOfSites());
                tl->setData(*sites);
                tl->initialize();
                logL = tl->getValue();
                if (std::isinf(logL)) {
                    throw Exception("Likelihood is still 0 after saturated sites are removed! Looks like a bug...");
                }
                ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
            }
        }


        /////////////////////////
        // OPTIMISE PARAMETERS (numerical + topology) according to user parameters
        // Optimise parameters automatically following standard pipeline

        ApplicationTools::displayMessage("\n[Executing numerical parameters and topology optimization]");

        tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(Optimizators::optimizeParameters(tl,
                                                                                                tl->getParameters(),
                                                                                                castorapp.getParams(),
                                                                                                "",
                                                                                                true,
                                                                                                true,
                                                                                                0));

        logL = tl->getLogLikelihood();


        /////////////////////////
        // OUTPUT

        // Export final alignment
//        if (PAR_output_file_msa.find("none") == std::string::npos) {
//            ApplicationTools::displayResult("\n\nOutput alignment to file", PAR_output_file_msa);
//            DLOG(INFO) << "[Output alignment]\t The final alignment can be found in " << PAR_output_file_msa;
//            bpp::Fasta seqWriter;
//            seqWriter.writeAlignment(PAR_output_file_msa, *sites, true);
//        }

        delete sequences;

        // Export final tree (if nexus is required, then our re-implementation of the the nexus writer is called)
        tree = new TreeTemplate<Node>(tl->getTree());

        std::string PAR_output_tree_format = ApplicationTools::getStringParameter("output.tree.format",
                                                                                  castorapp.getParams(),
                                                                                  "Newick",
                                                                                  "",
                                                                                  true,
                                                                                  true);
        if (PAR_output_tree_format.find("Nexus") != std::string::npos) {
            std::vector<Tree *> tmp;
            tmp.push_back(tree);
            OutputUtils::writeNexusMetaTree(tmp, castorapp.getParams());
        } else {
            PhylogeneticsApplicationTools::writeTree(*tree, castorapp.getParams());
        }


        // Export annotation file (tab separated values)
        std::string PAR_output_annotation_file = ApplicationTools::getAFilePath("output.annotation.file",
                                                                                castorapp.getParams(),
                                                                                false,
                                                                                false,
                                                                                "",
                                                                                true,
                                                                                "",
                                                                                1);
        if (PAR_output_annotation_file.find("none") == std::string::npos) {
            ApplicationTools::displayResult("Output annotation to file", PAR_output_annotation_file);
            OutputUtils::writeTreeAnnotations2TSV(tree, PAR_output_annotation_file);

        }

        // Write parameters to screen:
        ApplicationTools::displayResult("Final Log likelihood", TextTools::toString(logL, 15));

        parameters = tl->getSubstitutionModelParameters();
        for (size_t i = 0; i < parameters.size(); i++) {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }
        if (PAR_model_indels) {
            double pip_intensity =
                    parameters.getParameter("PIP.lambda").getValue() * parameters.getParameter("PIP.mu").getValue();
            ApplicationTools::displayResult("PIP.intensity", TextTools::toString(pip_intensity));
        }
        parameters = tl->getRateDistributionParameters();
        for (size_t i = 0; i < parameters.size(); i++) {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }

        // Checking convergence:
        PhylogeneticsApplicationTools::checkEstimatedParameters(tl->getParameters());

        // Write parameters to file (according to arguments)
        OutputUtils::exportOutput(tl, sites, castorapp.getParams());


        // Compute support measures
        std::string PAR_support = ApplicationTools::getStringParameter("support", castorapp.getParams(), "", "", true,
                                                                       true);
        if (PAR_support == "bootstrap") {
            ApplicationTools::displayMessage("\n[Tree support measures]");

            bpp::Bootstrap(tl, *sites, rDist, utree, &tm, castorapp.getParams(), "support.");
        }

        // Delete objects and free memory
        delete alphabet;
        delete alphabetNoGaps;
        delete sites;
        delete rDist;
        delete tl;
        delete tree;

        castorapp.done();
        google::ShutdownGoogleLogging();
        exit(0);

    } catch (exception &e) {
        cout << e.what() << endl;
        google::ShutdownGoogleLogging();
        exit(1);
    }
}
