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

#include <glog/logging.h>
#include "Version.hpp"
#include "progressivePIP.hpp"
#include "CastorApplication.hpp"

int main(int argc, char *argv[]) {

    /* ***************************************************
    * Standard workflow
    * [INPUT]
    * 1. tree + alignment => (1.1) just parse everything
    * 2. alignment  => (2.1) parse alignment => (2.2) generate tree using bioNJ
    * 3. sequences  => (3.1) parse sequences => (3.2) generate tree using bioNJ => (3.3) perform alignment
    * (it should not be supported in production)
    * 4. sequences + tree => (4.1) parse sequence => (4.2) parse tree => (4.3) perform alignment
    */

    //FLAGS_logtostderr = 1;
    FLAGS_log_dir = ".";
    google::InitGoogleLogging(software::name.c_str());
    google::InstallFailureSignalHandler();

    std::map<std::string, std::string> modelMap;
    std::unique_ptr<GeneticCode> gCode;
    UtreeBppUtils::treemap tm;
    ParameterList parameters;

    bpp::Alphabet *alphabet = nullptr;
    bpp::Alphabet *alphabetNoGaps = nullptr;
    bpp::SequenceContainer *sequences = nullptr;
    bpp::SiteContainer *sites = nullptr;
    bpp::Tree *tree = nullptr;
    UtreeBppUtils::Utree *utree = nullptr;
    //bpp::SubstitutionModel *smodel = nullptr;
    bpp::TransitionModel *dmodel = nullptr;
    bpp::TransitionModel *model = nullptr;
    bpp::DiscreteDistribution *rDist = nullptr;
    bpp::AbstractHomogeneousTreeLikelihood *tl = nullptr;

    unique_ptr<TransitionModel> test;

    try {

        /////////////////////////////////////////
        // START THE APP

        bpp::CastorApplication castorapp(argc,
                                         argv,
                                         std::string(software::name + " " + software::version),
                                         std::string(software::releasegitbranch + " " + software::releasegitref),
                                         std::string(software::releasedate + ", " + software::releasetime));



        castorapp.checkInputArgument(argc);

        bpp::ApplicationTools::displayResult("Random seed set to", castorapp.getSeed());
        ApplicationTools::displayResult("Log files location", std::string("current execution path"));

        /////////////////////////////////////////
        // CLI ARGUMENTS

        modelMap = castorapp.getCLIarguments();

        /////////////////////////////////////////
        // ALPHABET

        // The alphabetNoGaps object contains the not-extended alphabet as requested by the user,
        // while alphabet contains the extended version of the same alphabet.
        alphabetNoGaps = castorapp.getAlphabetNoGaps();

        alphabet = castorapp.getAlphabet();

        // Alphabet used for codon models
        gCode = castorapp.getGcode(alphabetNoGaps);

        ApplicationTools::displayResult("Alphabet", TextTools::toString(alphabetNoGaps->getAlphabetType()));

        ApplicationTools::displayBooleanResult("Allow gaps as extra character", castorapp.PAR_model_indels_);

        DLOG(INFO) << "alphabet:  " << castorapp.PAR_alphabet_ << " | gap-extention " << (int) castorapp.PAR_model_indels_;

        /////////////////////////////////////////
        // GET DATA

        ApplicationTools::displayMessage("\n[Preparing input data]");

        ApplicationTools::displayBooleanResult("Aligned sequences", !castorapp.PAR_alignment_);

        sequences = castorapp.getSequences(alphabet);

        sites = castorapp.getSites(alphabet);

        ApplicationTools::displayResult("Number of sequences",TextTools::toString(sites->getNumberOfSequences()));











        /////////////////////////////////////////
        /////////////////////////////////////////
        /////////////////////////////////////////
        /////////////////////////////////////////
        /////////////////////////////////////////
        /////////////////////////////////////////


        bool flag_infere_tree = false;
        bool flag_infere_indel_rates = false;

        if (castorapp.PAR_init_tree_opt_ == "user") {
            flag_infere_tree = false;
        } else if (castorapp.PAR_init_tree_opt_ == "random") {
            flag_infere_tree = false;
        } else if (castorapp.PAR_init_tree_opt_ == "distance") {
            flag_infere_tree = true;
        } else{
            throw Exception("Unknown init tree method.");
        }

        if (castorapp.PAR_model_indels_) {
            if ((modelMap.find("lambda") != modelMap.end()) && (modelMap.find("mu") != modelMap.end())) {
                flag_infere_indel_rates = false;
            } else {
                flag_infere_indel_rates = true;
            }
        }else{

            flag_infere_indel_rates = false;

        }

        if (castorapp.PAR_model_indels_) {

            //********************************************************//
            // indel model
            //********************************************************//

            if (!flag_infere_indel_rates && !flag_infere_tree) {
                //========================================================//
                // inde rates : given | tree : given
                //========================================================//

                // Get Indel Rates
                castorapp.lambda_ = (modelMap.find("lambda") == modelMap.end()) ? 1.0 : std::stod(modelMap["lambda"]);
                castorapp.mu_ = (modelMap.find("mu") == modelMap.end()) ? 0.1 : std::stod(modelMap["mu"]);

                // Get Model
                // Instantiate a substitution model and extend it with PIP
                castorapp.smodel = castorapp.getSubstitutionModelIndel(modelMap,gCode,alphabet,alphabetNoGaps,sites,sequences);
                DLOG(INFO) << "[Substitution model] Number of states: " << (int) castorapp.smodel->getNumberOfStates();
                ApplicationTools::displayResult("Substitution model", castorapp.smodel->getName());


                //dmodel = castorapp.getTransitionModel(smodel,alphabetDistMethod,gCode,sitesDistMethod,parmap);
                test.reset(castorapp.smodel);
                dmodel = test.release();

                // AMONG-SITE-RATE-VARIATION
                rDist = castorapp.getASVR(castorapp.smodel);

                // Get Tree
                ApplicationTools::displayMessage("\n[Preparing initial tree]");
                castorapp.PAR_init_tree_opt_ = ApplicationTools::getStringParameter("init.tree", castorapp.getParams(), "user", "", false,1);
                ApplicationTools::displayResult("Initial tree", castorapp.PAR_init_tree_opt_);
                if (castorapp.PAR_init_tree_opt_ == "user") {
                    tree = castorapp.getUserTree();
                } else if (castorapp.PAR_init_tree_opt_ == "random") {
                    tree = castorapp.getRandomTree(sites);
                } else{
                    throw Exception("Unknown init tree method.");
                }
                // If the tree has multifurcation, then resolve it with midpoint rooting
                tree = castorapp.resolveMultifurcation(tree);
                // Rename internal nodes with standard Vxx * where xx is a progressive number
                castorapp.renameInternalNodes(tree);
                // Try to write the current tree to file. This will be overwritten by the optimized tree,
                // but allow to check file existence before running optimization!
                PhylogeneticsApplicationTools::writeTree(*tree, castorapp.getParams());
                // Setting branch lengths?
                castorapp.initBranchLength(tree);
                utree = castorapp.getUtree(tree,sites,sequences,tm);
                DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(tree);


            } else if (!flag_infere_indel_rates && flag_infere_tree) {
                //========================================================//
                // inde rates : given | tree : infere
                //========================================================//

                // Get Indel Rates
                castorapp.lambda_ = (modelMap.find("lambda") == modelMap.end()) ? 1.0 : std::stod(modelMap["lambda"]);
                castorapp.mu_ = (modelMap.find("mu") == modelMap.end()) ? 0.1 : std::stod(modelMap["mu"]);


                // Get Model
                // Instantiate a substitution model and extend it with PIP
                castorapp.smodel = castorapp.getSubstitutionModelIndel(modelMap,gCode,alphabet,alphabetNoGaps,sites,sequences);
                DLOG(INFO) << "[Substitution model] Number of states: " << (int) castorapp.smodel->getNumberOfStates();
                ApplicationTools::displayResult("Substitution model", castorapp.smodel->getName());

                //dmodel = castorapp.getTransitionModel(smodel,alphabetDistMethod,gCode,sitesDistMethod,parmap);
                test.reset(castorapp.smodel);
                dmodel = test.release();

                // AMONG-SITE-RATE-VARIATION
                rDist = castorapp.getASVR(castorapp.smodel);

                // Get Tree
                ApplicationTools::displayMessage("\n[Preparing initial tree]");
                castorapp.PAR_init_tree_opt_ = ApplicationTools::getStringParameter("init.tree", castorapp.getParams(), "user", "", false,1);
                ApplicationTools::displayResult("Initial tree", castorapp.PAR_init_tree_opt_);
                if (castorapp.PAR_init_tree_opt_ == "distance") {

                    bpp::TransitionModel *dmodel_copy = dmodel->clone(); /// ??? why ????
                    bpp::DiscreteDistribution *rDist_copy  = rDist->clone(); /// ??? why ????

                    tree = castorapp.getDistanceTree(dmodel_copy,sites,sequences,alphabet,alphabetNoGaps,modelMap,gCode,rDist_copy);

                    //delete dmodel_copy;
                    //delete rDist_copy;

                } else{
                    throw Exception("Unknown init tree method.");
                }
                // If the tree has multifurcation, then resolve it with midpoint rooting
                tree = castorapp.resolveMultifurcation(tree);
                // Rename internal nodes with standard Vxx * where xx is a progressive number
                castorapp.renameInternalNodes(tree);
                // Try to write the current tree to file. This will be overwritten by the optimized tree,
                // but allow to check file existence before running optimization!
                PhylogeneticsApplicationTools::writeTree(*tree, castorapp.getParams());
                // Setting branch lengths?
                castorapp.initBranchLength(tree);
                utree = castorapp.getUtree(tree,sites,sequences,tm);
                DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(tree);



            } else if (flag_infere_indel_rates && !flag_infere_tree) {
                //========================================================//
                // inde rates : infere | tree : given
                //========================================================//

                // Get Tree
                ApplicationTools::displayMessage("\n[Preparing initial tree]");
                castorapp.PAR_init_tree_opt_ = ApplicationTools::getStringParameter("init.tree", castorapp.getParams(), "user", "", false,1);
                ApplicationTools::displayResult("Initial tree", castorapp.PAR_init_tree_opt_);
                if (castorapp.PAR_init_tree_opt_ == "user") {
                    tree = castorapp.getUserTree();
                } else if (castorapp.PAR_init_tree_opt_ == "random") {
                    tree = castorapp.getRandomTree(sites);
                } else{
                    throw Exception("Unknown init tree method.");
                }
                // If the tree has multifurcation, then resolve it with midpoint rooting
                tree = castorapp.resolveMultifurcation(tree);
                // Rename internal nodes with standard Vxx * where xx is a progressive number
                castorapp.renameInternalNodes(tree);
                // Try to write the current tree to file. This will be overwritten by the optimized tree,
                // but allow to check file existence before running optimization!
                PhylogeneticsApplicationTools::writeTree(*tree, castorapp.getParams());
                // Setting branch lengths?
                castorapp.initBranchLength(tree);
                utree = castorapp.getUtree(tree,sites,sequences,tm);
                DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(tree);


                // Get Indel Rates
                inference_indel_rates::infere_indel_rates_from_sequences(castorapp.PAR_input_sequences_,castorapp.PAR_alphabet_,tree,castorapp.lambda_,castorapp.mu_,gCode.get(),modelMap);
                DLOG(INFO) << "[PIP model] Estimated PIP parameters from data using input sequences (lambda=" <<castorapp.lambda_ << ",mu=" << castorapp.mu_ << "," "I=" << castorapp.lambda_ * castorapp.mu_ << ")";


                // Get Model
                castorapp.smodel = castorapp.getSubstitutionModelIndel(modelMap,gCode,alphabet,alphabetNoGaps,sites,sequences);

                //dmodel = castorapp.getTransitionModel(smodel,alphabetDistMethod,gCode,sitesDistMethod,parmap);
                test.reset(castorapp.smodel);
                dmodel = test.release();

                // AMONG-SITE-RATE-VARIATION
                rDist = castorapp.getASVR(castorapp.smodel);

            } else if (flag_infere_indel_rates && flag_infere_tree) {
                //========================================================//
                // inde rates : infere | tree : infere
                //========================================================//

                // Get Tree
                ApplicationTools::displayMessage("\n[Preparing initial tree]");
                castorapp.PAR_init_tree_opt_ = ApplicationTools::getStringParameter("init.tree", castorapp.getParams(), "user", "", false,1);
                ApplicationTools::displayResult("Initial tree", castorapp.PAR_init_tree_opt_);
                if (castorapp.PAR_init_tree_opt_ == "distance") {


                    //????
                    //tree = castorapp.getDistanceTree(sites,sequences,alphabet,alphabetNoGaps,modelMap,gCode);




                } else{
                    throw Exception("Unknown init tree method.");
                }
                // If the tree has multifurcation, then resolve it with midpoint rooting
                tree = castorapp.resolveMultifurcation(tree);
                // Rename internal nodes with standard Vxx * where xx is a progressive number
                castorapp.renameInternalNodes(tree);
                // Try to write the current tree to file. This will be overwritten by the optimized tree,
                // but allow to check file existence before running optimization!
                PhylogeneticsApplicationTools::writeTree(*tree, castorapp.getParams());
                // Setting branch lengths?
                castorapp.initBranchLength(tree);
                utree = castorapp.getUtree(tree,sites,sequences,tm);
                DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(tree);

                // Get Model
                // Instantiate a substitution model and extend it with PIP
                castorapp.smodel = castorapp.getSubstitutionModelIndel(modelMap,gCode,alphabet,alphabetNoGaps,sites,sequences);
                DLOG(INFO) << "[Substitution model] Number of states: " << (int) castorapp.smodel->getNumberOfStates();
                ApplicationTools::displayResult("Substitution model", castorapp.smodel->getName());

                //dmodel = castorapp.getTransitionModel(smodel,alphabetDistMethod,gCode,sitesDistMethod,parmap);
                test.reset(castorapp.smodel);
                dmodel = test.release();

                // AMONG-SITE-RATE-VARIATION
                rDist = castorapp.getASVR(castorapp.smodel);

            } else {
                DLOG(WARNING) << "Problem encountered, exit!";
            }


        }else{

            //********************************************************//
            // no indel model
            //********************************************************//

            // Get Model
            castorapp.smodel = castorapp.getSubstitutionModelNoIndel(gCode,sites,alphabet);



            //????????
            //dmodel = castorapp.getTransitionModel(smodel,alphabetDistMethod,gCode,sitesDistMethod,parmap);





            // AMONG-SITE-RATE-VARIATION
            rDist = castorapp.getASVR(castorapp.smodel);

            // Get Tree
            ApplicationTools::displayMessage("\n[Preparing initial tree]");
            castorapp.PAR_init_tree_opt_ = ApplicationTools::getStringParameter("init.tree", castorapp.getParams(), "user", "", false,1);
            ApplicationTools::displayResult("Initial tree", castorapp.PAR_init_tree_opt_);
            if (castorapp.PAR_init_tree_opt_ == "user") {
                tree = castorapp.getUserTree();
            } else if (castorapp.PAR_init_tree_opt_ == "random") {
                tree = castorapp.getRandomTree(sites);
            } else if (castorapp.PAR_init_tree_opt_ == "distance") {
                tree = castorapp.getDistanceTree(dmodel,sites,sequences,alphabet,alphabetNoGaps,modelMap,gCode,rDist);
            } else{
                throw Exception("Unknown init tree method.");
            }
            // If the tree has multifurcation, then resolve it with midpoint rooting
            tree = castorapp.resolveMultifurcation(tree);
            // Rename internal nodes with standard Vxx * where xx is a progressive number
            castorapp.renameInternalNodes(tree);
            // Try to write the current tree to file. This will be overwritten by the optimized tree,
            // but allow to check file existence before running optimization!
            PhylogeneticsApplicationTools::writeTree(*tree, castorapp.getParams());
            // Setting branch lengths?
            castorapp.initBranchLength(tree);
            utree = castorapp.getUtree(tree,sites,sequences,tm);
            DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(tree);

        }


        /////////////////////////////////////////
        /////////////////////////////////////////
        /////////////////////////////////////////
        /////////////////////////////////////////
        /////////////////////////////////////////
        /////////////////////////////////////////






        /*
        /////////////////////////////////////////
        // INITIAL TREE
        ApplicationTools::displayMessage("\n[Preparing initial tree]");

        tree = castorapp.getTree(smodel,sites,sequences,alphabet,alphabetNoGaps,modelMap,gCode);

        utree = castorapp.getUtree(tree,sites,sequences,tm);

        DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(tree);

        ApplicationTools::displayResult("Initial tree total length", TextTools::toString(tree->getTotalLength(), 6));

        /////////////////////////////////////////
        // SUBSTITUTION MODEL

        ApplicationTools::displayMessage("\n[Setting up substitution model]");


        smodel = castorapp.getSubstitutionModel(modelMap,gCode,alphabet,alphabetNoGaps,sites,sequences,tree);
        */








        /////////////////////////////////////////
        // GET PARAMETERS

        parameters = castorapp.getParameters(castorapp.smodel);

        CastorApplicationUtils::printParameters(parameters,castorapp.smodel);

        /*
        /////////////////////////////////////////
        // AMONG-SITE-RATE-VARIATION

        rDist = castorapp.getASVR(smodel);
        */

        /////////////////////////////////////////
        // COMPUTE ALIGNMENT USING PROGRESSIVE-PIP

        if (castorapp.PAR_alignment_) {

            sites = castorapp.getMSA(sequences,rDist,castorapp.smodel,tm,tree,utree);

        }

        /////////////////////////////////////////
        // HOMOGENEOUS MODELING - INITIALIZATION LIKELIHOOD FUNCTIONS

        ApplicationTools::displayMessage("\n[Setting up likelihood functions]");

        tl = castorapp.initLK(model,castorapp.smodel,tm,rDist,gCode,alphabet,sites,tree,utree);

        /////////////////////////////////////////
        // PARAMETER SANITY CHECK

        ApplicationTools::displayMessage("\n[Parameter sanity check]");

        castorapp.getParSanityCheck(tl,sites,gCode);

        /////////////////////////////////////////
        // OPTIMISE PARAMETERS

        ApplicationTools::displayMessage("\n[Executing numerical parameters and topology optimization]");

        // Optimise parameters (numerical + topology) automatically following standard pipeline according to user parameters
        castorapp.getOptParams(tl);

        /////////////////////////////////////////
        // BOOTSTRAP

        castorapp.getBootstrap(utree,tl,sites,tm,rDist);

        /////////////////////////////////////////
        // OUTPUT

        castorapp.output(tree,utree,sites,tl,tm,parameters,rDist);

        /////////////////////////////////////////
        // DELETE OBJECTS AND FREE MEMORY

        delete castorapp.smodel;

        delete sequences;
        delete alphabet;
        delete alphabetNoGaps;
        delete sites;
        delete rDist;
        delete tl;
        delete tree;
        delete utree;
        delete model;

        /////////////////////////////////////////
        // CLOSE APP
        castorapp.done();
        google::ShutdownGoogleLogging();
        exit(EXIT_SUCCESS);

    } catch (exception &e) {
        cout << e.what() << endl;
        google::ShutdownGoogleLogging();
        exit(EXIT_FAILURE);
    }
}
