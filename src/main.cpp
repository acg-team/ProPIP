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
    bpp::SubstitutionModel *smodel = nullptr;
    bpp::TransitionModel *model = nullptr;
    bpp::DiscreteDistribution *rDist = nullptr;
    bpp::AbstractHomogeneousTreeLikelihood *tl = nullptr;

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
        // INITIAL TREE
        ApplicationTools::displayMessage("\n[Preparing initial tree]");

        tree = castorapp.getTree(sites,sequences,alphabet,alphabetNoGaps,modelMap,gCode);

        utree = castorapp.getUtree(tree,sites,sequences,tm);

        DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(tree);

        ApplicationTools::displayResult("Initial tree total length", TextTools::toString(tree->getTotalLength(), 6));

        /////////////////////////////////////////
        // SUBSTITUTION MODEL

        ApplicationTools::displayMessage("\n[Setting up substitution model]");

        smodel = castorapp.getSubstitutionModel(modelMap,gCode,alphabet,alphabetNoGaps,sites,sequences,tree);

        /////////////////////////////////////////
        // GET PARAMETERS

        parameters = castorapp.getParameters(smodel);

        CastorApplicationUtils::printParameters(parameters,smodel);

        /////////////////////////////////////////
        // AMONG-SITE-RATE-VARIATION

        rDist = castorapp.getASVR(smodel);

        /////////////////////////////////////////
        // COMPUTE ALIGNMENT USING PROGRESSIVE-PIP

        if (castorapp.PAR_alignment_) {

            sites = castorapp.getMSA(sequences,rDist,smodel,tm,tree,utree);

        }

        /////////////////////////////////////////
        // HOMOGENEOUS MODELING - INITIALIZATION LIKELIHOOD FUNCTIONS

        ApplicationTools::displayMessage("\n[Setting up likelihood functions]");

        tl = castorapp.initLK(model,smodel,tm,rDist,gCode,alphabet,sites,tree,utree);

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

        delete sequences;
        delete alphabet;
        delete alphabetNoGaps;
        delete sites;
        delete rDist;
        delete tl;
        delete tree;
        delete utree;
        delete smodel;
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
