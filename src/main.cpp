/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2019 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of ProPIP
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

/*
* From miniJati:
*/
#include "CastorApplication.hpp"
#include "Version.hpp"
#include "PIP.hpp"

using namespace tshlib;

int main(int argc, char *argv[]) {

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

        bpp::ApplicationTools::displayResult("Alphabet", bpp::TextTools::toString(castorapp.alphabetNoGaps->getAlphabetType()));
        bpp::ApplicationTools::displayBooleanResult("Allow gaps as extra character", castorapp.PAR_model_indels);
        bpp::ApplicationTools::displayResult("Genetic Code", castorapp.codeDesc);

        DLOG(INFO) << "alphabet:  " << castorapp.PAR_Alphabet << " | gap-extention " << (int) castorapp.PAR_model_indels;

        //////////////////////////////////////////////
        // DATA
        bpp::ApplicationTools::displayMessage("\n[Preparing input data]");

        castorapp.getData();

        bpp::ApplicationTools::displayBooleanResult("Aligned sequences", !castorapp.PAR_alignment);
        if (castorapp.PAR_alignment) {
            bpp::ApplicationTools::displayResult("Number of sequences",bpp::TextTools::toString(castorapp.sequences->getNumberOfSequences()));
        } else {
            bpp::ApplicationTools::displayResult("Number of sequences",bpp::TextTools::toString(castorapp.sites->getNumberOfSequences()));
            bpp::ApplicationTools::displayResult("Number of sites", bpp::TextTools::toString(castorapp.sites->getNumberOfSites()));
        }

        /////////////////////////////////////////
        // INITIAL TREE
        bpp::ApplicationTools::displayMessage("\n[Preparing initial tree]");

        castorapp.getTree();

        bpp::ApplicationTools::displayResult("Initial tree total length", bpp::TextTools::toString(castorapp.tree->getTotalLength(), 6));

        /////////////////////////
        // SUBSTITUTION MODEL

        bpp::ApplicationTools::displayMessage("\n[Setting up substitution model]");
        
        castorapp.getSubstitutionModel();

        /////////////////////////
        // AMONG-SITE-RATE-VARIATION

        bpp::ApplicationTools::displayMessage("\n[Getting ASRV distribution]");

        castorapp.getASRV();

        /////////////////////////
        // ALIGN SEQUENCES

        if (castorapp.PAR_alignment) {

            bpp::ApplicationTools::displayMessage("\n[Computing the alignment]");

            castorapp.computeMSA();

            bpp::ApplicationTools::displayResult("Aligner optimised for:", castorapp.PAR_alignment_version);

            bpp::ApplicationTools::displayBooleanResult("Stochastic backtracking active", castorapp.PAR_alignment_sbsolutions > 1);

            if (castorapp.PAR_alignment_sbsolutions > 1) {
                bpp::ApplicationTools::displayResult("Number of stochastic solutions:",bpp::TextTools::toString(castorapp.PAR_alignment_sbsolutions));
            }

            bpp::ApplicationTools::displayResult("Alignment log likelihood", bpp::TextTools::toString(castorapp.proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getScore(), 15));
        }

        /////////////////////////
        // HOMOGENEOUS MODELING - initialization likelihood functions

        castorapp.initLkFun();

        /////////////////////////
        // PARAMETER SANITY CHECK

        bpp::ApplicationTools::displayMessage("\n[Parameter sanity check]");

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

        castorapp.bootstrapping();

        /////////////////////////
        // CLOSE APPLICATION
        castorapp.done();

        google::ShutdownGoogleLogging();

        exit(EXIT_SUCCESS);

    } catch (exception &e) {

        std::cout << e.what() << std::endl;

        google::ShutdownGoogleLogging();

        exit(EXIT_SUCCESS);
    }
}
