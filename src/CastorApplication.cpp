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
 * @file CastorApplication.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 06 02 2018
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
#include "CastorApplication.hpp"

#include <ctime>
#include <iostream>
#include <Bpp/Utils/AttributesTools.h>
#include <glog/logging.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/AlphabetIndex/DefaultNucleotideScore.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>

#include <Bpp/Seq/Io/Fasta.h>
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

#include "ExtendedAlphabet.hpp"
#include "Utils.hpp"
#include "UnifiedTSHTopologySearch.hpp"
#include "PIP.hpp"
#include "ExtendedAlphabet.hpp"
#include "RHomogeneousTreeLikelihood_PIP.hpp"
#include "RHomogeneousTreeLikelihood_Generic.hpp"
#include "Optimizators.hpp"
#include "SupportMeasures.hpp"
#include "UnifiedDistanceEstimation.hpp"

using namespace bpp;

CastorApplication::CastorApplication(int argc, char *argv[], const std::string &name, const std::string &strVersion, const std::string &build_date) :
        appName_(name), appBuild_(build_date), appVersion_(strVersion), params_(), timerStarted_(false) {

    params_ = bpp::AttributesTools::parseOptions(argc, argv);
    bool showversion = bpp::ApplicationTools::getBooleanParameter("version", params_, false, "", true, 3);
    bpp::ApplicationTools::warningLevel = bpp::ApplicationTools::getIntParameter("warning", params_, 0, "", true, 3);
    bool noint = bpp::ApplicationTools::getBooleanParameter("noninteractive", params_, false, "", true, 3);
    bpp::ApplicationTools::interactive = !noint;

    seed_ = bpp::ApplicationTools::getParameter<long>("seed", params_, -1, "", true, 3);
    if (seed_ >= 0) {
        bpp::RandomTools::setSeed(seed_);

    }else{
        unsigned int time_ui = (unsigned int) std::time(NULL);
        seed_ = time_ui;
        bpp::RandomTools::setSeed(seed_);
    }

    // Print version header

    if (showversion) {
        this->version();
        exit(0);
    }
}

void CastorApplication::startTimer() {
    ApplicationTools::startTimer();
    timerStarted_ = true;
}

void CastorApplication::done() {
   DLOG(INFO) << appName_ << "'s done. Bye.";
    if (timerStarted_)
        bpp::ApplicationTools::displayTime("Total execution time:");
}

void CastorApplication::getCLIarguments(){

    this->PAR_model_substitution = ApplicationTools::getStringParameter("model", this->getParams(), "JC69","", true, true);

    this->PAR_alignment = ApplicationTools::getBooleanParameter("alignment", this->getParams(), false);

    // Split model string description and test if PIP is required
    KeyvalTools::parseProcedure(this->PAR_model_substitution, this->modelStringName, this->modelMap);
    this->PAR_model_indels = (modelStringName == "PIP");

}

void CastorApplication::start(int argc){

    if (argc < 2) {
        this->help();
        exit(EXIT_SUCCESS);
    } else {
        this->banner();
        this->startTimer();
    }

}

void CastorApplication::getAlphabetIndel(){

    if (this->PAR_Alphabet.find("DNA") != std::string::npos) {
        this->alphabet = new bpp::DNA_EXTENDED();
        this->codonAlphabet = false;
    } else if(this->PAR_Alphabet.find("Codon") != std::string::npos){
        this->alphabet = new bpp::DNA_EXTENDED();
        // This is additional to the alphabet instance
        this->alphabet = new CodonAlphabet_Extended(dynamic_cast<bpp::NucleicAlphabet *>(this->alphabet));
        this->codonAlphabet = true;

        this->codeDesc = ApplicationTools::getStringParameter("genetic_code", this->getParams(), "Standard","", true, true);
        this->gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(dynamic_cast<bpp::CodonAlphabet *>(this->alphabetNoGaps)->getNucleicAlphabet(),this->codeDesc));

    } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
        this->alphabet = new bpp::ProteicAlphabet_Extended();
        this->codonAlphabet = false;
    }

}

void CastorApplication::getAlphabetNoIndel(){

    if (this->PAR_Alphabet.find("DNA") != std::string::npos) {
        this->alphabet = new bpp::DNA();
        this->codonAlphabet = false;
    } else if (this->PAR_Alphabet.find("Codon") != std::string::npos){
        this->alphabet = new bpp::DNA();
        // This is additional to the alphabet instance
        this->alphabet = new CodonAlphabet(dynamic_cast<bpp::NucleicAlphabet *>(this->alphabet));
        this->codonAlphabet = true;

        this->codeDesc = ApplicationTools::getStringParameter("genetic_code", this->getParams(), "Standard","", true, true);
        this->gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(dynamic_cast<bpp::CodonAlphabet *>(this->alphabetNoGaps)->getNucleicAlphabet(),this->codeDesc));

    } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
        this->alphabet = new bpp::ProteicAlphabet();
        this->codonAlphabet = false;
    } else {
        exit(EXIT_FAILURE);
    }

}

void CastorApplication::getAlphabet(){

    this->PAR_Alphabet = ApplicationTools::getStringParameter("alphabet", this->getParams(), "DNA", "",true, true);

    // Alphabet without gaps
    this->alphabetNoGaps = bpp::SequenceApplicationTools::getAlphabet(this->getParams(), "", false,false);

    // Alphabet used for all the computational steps (it can allows for gap extension)
    if (this->PAR_model_indels) {

        this->getAlphabetIndel();

    } else {

        this->getAlphabetNoIndel();

    }

}

void CastorApplication::getUnalignedSequences(){

    // If the user requires the computation of an alignment, then the input file is made of unaligned sequences
    bpp::Fasta seqReader;
    this->sequences = seqReader.readSequences(this->PAR_input_sequences, this->alphabet);

}

void CastorApplication::getAlignedSequences(){

    bpp::VectorSiteContainer *allSites = SequenceApplicationTools::getSiteContainer(this->alphabet,this->getParams());

    this->sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, this->getParams(), "", true,!this->PAR_model_indels, true, 1);

    delete allSites;

    AlignmentUtils::checkAlignmentConsistency(*this->sites);

}

void CastorApplication::getData(){

    this->PAR_input_sequences = ApplicationTools::getAFilePath("input.sequence.file", this->getParams(),true, true, "", false, "", 1);

    this->sequences = nullptr;
    this->sites = nullptr;

    try {

        if (this->PAR_alignment) {

            this->getUnalignedSequences();

        } else {

            this->getAlignedSequences();

        }

    } catch (bpp::Exception &e) {
        LOG(FATAL) << "Error when reading sequence file due to: " << e.message();
    }

}

void CastorApplication::getASRV(bpp::SubstitutionModel *smodel){
    // Among site rate variation (ASVR)

    this->rDist = nullptr;

    if (smodel->getNumberOfStates() >= 2 * smodel->getAlphabet()->getSize()) {
        // Markov-modulated Markov model!
        this->rDist = new ConstantRateDistribution();
    } else {
        this->rDist = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());
    }

}

void CastorApplication::computeMSA(bpp::SubstitutionModel *smodel,bpp::Tree *tree,tshlib::Utree *utree,UtreeBppUtils::treemap &tm){

    this->alignment = nullptr;
    this->proPIP = nullptr;

    this->PAR_output_file_msa = ApplicationTools::getAFilePath("output.msa.file", this->getParams(), false,false, "", true, "", 1);

    this->PAR_alignment_version = ApplicationTools::getStringParameter("alignment.version",this->getParams(), "cpu", "",true, 0);

    this->PAR_alignment_sbsolutions = ApplicationTools::getIntParameter("alignment.sb_solutions",this->getParams(), 1, "", true, 0);

    this->PAR_alignment_sbtemperature = ApplicationTools::getDoubleParameter("alignment.sb_temperature",this->getParams(), 1.0, "", true, 0);


    bpp::ApplicationTools::displayResult("Aligner optimised for:", this->PAR_alignment_version);

    bpp::ApplicationTools::displayBooleanResult("Stochastic backtracking active", this->PAR_alignment_sbsolutions > 1);
    if (PAR_alignment_sbsolutions > 1) {
        ApplicationTools::displayResult("Number of stochastic solutions:",TextTools::toString(this->PAR_alignment_sbsolutions));
    }

    DLOG(INFO) << "[Alignment sequences] Starting MSA_t inference using Pro-PIP...";

    // Execute alignment on post-order node list
    std::vector<tshlib::VirtualNode *> ftn = utree->getPostOrderNodeList();


    int num_sb = 1;         // number of sub-optimal MSAs
    double temperature = 0; // temperature for SB version

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
        ApplicationTools::displayError("The user specified an unknown alignment.version. The execution will not continue.");
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    proPIP = new bpp::progressivePIP(utree,               // tshlib tree
                                     tree,                // bpp tree
                                     smodel,              // substitution model
                                     tm,                  // tree-map
                                     this->sequences,           // un-aligned input sequences
                                     this->rDist,               // rate-variation among site distribution
                                     this->getSeed());  // seed for random number generation

    proPIP->_initializePIP(ftn,          // list of tshlib nodes in on post-order (correct order of execution)
                           DPversion,    // version of the alignment algorithm
                           num_sb,       // number of suboptimal MSAs
                           temperature); // to tune the greedyness of the sub-optimal solution

    proPIP->PIPnodeAlign(); // align input sequences with a DP algorithm under PIP

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "\nAlignment elapsed time (msec): " << duration << std::endl;
    ApplicationTools::displayResult("\nAlignment elapsed time (msec):",TextTools::toString((double) duration, 4));

    // convert PIPmsa into a sites objects
    this->sites = PIPmsaUtils::PIPmsa2Sites(proPIP->alphabet_,
                                                *(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getseqNames()),
                                                *(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getMSA()));


    // Export alignment to file
    if (PAR_output_file_msa.find("none") == std::string::npos) {
        DLOG(INFO) << "[Alignment sequences]\t The final alignment can be found in " << PAR_output_file_msa;
        bpp::Fasta seqWriter;
        seqWriter.writeAlignment(TextUtils::appendToFilePath(PAR_output_file_msa, "initial"), *this->sites, true);
    }

    double MSAscore = proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getScore();

    ApplicationTools::displayResult("Alignment log likelihood", TextTools::toString(MSAscore, 15));

}