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

void CastorApplication::getASRV(){
    // Among site rate variation (ASVR)

    this->rDist = nullptr;

    if (this->smodel->getNumberOfStates() >= 2 * smodel->getAlphabet()->getSize()) {
        // Markov-modulated Markov model!
        this->rDist = new ConstantRateDistribution();
    } else {
        this->rDist = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());
    }

}

void CastorApplication::computeMSA(bpp::Tree *tree,tshlib::Utree *utree,UtreeBppUtils::treemap &tm){

    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::vector<tshlib::VirtualNode *> ftn;
    enumDP3Dversion DPversion = CPU;
    double duration = 0.0;

    this->alignment = nullptr;
    this->proPIP = nullptr;

    this->PAR_output_file_msa = ApplicationTools::getAFilePath("output.msa.file", this->getParams(), false,false, "", true, "", 1);

    this->PAR_alignment_version = ApplicationTools::getStringParameter("alignment.version",this->getParams(), "cpu", "",true, 0);

    this->PAR_alignment_sbsolutions = ApplicationTools::getIntParameter("alignment.sb_solutions",this->getParams(), 1, "", true, 0);

    this->PAR_alignment_sbtemperature = ApplicationTools::getDoubleParameter("alignment.sb_temperature",this->getParams(), 1.0, "", true, 0);

    // Execute alignment on post-order node list
    ftn = utree->getPostOrderNodeList();

    // DP3D version
    if (this->PAR_alignment_version.find("cpu") != std::string::npos) {
        DPversion = CPU; // slower but uses less memory
        this->num_sb = 1;
        this->temperature = 0;
    } else if (this->PAR_alignment_version.find("ram") != std::string::npos) {
        DPversion = RAM; // faster but uses more memory
        this->num_sb = 1;
        this->temperature = 0;
    } else if (this->PAR_alignment_version.find("sb") != std::string::npos) {
        DPversion = SB;  // stochastic backtracking version
        this->num_sb = PAR_alignment_sbsolutions;
        this->temperature = PAR_alignment_sbtemperature;
    } else {
        ApplicationTools::displayError("The user specified an unknown alignment.version. The execution will not continue.");
    }

    t1 = std::chrono::high_resolution_clock::now();

    this->proPIP = new bpp::progressivePIP(utree,               // tshlib tree
                                     tree,                // bpp tree
                                     this->smodel,              // substitution model
                                     tm,                  // tree-map
                                     this->sequences,           // un-aligned input sequences
                                     this->rDist,               // rate-variation among site distribution
                                     this->getSeed());  // seed for random number generation

    this->proPIP->_initializePIP(ftn,          // list of tshlib nodes in on post-order (correct order of execution)
                           DPversion,    // version of the alignment algorithm
                           num_sb,       // number of suboptimal MSAs
                           temperature); // to tune the greedyness of the sub-optimal solution

    this->proPIP->PIPnodeAlign(); // align input sequences with a DP algorithm under PIP

    t2 = std::chrono::high_resolution_clock::now();

    //auto
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    bpp::ApplicationTools::displayResult("\nAlignment elapsed time (msec): ",duration);
    bpp::ApplicationTools::displayResult("\nAlignment elapsed time (msec):",TextTools::toString((double) duration, 4));

    // convert PIPmsa into a sites objects
    this->sites = PIPmsaUtils::PIPmsa2Sites(proPIP->alphabet_,
                                                *(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getseqNames()),
                                                *(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getMSA()));

    // Export alignment to file
    if (this->PAR_output_file_msa.find("none") == std::string::npos) {

        bpp::Fasta seqWriter;
        seqWriter.writeAlignment(TextUtils::appendToFilePath(this->PAR_output_file_msa, "initial"), *this->sites, true);

        DLOG(INFO) << "[Alignment sequences]\t The final alignment can be found in " << this->PAR_output_file_msa;

    }

}

void CastorApplication::initLkFun(bpp::Tree *tree,tshlib::Utree *utree,UtreeBppUtils::treemap &tm){

    // Initialization likelihood functions
    this->tl = nullptr;

    // Get transition model from substitution model
    if (!this->PAR_model_indels) {
        this->model = bpp::PhylogeneticsApplicationTools::getTransitionModel(this->alphabet, this->gCode.get(), this->sites,this->getParams(), "", true, false, 0);
    } else {
        std::unique_ptr<bpp::TransitionModel> test;
        test.reset(this->smodel);
        this->model = test.release();
    }

    // Initialise likelihood functions
    if (!this->PAR_model_indels) {

        tl = new bpp::UnifiedTSHomogeneousTreeLikelihood(*tree, *this->sites, this->model, this->rDist, utree, &tm, true,this->getParams(), "", false, false,false);

    } else {

        tl = new bpp::UnifiedTSHomogeneousTreeLikelihood_PIP(*tree, *this->sites, this->model, this->rDist, utree, &tm, true,this->getParams(), "", false,false, false);

    }

    tl->initialize();

}

void CastorApplication::optimizeParameters(){

    this->tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(Optimizators::optimizeParameters(this->tl,this->tl->getParameters(),this->getParams(),"",true,true,0));

    this->logL = this->tl->getLogLikelihood();

}

void CastorApplication::parameterSanityCheck(){

    bpp::ParameterList pl;
    bool f = false;
    bool removeSaturated = false;
    size_t s = 0;

    //Listing parameters
    this->paramNameFile = ApplicationTools::getAFilePath("output.parameter_names.file", this->getParams(), false,false, "", true, "none", 1);

    if (this->paramNameFile != "none") {

        bpp::ApplicationTools::displayResult("List parameters to", this->paramNameFile);

        std::ofstream pnfile(this->paramNameFile.c_str(), ios::out);

        pl = this->tl->getParameters();

        for (size_t i = 0; i < pl.size(); ++i) {
            pnfile << pl[i].getName() << endl;
        }

        pnfile.close();
    }

    //Check initial likelihood:
    this->logL = this->tl->getValue();

    if (std::isinf(this->logL)) {

        // This may be due to null branch lengths, leading to null likelihood!
        bpp::ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
        bpp::ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
        bpp::ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");

        pl = this->tl->getBranchLengthsParameters();

        for (unsigned int i = 0; i < pl.size(); i++) {
            if (pl[i].getValue() < 0.000001){
                pl[i].setValue(0.000001);
            }
        }

        this->tl->matchParametersValues(pl);

        this->logL = this->tl->getLogLikelihood();

    }

    bpp::ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-this->logL, 15));

    if (std::isinf(this->logL)) {

        bpp::ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");

        if (this->codonAlphabet) {

            f = false;

            for (size_t i = 0; i < this->sites->getNumberOfSites(); i++) {

                if (std::isinf(this->tl->getLogLikelihoodForASite(i))) {

                    const Site &site = this->sites->getSite(i);

                    s = site.size();

                    for (size_t j = 0; j < s; j++) {

                        if (this->gCode->isStop(site.getValue(j))) {

                            (*bpp::ApplicationTools::error << "Stop Codon at site " << site.getPosition()<< " in sequence "<< this->sites->getSequence(j).getName()).endLine();

                            f = true;
                        }

                    }

                }

            }

            if (f){
                exit(EXIT_FAILURE);
            }

        }

        removeSaturated = bpp::ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites",this->getParams(), false, "",true, 1);

        if (!removeSaturated) {

            std::ofstream debug("DEBUG_likelihoods.txt", ios::out);

            for (size_t i = 0; i < this->sites->getNumberOfSites(); i++) {
                debug << "Position " << this->sites->getSite(i).getPosition() << " = " << this->tl->getLogLikelihoodForASite(i)<< endl;
            }

            debug.close();

            bpp::ApplicationTools::displayError("!!! Site-specific likelihood have been written in file DEBUG_likelihoods.txt .");
            bpp::ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
            bpp::ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");

            exit(EXIT_FAILURE);

        } else {

            bpp::ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);

            for (size_t i = this->sites->getNumberOfSites(); i > 0; --i) {

                if (std::isinf(this->tl->getLogLikelihoodForASite(i - 1))) {

                    bpp::ApplicationTools::displayResult("Ignore saturated site", this->sites->getSite(i - 1).getPosition());

                    this->sites->deleteSite(i - 1);
                }
            }

            bpp::ApplicationTools::displayResult("Number of sites retained", this->sites->getNumberOfSites());

            this->tl->setData(*this->sites);

            this->tl->initialize();

            this->logL = this->tl->getValue();

            if (std::isinf(this->logL)) {
                throw Exception("Likelihood is still 0 after saturated sites are removed! Looks like a bug...");
            }

            bpp::ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-this->logL, 15));

        }
    }

}

void CastorApplication::output(){

    // Export final tree (if nexus is required, then our re-implementation of the the nexus writer is called)
    bpp::Tree *tree = nullptr;
    bpp::ParameterList parameters;

    tree = new TreeTemplate<Node>(this->tl->getTree());

    this->PAR_output_tree_format = bpp::ApplicationTools::getStringParameter("output.tree.format",this->getParams(),"Newick","",true,true);

    if (this->PAR_output_tree_format.find("Nexus") != std::string::npos) {
        std::vector<Tree *> tmp;
        tmp.push_back(tree);
        OutputUtils::writeNexusMetaTree(tmp, this->getParams());
    } else {
        PhylogeneticsApplicationTools::writeTree(*tree, this->getParams());
    }

    // Export annotation file (tab separated values)
    this->PAR_output_annotation_file = ApplicationTools::getAFilePath("output.annotation.file",this->getParams(),false,false,"",true,"",1);

    if (this->PAR_output_annotation_file.find("none") == std::string::npos) {
        bpp::ApplicationTools::displayResult("Output annotation to file", this->PAR_output_annotation_file);
        OutputUtils::writeTreeAnnotations2TSV(tree, this->PAR_output_annotation_file);
    }

    // Write parameters to screen:
    bpp::ApplicationTools::displayResult("Final Log likelihood", TextTools::toString(this->logL, 15));

    parameters = this->tl->getSubstitutionModelParameters();
    for (size_t i = 0; i < parameters.size(); i++) {
        bpp::ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }

    parameters = this->tl->getRateDistributionParameters();
    for (size_t i = 0; i < parameters.size(); i++) {
        bpp::ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }

    // Checking convergence:
    PhylogeneticsApplicationTools::checkEstimatedParameters(this->tl->getParameters());

    // Write parameters to file (according to arguments)
    OutputUtils::exportOutput(this->tl, this->sites, this->getParams());

}

void CastorApplication::bootstrapping(tshlib::Utree *utree,UtreeBppUtils::treemap &tm){

    // Compute support measures
    this->PAR_support = bpp::ApplicationTools::getStringParameter("support", this->getParams(), "", "", true,true);

    if (this->PAR_support == "bootstrap") {

        bpp::ApplicationTools::displayMessage("\n[Tree support measures]");

        bpp::Bootstrap(this->tl, *this->sites, this->rDist, utree, &tm, this->getParams(), "support.");
    }

}