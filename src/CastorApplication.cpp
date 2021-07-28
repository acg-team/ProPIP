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
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Distance/PGMA.h>

#include "ExtendedAlphabet.hpp"
#include "Utils.hpp"
#include "UnifiedTSHTopologySearch.hpp"
#include "PIP.hpp"
#include "RHomogeneousTreeLikelihood_Generic.hpp"
#include "Optimizators.hpp"
#include "SupportMeasures.hpp"
#include "UnifiedDistanceEstimation.hpp"
#include "inference_indel_rates.hpp"

#include "DistanceFactory.hpp"
#include "DistanceFactoryAngle.hpp"

#include <random>

using namespace bpp;

void CastorApplication::init() {

    this->timerStarted_ = false;
    this->codonAlphabet_ = false;
    this->estimatePIPparameters = false;
    this->computeFrequenciesFromData = false;
    this->PAR_model_indels = false;
    this->PAR_alignment = false;

    this->lambda = 0.0;
    this->mu = 0.0;
    this->logL = 0.0;
    this->temperature = 0.0;
    this->PAR_alignment_sbtemperature = 0.0;
    this->seed_ = 0.0;

    this->num_sb = 0;
    this->PAR_alignment_sbsolutions = 0;

    this->DPversion = CPU;

    this->distances = nullptr;
    this->distMethod = nullptr;
    this->smodel = nullptr;
    this->model = nullptr;
    this->utree = nullptr;
    this->alignment = nullptr;
    this->alphabetNoGaps = nullptr;
    this->alphabet = nullptr;
    this->tl = nullptr;
    this->sequences = nullptr;
    this->sites = nullptr;
    this->rDist = nullptr;
    this->tree = nullptr;
    this->proPIP = nullptr;

}

CastorApplication::CastorApplication(int argc, char *argv[], const std::string &name, const std::string &strVersion, const std::string &build_date) :
        appName_(name), appBuild_(build_date), appVersion_(strVersion), params_(), timerStarted_(false) {


    init();

    bool showversion = false;
    bool noint = false;

    params_ = bpp::AttributesTools::parseOptions(argc, argv);

    showversion = bpp::ApplicationTools::getBooleanParameter("version", params_, false, "", true, 3);

    bpp::ApplicationTools::warningLevel = bpp::ApplicationTools::getIntParameter("warning", params_, 0, "", true, 3);

    noint = bpp::ApplicationTools::getBooleanParameter("noninteractive", params_, false, "", true, 3);

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
        exit(EXIT_SUCCESS);
    }
}

void CastorApplication::version() {
    std::cout << appName_ << std::endl;
    std::cout << appVersion_ << std::endl;
    std::cout << appBuild_ << std::endl;
}

void CastorApplication::banner() {

    auto host_name = boost::asio::ip::host_name();

    bpp::ApplicationTools::displayMessage("------------------------------------------------------------------------------");
    bpp::ApplicationTools::displayMessage(appName_);
    bpp::ApplicationTools::displayMessage("Phylogenetic Tree Inference and Multiple Sequence Alignment under Indel models");
    bpp::ApplicationTools::displayMessage("Authors: Lorenzo Gatti & Massimo Maiolo");
    bpp::ApplicationTools::displayMessage("Build on commit: " + appVersion_);
    bpp::ApplicationTools::displayMessage("On date: "+ appBuild_);
    bpp::ApplicationTools::displayMessage("------------------------------------------------------------------------------");
    bpp::ApplicationTools::displayResult("Execution started on:", host_name);

}

void CastorApplication::help() {
    std::cout << appName_ << std::endl << std::endl;
    std::cout << "Usage: Castor [arguments] or [params=file.txt]" << std::endl;
    std::cout << "Documentation can be found at https://bitbucket.org/lorenzogatti89/castor/" << std::endl;
}

const std::string &CastorApplication::getParam(const std::string &name) const {
    if (params_.find(name) == params_.end()) throw bpp::Exception("BppApplication::getParam(). Parameter '" + name + "' not found.");
    return params_[name];
}

void CastorApplication::startTimer() {
    ApplicationTools::startTimer();
    timerStarted_ = true;
}

void CastorApplication::done() {

    // Delete objects and free memory
    delete this->alphabet;
    delete this->alphabetNoGaps;
    delete this->sites;
    delete this->rDist;
    delete this->tl;
    delete this->sequences;
    delete this->tree;

    DLOG(INFO) << appName_ << "'s done. Bye.";

    if (timerStarted_){
        bpp::ApplicationTools::displayTime("Total execution time:");
    }
}

void CastorApplication::getCLIarguments(){

    this->PAR_model_substitution = bpp::ApplicationTools::getStringParameter("model", this->getParams(), "JC69","", true, true);

    this->PAR_alignment = bpp::ApplicationTools::getBooleanParameter("alignment", this->getParams(), false);

    // Split model string description and test if PIP is required
    bpp::KeyvalTools::parseProcedure(this->PAR_model_substitution, this->modelStringName, this->modelMap);

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
        this->codonAlphabet_ = false;
    } else if(this->PAR_Alphabet.find("Codon") != std::string::npos){
        this->alphabet = new bpp::DNA_EXTENDED();

        // This is additional to the alphabet instance
        this->alphabet = new CodonAlphabet_Extended(dynamic_cast<bpp::NucleicAlphabet *>(this->alphabet));
        this->codonAlphabet_ = true;

        this->codeDesc = bpp::ApplicationTools::getStringParameter("genetic_code", this->getParams(), "Standard","", true, true);
        this->gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(dynamic_cast<bpp::CodonAlphabet *>(this->alphabetNoGaps)->getNucleicAlphabet(),this->codeDesc));

    } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
        this->alphabet = new bpp::ProteicAlphabet_Extended();
        this->codonAlphabet_ = false;
    }

}

void CastorApplication::getAlphabetNoIndel(){

    if (this->PAR_Alphabet.find("DNA") != std::string::npos) {
        this->alphabet = new bpp::DNA();
        this->codonAlphabet_ = false;
    } else if (this->PAR_Alphabet.find("Codon") != std::string::npos){
        this->alphabet = new bpp::DNA();

        // This is additional to the alphabet instance
        this->alphabet = new CodonAlphabet(dynamic_cast<bpp::NucleicAlphabet *>(this->alphabet));
        this->codonAlphabet_ = true;

        this->codeDesc = bpp::ApplicationTools::getStringParameter("genetic_code", this->getParams(), "Standard","", true, true);
        this->gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(dynamic_cast<bpp::CodonAlphabet *>(this->alphabetNoGaps)->getNucleicAlphabet(),this->codeDesc));

    } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
        this->alphabet = new bpp::ProteicAlphabet();
        this->codonAlphabet_ = false;
    } else {
        exit(EXIT_FAILURE);
    }

}

void CastorApplication::getAlphabet(){

    this->PAR_Alphabet = bpp::ApplicationTools::getStringParameter("alphabet", this->getParams(), "DNA", "",true, true);

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
    if (this->smodel->getNumberOfStates() >= 2 * smodel->getAlphabet()->getSize()) {
        // Markov-modulated Markov model!
        this->rDist = new ConstantRateDistribution();
    } else {
        this->rDist = bpp::PhylogeneticsApplicationTools::getRateDistribution(this->getParams());
    }

}

bpp::TransitionModel * CastorApplication::getTransitionModelFromSubsModel(bool PAR_model_indels,
                                                                          bpp::SubstitutionModel *smodel,
                                                                          const Alphabet* alphabet,
                                                                          const GeneticCode* gCode,
                                                                          const SiteContainer* data,
                                                                          std::map<std::string, std::string>& params,
                                                                          const std::string& suffix,
                                                                          bool suffixIsOptional,
                                                                          bool verbose,
                                                                          int warn){

    bpp::TransitionModel *transition_model = nullptr;

    // Get transition model from substitution model
    if (!PAR_model_indels) {
        transition_model = bpp::PhylogeneticsApplicationTools::getTransitionModel(alphabet,gCode,data,params,suffix,suffixIsOptional,verbose,warn);
    } else {
        std::unique_ptr<bpp::TransitionModel> test;
        test.reset(smodel);
        transition_model = test.release();
    }

    return transition_model;
}

void CastorApplication::computeMSA(){

    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::vector<tshlib::VirtualNode *> ftn;
    double duration = 0.0;

    this->PAR_output_file_msa = bpp::ApplicationTools::getAFilePath("output.msa.file", this->getParams(), false,false, "", true, "", 1);

    this->PAR_alignment_version = bpp::ApplicationTools::getStringParameter("alignment.version",this->getParams(), "cpu", "",true, 0);

    this->PAR_alignment_sbsolutions = bpp::ApplicationTools::getIntParameter("alignment.sb_solutions",this->getParams(), 1, "", true, 0);

    this->PAR_alignment_sbtemperature = bpp::ApplicationTools::getDoubleParameter("alignment.sb_temperature",this->getParams(), 1.0, "", true, 0);

    // Execute alignment on post-order node list
    ftn = utree->getPostOrderNodeList();

    // DP3D version
    if (this->PAR_alignment_version.find("cpu") != std::string::npos) {
        this->DPversion = CPU; // slower but uses less memory
        this->num_sb = 1;
        this->temperature = 0;
    } else if (this->PAR_alignment_version.find("ram") != std::string::npos) {
        this->DPversion = RAM; // faster but uses more memory
        this->num_sb = 1;
        this->temperature = 0;
    } else if (this->PAR_alignment_version.find("sb") != std::string::npos) {
        this->DPversion = SB;  // stochastic backtracking version
        this->num_sb = this->PAR_alignment_sbsolutions;
        this->temperature = this->PAR_alignment_sbtemperature;
    } else {
        bpp::ApplicationTools::displayError("The user specified an unknown alignment.version. The execution will not continue.");
    }

    t1 = std::chrono::high_resolution_clock::now();

    this->proPIP = new bpp::progressivePIP(this->utree,               // tshlib tree
                                           this->tree,                // bpp tree
                                           this->smodel,              // substitution model
                                           this->tm,                  // tree-map
                                           this->sequences,           // un-aligned input sequences
                                           this->rDist,               // rate-variation among site distribution
                                           this->getSeed());          // seed for random number generation

    this->proPIP->_initializePIP(ftn,          // list of tshlib nodes in on post-order (correct order of execution)
                                 this->DPversion,    // version of the alignment algorithm
                                 this->num_sb,       // number of suboptimal MSAs
                                 this->temperature); // to tune the greedyness of the sub-optimal solution

    this->proPIP->PIPnodeAlign(); // align input sequences with a DP algorithm under PIP

    t2 = std::chrono::high_resolution_clock::now();

    //auto
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    bpp::ApplicationTools::displayResult("\nAlignment elapsed time (msec): ",duration);
    bpp::ApplicationTools::displayResult("\nAlignment elapsed time (msec):",bpp::TextTools::toString((double) duration, 4));

    // convert PIPmsa into a sites objects
    this->sites = PIPmsaUtils::PIPmsa2Sites(this->proPIP->alphabet_,
                                            *(this->proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getseqNames()),
                                            *(this->proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getMSA()));

    // Export alignment to file
    if (this->PAR_output_file_msa.find("none") == std::string::npos) {

        bpp::Fasta seqWriter;
        seqWriter.writeAlignment(TextUtils::appendToFilePath(this->PAR_output_file_msa, "initial"), *this->sites, true);

        DLOG(INFO) << "[Alignment sequences]\t The final alignment can be found in " << this->PAR_output_file_msa;

    }

}

void CastorApplication::initLkFun(){

    // Initialization likelihood functions
    this->tl = nullptr;

    // Get transition model from substitution model
    this->model=getTransitionModelFromSubsModel(this->PAR_model_indels,
                                                this->smodel,
                                                this->alphabet,
                                                this->gCode.get(),
                                                this->sites,
                                                this->getParams(),
                                                "",
                                                true,
                                                false,
                                                0);

    // Initialise likelihood functions
    if (!this->PAR_model_indels) {

        this->tl = new bpp::UnifiedTSHomogeneousTreeLikelihood(*this->tree, *this->sites, this->model, this->rDist, this->utree, &this->tm, true,this->getParams(), "", false, false,false);

    } else {

        this->tl = new bpp::UnifiedTSHomogeneousTreeLikelihood_PIP(*this->tree, *this->sites, this->model, this->rDist, this->utree, &this->tm, true,this->getParams(), "", false,false, false);

    }

    this->tl->initialize();

}

void CastorApplication::optimizeParameters(){

    this->tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(Optimizators::optimizeParameters(this->tl,this->tl->getParameters(),this->getParams(),this->tm,this->sites,"",true,true,0));

    this->logL = this->tl->getLogLikelihood();

}

bpp::ParameterList CastorApplication::getParametersList(){

    bpp::ParameterList pl;

    //Listing parameters
    this->paramNameFile = bpp::ApplicationTools::getAFilePath("output.parameter_names.file", this->getParams(), false,false, "", true, "none", 1);

    if (this->paramNameFile != "none") {

        bpp::ApplicationTools::displayResult("List parameters to", this->paramNameFile);

        std::ofstream pnfile(this->paramNameFile.c_str(), ios::out);

        pl = this->tl->getParameters();

        for (size_t i = 0; i < pl.size(); ++i) {
            pnfile << pl[i].getName() << endl;
        }

        pnfile.close();
    }

    return pl;
}

double CastorApplication::checkLkValue(bpp::ParameterList &pl){

    double logL = 0.0;

    //Check initial likelihood:
    logL = this->tl->getValue();

    if (std::isinf(logL)) {

        // This may be due to null branch lengths, leading to null likelihood!
        bpp::ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
        bpp::ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
        bpp::ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");

        pl = this->tl->getBranchLengthsParameters();

        for (unsigned int i = 0; i < pl.size(); i++) {
            if (pl[i].getValue() < MIN_BRANCH_LEN){
                pl[i].setValue(MIN_BRANCH_LEN);
            }
        }

        this->tl->matchParametersValues(pl);

        logL = this->tl->getLogLikelihood();

    }

    return logL;
}

void CastorApplication::checkStopCodon(){

    bool f = false;
    size_t s = 0;

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

void CastorApplication::removeSaturatedSite(){

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

    bpp::ApplicationTools::displayResult("Initial log likelihood", bpp::TextTools::toString(-this->logL, 15));

}

void CastorApplication::resolveZeroLKValue(){

    bool removeSaturated = false;

    bpp::ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");

    if (this->codonAlphabet_) {

        this->checkStopCodon();

    }

    removeSaturated = bpp::ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites",this->getParams(), false, "",true, 1);

    if (removeSaturated) {

        this->removeSaturatedSite();

    } else {

        std::ofstream debug("DEBUG_likelihoods.txt", ios::out);

        for (size_t i = 0; i < this->sites->getNumberOfSites(); i++) {
            debug << "Position " << this->sites->getSite(i).getPosition() << " = " << this->tl->getLogLikelihoodForASite(i)<< endl;
        }

        debug.close();

        bpp::ApplicationTools::displayError("!!! Site-specific likelihood have been written in file DEBUG_likelihoods.txt .");
        bpp::ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
        bpp::ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");

        exit(EXIT_FAILURE);

    }

}

void CastorApplication::parameterSanityCheck(){

    bpp::ParameterList pl;

    //Listing parameters
    pl=getParametersList();

    //Check initial likelihood
    this->logL=checkLkValue(pl);

    bpp::ApplicationTools::displayResult("Initial log likelihood", bpp::TextTools::toString(-this->logL, 15));

    if (std::isinf(this->logL)) {
        this->resolveZeroLKValue();
    }

}

void CastorApplication::output(){

    // Export final tree (if nexus is required, then our re-implementation of the the nexus writer is called)
    bpp::Tree *local_tree = nullptr;
    bpp::ParameterList local_parameters;

    local_tree = new TreeTemplate<Node>(this->tl->getTree());

    this->PAR_output_tree_format = bpp::ApplicationTools::getStringParameter("output.tree.format",this->getParams(),"Newick","",true,true);

    if (this->PAR_output_tree_format.find("Nexus") != std::string::npos) {
        std::vector<Tree *> tmp_vector_tree;
        tmp_vector_tree.push_back(local_tree);
        OutputUtils::writeNexusMetaTree(tmp_vector_tree, this->getParams());
    } else {
        bpp::PhylogeneticsApplicationTools::writeTree(*local_tree, this->getParams());
    }

    // Export annotation file (tab separated values)
    this->PAR_output_annotation_file = bpp::ApplicationTools::getAFilePath("output.annotation.file",this->getParams(),false,false,"",true,"",1);

    if (this->PAR_output_annotation_file.find("none") == std::string::npos) {
        bpp::ApplicationTools::displayResult("Output annotation to file", this->PAR_output_annotation_file);
        OutputUtils::writeTreeAnnotations2TSV(local_tree, this->PAR_output_annotation_file);
    }

    // Write local_parameters to screen:
    bpp::ApplicationTools::displayResult("Final Log likelihood", bpp::TextTools::toString(this->logL, 15));

    local_parameters = this->tl->getSubstitutionModelParameters();
    for (size_t i = 0; i < local_parameters.size(); i++) {
        bpp::ApplicationTools::displayResult(local_parameters[i].getName(), bpp::TextTools::toString(local_parameters[i].getValue()));
    }

    local_parameters = this->tl->getRateDistributionParameters();
    for (size_t i = 0; i < local_parameters.size(); i++) {
        bpp::ApplicationTools::displayResult(local_parameters[i].getName(), bpp::TextTools::toString(local_parameters[i].getValue()));
    }

    // Checking convergence:
    bpp::PhylogeneticsApplicationTools::checkEstimatedParameters(this->tl->getParameters());

    // Write local_parameters to file (according to arguments)
    OutputUtils::exportOutput(this->tl, this->sites, this->getParams());

}

void CastorApplication::bootstrapping(){

    // Compute support measures
    this->PAR_support = bpp::ApplicationTools::getStringParameter("support", this->getParams(), "", "", true,true);

    if (this->PAR_support == "bootstrap") {

        bpp::ApplicationTools::displayMessage("\n[Tree support measures]");

        bpp::Bootstrap(this->tl, *this->sites, this->rDist, this->utree, &this->tm, this->getParams(), "support.");
    }

}


void CastorApplication::initCanonicalSubstitutionModel(){

    // Instantiation of the canonical substitution model
    if (this->PAR_Alphabet.find("Codon") != std::string::npos) {

        this->smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabetNoGaps,this->gCode.get(), this->sites,this->modelMap, "", true,false, 0);

    }else if (this->PAR_Alphabet.find("Protein") != std::string::npos){

        this->smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabetNoGaps, this->gCode.get(), this->sites,this->modelMap, "", true, false, 0);

    } else {

        this->smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabet, this->gCode.get(), this->sites,this->modelMap, "", true, false, 0);

    }

}

void CastorApplication::extendSubstitutionModelWithPIP(const SequenceContainer *data){

    if (this->PAR_Alphabet.find("DNA") != std::string::npos && this->PAR_Alphabet.find("Codon") == std::string::npos) {

        this->smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(this->alphabet), this->smodel, *data, this->lambda, this->mu,computeFrequenciesFromData);

    } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {

        this->smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(this->alphabet), this->smodel, *data, this->lambda, this->mu,computeFrequenciesFromData);

    } else if (this->PAR_Alphabet.find("Codon") != std::string::npos) {

        this->smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(this->alphabet), this->gCode.get(), this->smodel,*data, lambda, mu,computeFrequenciesFromData);

        bpp::ApplicationTools::displayWarning("Codon models are experimental in the current version... use with caution!");

        DLOG(WARNING) << "CODONS activated byt the program is not fully tested under these settings!";
    }

}

void CastorApplication::getIndelRates(){

    // If PIP, then check if lambda/mu initial values are estimated from the data
    if (this->modelMap.find("estimated") != this->modelMap.end()) {

        bpp::ApplicationTools::displayError("The use of the tag [observed] is obsolete. Use the tag [initFreqs] instead");

        exit(EXIT_FAILURE);

    } else if (this->modelMap.find("initFreqs") != this->modelMap.end()) {

        if (this->modelMap["initFreqs"] == "observed") {

            this->estimatePIPparameters = true;

        } else{

            bpp::ApplicationTools::displayError("Unknow option");

            exit(EXIT_FAILURE);
        }

    } else if( (this->modelMap.find("lambda") == this->modelMap.end()) || (this->modelMap.find("mu") == this->modelMap.end())) {

        this->estimatePIPparameters = true;

    }

    if (this->estimatePIPparameters) {

        inference_indel_rates::infere_indel_rates_from_sequences(this->PAR_input_sequences,
                                                                 this->PAR_Alphabet,
                                                                 this->PAR_alignment,
                                                                 this->PAR_model_indels,
                                                                 this->getParams(),
                                                                 this->tree,
                                                                 this->lambda,
                                                                 this->mu,
                                                                 this->gCode.get(),
                                                                 this->modelMap);

    } else {

        this->lambda = (this->modelMap.find("lambda") == this->modelMap.end()) ? 0.1 : std::stod(this->modelMap["lambda"]);
        this->mu = (this->modelMap.find("mu") == this->modelMap.end()) ? 0.2 : std::stod(this->modelMap["mu"]);

    }

    DLOG(INFO) << "[PIP model] Fixed PIP parameters to (lambda=" << this->lambda << ",mu=" << this->mu << "," "I="<< this->lambda * this->mu << ")";

}

void CastorApplication::getBackgroundFrequencies(){

    std::vector<std::string> keys;

    this->computeFrequenciesFromData = false;

    // If frequencies are estimated from the data, but there is no alignment, then flag it.
    bpp::KeyvalTools::parseProcedure(this->modelMap["model"], this->baseModel, this->basemodelMap);

    for (auto it = this->basemodelMap.begin(); it != this->basemodelMap.end(); ++it){
        keys.push_back(it->first);
    }

    if (!keys.empty()) {
        this->baseModel += "(";

        for (auto &key:keys) {

            if (key != "initFreqs") {
                this->baseModel += key + "=" + this->basemodelMap[key];
            } else {
                if (this->basemodelMap[key] == "observed") {
                    this->computeFrequenciesFromData = true;
                }
            }

            this->baseModel += ",";
        }

        this->baseModel.pop_back();
        this->baseModel += ")";
        this->modelMap["model"] = this->baseModel;
    }

}

void CastorApplication::getSubstitutionIndelModel(){

    const SequenceContainer *data;

    this->getBackgroundFrequencies();

    // Instantiation of the canonical substitution model
    this->initCanonicalSubstitutionModel();

    this->getIndelRates();

    // Instantiate the corrisponding PIP model given the alphabet
    if (this->PAR_alignment) {
        data = this->sequences;
    } else {
        data = this->sites;
    }

    this->extendSubstitutionModelWithPIP(data);

}

void CastorApplication::getSubstitutionNoIndelModel(){

    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*this->sites);

    this->smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabet, this->gCode.get(), this->sites,this->getParams(), "", true, false, 0);

}

void CastorApplication::getSubstitutionModel(){

    bpp::StdStr s1;
    ParameterList local_parameters;

    // Instantiate a substitution model and extend it with PIP
    if (this->PAR_model_indels) {
        this->getSubstitutionIndelModel();
    } else {
        this->getSubstitutionNoIndelModel();
    }

    bpp::ApplicationTools::displayResult("Substitution model", this->smodel->getName());

    if (this->PAR_model_indels){
        bpp::ApplicationTools::displayResult("Indel parameter initial value",(this->estimatePIPparameters) ? "estimated" : "fixed");
    }

    local_parameters = this->smodel->getParameters();

    for (size_t i = 0; i < local_parameters.size(); i++) {
        bpp::ApplicationTools::displayResult(local_parameters[i].getName(), bpp::TextTools::toString(local_parameters[i].getValue()));
    }

    for (size_t i = 0; i < this->smodel->getFrequencies().size(); i++) {
        bpp::ApplicationTools::displayResult("eq.freq(" + this->smodel->getAlphabet()->getName(i) + ")",bpp::TextTools::toString(this->smodel->getFrequencies()[i], 4));
    }

    DLOG(INFO) << "[Substitution model] Number of states: " << (int) this->smodel->getNumberOfStates();

    bpp::PhylogeneticsApplicationTools::printParameters(this->smodel, s1, 1, true);

    DLOG(INFO) << s1.str();

}

void CastorApplication::initTreeMethodUser(){

    this->tree = bpp::PhylogeneticsApplicationTools::getTree(this->getParams());

    DLOG(INFO) << "[Input tree parser] Number of leaves" << this->tree->getNumberOfLeaves();

}

void CastorApplication::initTreeMethodRandom(){

    std::vector<std::string> names;

    names = this->sites->getSequencesNames();

    this->tree = TreeTemplateTools::getRandomTree(names);

}

void CastorApplication::initTreeMethodDistanceWPGMA(){

    auto *wpgma = new PGMA(true);

    this->distMethod = wpgma;

}

void CastorApplication::initTreeMethodDistanceUPGMA(){

    auto *upgma = new PGMA(false);

    this->distMethod = upgma;

}

void CastorApplication::initTreeMethodDistanceNJ(){

    auto *nj = new NeighborJoining();

    nj->outputPositiveLengths(true);

    this->distMethod = nj;

}

void CastorApplication::initTreeMethodDistanceBIONJ(){

    auto *bionj = new BioNJ();

    bionj->outputPositiveLengths(true);

    this->distMethod = bionj;

}

void CastorApplication::initTreeMethodDistanceDistmatrix(){

    // Use a distance matrix provided by the user
    bpp::ApplicationTools::displayResult("Initial tree method", std::string("LZ compression"));

    try {
        this->PAR_distance_matrix = bpp::ApplicationTools::getAFilePath("init.distance.matrix.file",this->getParams(), true, true, "", false,"",0);
    } catch (bpp::Exception &e) {
        LOG(FATAL) << "Error when reading distance matrix file: " << e.message();
    }

    DLOG(INFO) << "initial tree method from LZ compression from matrix file" << this->PAR_distance_matrix;

    this->distances = InputUtils::parseDistanceMatrix(this->PAR_distance_matrix);

    bpp::BioNJ bionj(*this->distances, true, true, false);

    this->tree = bionj.getTree();

}

void CastorApplication::initTreeMethodDistanceInfereDistanceMatrix(){

    int ALPHABET_DIM = 0;
    int K = 0;
    bool mldist_flag = true;
    bool mldist_gap_flag = false;
    double cutoff_dist = 1.0;
    double indel_rate = 1.0;
    bpp::DistanceMethod *local_distMethod = nullptr;
    DistanceFactoryPrographMSA::DistanceFactory *dist_factory_angle = nullptr;
    bpp::DistanceMatrix *local_dist = nullptr;

    auto *bionj = new BioNJ(true, true, false);

    if (this->PAR_Alphabet.find("DNA") != std::string::npos) {
        ALPHABET_DIM = 4;
        K = 6;
    } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
        ALPHABET_DIM = 20;
        K = 2;
    } else if (this->PAR_Alphabet.find("Codon") != std::string::npos) {
        ALPHABET_DIM = 60;
        K = 2;
    }

    if (!this->sequences) {
        bpp::Fasta seqReader;
        this->sequences = seqReader.readSequences(this->PAR_input_sequences, this->alphabet);
    }

    dist_factory_angle = new DistanceFactoryPrographMSA::DistanceFactoryAngle(ALPHABET_DIM, K);

    DistanceFactoryPrographMSA::DistanceMatrix dist_ml = dist_factory_angle->computePwDistances(this->sequences,ALPHABET_DIM,K,mldist_flag,mldist_gap_flag,cutoff_dist,indel_rate);

    local_dist = new DistanceMatrix(this->sequences->getSequencesNames());

    for (int iii = 0; iii < this->sequences->getNumberOfSequences(); iii++) {
        for (int jjj = 0; jjj < this->sequences->getNumberOfSequences(); jjj++) {
            (*local_dist)(iii, jjj) = (abs(dist_ml.distances(iii, jjj)) < DISTCUTOFF ? 0.0 : abs(
                    dist_ml.distances(iii, jjj)));
        }
    }

    bionj->outputPositiveLengths(true);

    local_distMethod = bionj;
    local_distMethod->setDistanceMatrix(*local_dist);
    local_distMethod->computeTree();

    this->tree = local_distMethod->getTree();

    delete dist_factory_angle;
    delete local_dist;
    delete bionj;

}

void CastorApplication::infereDistanceTreeFast(bpp::TransitionModel *local_dmodel,
                                               bpp::VectorSiteContainer *local_sitesDistMethod,
                                               DiscreteDistribution *local_rDist){

    DistanceMatrix *local_distance_matrix = nullptr;

    UnifiedDistanceEstimation local_distEstimation(local_dmodel, local_rDist, local_sitesDistMethod, 1, false);

    local_distEstimation.computeMatrix();

    local_distance_matrix = local_distEstimation.getMatrix();

    this->distMethod->setDistanceMatrix((*local_distance_matrix));

    this->distMethod->computeTree();

    this->tree = distMethod->getTree();

}

void CastorApplication::infereDistanceTreeML(bpp::TransitionModel *local_dmodel,
                                             DiscreteDistribution *local_rDist,
                                             bpp::VectorSiteContainer *local_sitesDistMethod){

    UnifiedDistanceEstimation distEstimation(local_dmodel, local_rDist, local_sitesDistMethod, 1, false);

    this->PAR_optim_distance = ApplicationTools::getStringParameter("init.distance.optimization.method", this->getParams(),"init");

    bpp::ApplicationTools::displayResult("Initial tree model parameters estimation method",this->PAR_optim_distance);

    if (this->PAR_optim_distance == "init") this->PAR_optim_distance = Optimizators::DISTANCEMETHOD_INIT;
    else if (this->PAR_optim_distance == "pairwise")
        this->PAR_optim_distance = Optimizators::DISTANCEMETHOD_PAIRWISE;
    else if (this->PAR_optim_distance == "iterations")
        this->PAR_optim_distance = Optimizators::DISTANCEMETHOD_ITERATIONS;
    else throw Exception("Unknown parameter estimation procedure '" + this->PAR_optim_distance + "'.");

    // Optimisation method verbosity
    auto optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose",this->getParams(), 2);
    string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", this->getParams(),
                                                   false, false);
    auto *messenger = (mhPath == "none") ? nullptr : (mhPath == "std") ? bpp::ApplicationTools::message.get(): new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));

    bpp::ApplicationTools::displayResult("Initial tree optimization handler", mhPath);

    // Optimisation method profiler
    string prPath = ApplicationTools::getAFilePath("optimization.profiler", this->getParams(), false,false);

    auto *profiler = (prPath == "none") ? nullptr : (prPath == "std") ? bpp::ApplicationTools::message.get(): new StlOutputStream(new ofstream(prPath.c_str(), ios::out));

    if (profiler) profiler->setPrecision(20);
    bpp::ApplicationTools::displayResult("Initial tree optimization profiler", prPath);

    // Should I ignore some parameters?
    ParameterList allParameters = local_dmodel->getParameters();
    allParameters.addParameters(local_rDist->getParameters());

    ParameterList parametersToIgnore;
    string paramListDesc = ApplicationTools::getStringParameter(
            "init.distance.optimization.ignore_parameter", this->getParams(),
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
                                                                  this->getParams(), 1000000);
    ApplicationTools::displayResult("Initial tree optimization | max # ML evaluations",
                                    TextTools::toString(nbEvalMax));

    double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance",
                                                            this->getParams(), .000001);
    ApplicationTools::displayResult("Initial tree optimization | Tolerance",
                                    TextTools::toString(tolerance));

    //Here it is:
    this->tree = Optimizators::buildDistanceTreeGeneric(distEstimation, *distMethod, parametersToIgnore,
                                                        !ignoreBrLen, PAR_optim_distance,
                                                        tolerance, nbEvalMax, profiler, messenger,
                                                        optVerbose);

}

void CastorApplication::removeGaps(bpp::VectorSiteContainer *local_sitesDistMethod){

    // Remove gap characters since we are roughly estimating the initial topology
    if (!this->PAR_model_indels) {
        bpp::SiteContainerTools::changeGapsToUnknownCharacters(*local_sitesDistMethod);
    }

}

void CastorApplication::addASRVdistribution(DiscreteDistribution *local_rDist,bpp::TransitionModel *local_dmodel) {

    if (local_dmodel->getNumberOfStates() > local_dmodel->getAlphabet()->getSize()) {

        //Markov-modulated Markov model!
        local_rDist = new ConstantRateDistribution();

    } else {
        local_rDist = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());
    }

}

void CastorApplication::getUsersIndelRates(){

    this->lambda = std::stod(this->modelMap["lambda"]);
    this->mu = std::stod(this->modelMap["mu"]);

}

void CastorApplication::infereDistanceTree(){

    // Compute bioNJ tree using the GTR model
    std::map<std::string, std::string> local_parmap;
    bpp::VectorSiteContainer *local_allSites = nullptr;
    bpp::VectorSiteContainer *local_sitesDistMethod = nullptr;
    bpp::Alphabet *local_alphabetDistMethod = nullptr;
    bpp::SubstitutionModel *local_smodel = nullptr;
    bpp::TransitionModel *local_dmodel = nullptr;
    DiscreteDistribution *local_rDist = nullptr;

    if (true) {

        if (this->PAR_model_indels) {
            local_alphabetDistMethod = this->alphabet;
            local_parmap["model"] = this->modelMap["model"];
        } else {
            local_alphabetDistMethod = this->alphabetNoGaps;
            local_parmap["model"] = this->getParams()["model"];
        }

        local_allSites = SequenceApplicationTools::getSiteContainer(local_alphabetDistMethod, this->getParams());
        local_sitesDistMethod = SequenceApplicationTools::getSitesToAnalyse(*local_allSites, this->getParams());

        if (this->PAR_model_indels) {

            // Instantiation of the canonical substitution model
            if (this->PAR_Alphabet.find("Codon") != std::string::npos ||
                this->PAR_Alphabet.find("Protein") != std::string::npos) {
                local_smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabetNoGaps, this->gCode.get(),
                                                                                  local_sitesDistMethod, this->modelMap, "",
                                                                                  true,
                                                                                  false, 0);
            } else {
                local_smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(local_alphabetDistMethod,
                                                                                  this->gCode.get(), local_sitesDistMethod,
                                                                                  this->modelMap,
                                                                                  "", true,
                                                                                  false, 0);
            }

            if(this->modelMap.find("lambda") == this->modelMap.end() || this->modelMap.find("mu") == this->modelMap.end()){

                if(!this->tree){

                    bpp::SubstitutionModel *smodel_tmp = nullptr;
                    bpp::SubstitutionModel *smodel_copy = nullptr;
                    VectorSiteContainer *sitesDistMethod_tmp = nullptr;
                    bpp::Alphabet *alphabetDistMethod_tmp = nullptr;
                    TransitionModel *dmodel_tmp = nullptr;
                    DiscreteDistribution *rDist_tmp = nullptr;
                    bpp::DistanceMatrix *dm_tmp = nullptr;
                    bpp::DistanceMethod *distMethod_tmp = nullptr;

                    double lambda_tmp = 10.0;
                    double mu_tmp = 0.1;

                    auto *bionj_tmp = new BioNJ(true, true, false);

                    smodel_copy = local_smodel->clone();

                    sitesDistMethod_tmp = local_sitesDistMethod->clone();
                    alphabetDistMethod_tmp = local_alphabetDistMethod->clone();

                    // Instatiate the corrisponding PIP model given the alphabet
                    if (this->PAR_Alphabet.find("DNA") != std::string::npos &&
                        this->PAR_Alphabet.find("Codon") == std::string::npos) {
                        smodel_tmp = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabetDistMethod_tmp), smodel_copy,
                                                 *sitesDistMethod_tmp, lambda_tmp, mu_tmp, false);
                    } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
                        smodel_tmp = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabetDistMethod_tmp), smodel_copy,
                                                *sitesDistMethod_tmp, lambda_tmp, mu_tmp, false);
                    } else if (this->PAR_Alphabet.find("Codon") != std::string::npos) {
                        smodel_tmp = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabetDistMethod_tmp), this->gCode.get(),
                                                   smodel_copy, *sitesDistMethod_tmp,lambda_tmp,mu_tmp, false);
                        ApplicationTools::displayWarning(
                                "Codon models are experimental in the current version... use with caution!");
                        DLOG(WARNING) << "CODONS activated but the program is not fully tested under these settings!";
                    }

                    // Get transition model from substitution model
                    if (!this->PAR_model_indels) {
                        dmodel_tmp = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod_tmp,
                                                                                       this->gCode.get(),
                                                                                       sitesDistMethod_tmp,
                                                                                       local_parmap);
                    } else {
                        unique_ptr<TransitionModel> test;
                        test.reset(smodel_tmp);
                        dmodel_tmp = test.release();
                    }

                    // Add a ASRV distribution
                    if (dmodel_tmp->getNumberOfStates() > dmodel_tmp->getAlphabet()->getSize()) {
                        //Markov-modulated Markov model!
                        rDist_tmp = new ConstantRateDistribution();
                    } else {
                        rDist_tmp = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());
                    }

                    // Remove gap characters since we are roughly estimating the initial topology
                    if (!this->PAR_model_indels) {
                        bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sitesDistMethod_tmp);
                    }

                    UnifiedDistanceEstimation distEstimation_tmp(dmodel_tmp, rDist_tmp, sitesDistMethod_tmp, 1, false);

                    dm_tmp = bpp::SiteContainerTools::computeSimilarityMatrix(*sitesDistMethod_tmp,true,"no full gap",true);

                    bionj_tmp->outputPositiveLengths(true);
                    distMethod_tmp = bionj_tmp;
                    distMethod_tmp->setDistanceMatrix((*dm_tmp));
                    distMethod_tmp->computeTree();
                    tree = distMethod_tmp->getTree();


                }

                inference_indel_rates::infere_indel_rates_from_sequences(this->PAR_input_sequences,this->PAR_Alphabet,this->PAR_alignment,this->PAR_model_indels,this->getParams(),this->tree,this->lambda,this->mu,this->gCode.get(),this->modelMap);

            }else{
                this->getUsersIndelRates();
            }

            // Instatiate the corrisponding PIP model given the alphabet
            if (this->PAR_Alphabet.find("DNA") != std::string::npos){
                //&&this->PAR_Alphabet.find("Codon") == std::string::npos) {

                local_smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(local_alphabetDistMethod),
                                           local_smodel,
                                           *local_sitesDistMethod,
                                           this->lambda,
                                           this->mu,
                                           false);

            } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {

                local_smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(local_alphabetDistMethod),
                                          local_smodel,
                                          *local_sitesDistMethod,
                                          this->lambda,
                                          this->mu,
                                          false);

            } else if (this->PAR_Alphabet.find("Codon") != std::string::npos) {

                local_smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(local_alphabetDistMethod),
                                             this->gCode.get(),
                                             local_smodel,
                                             *local_sitesDistMethod,
                                             this->lambda,
                                             this->mu,
                                             false);

                bpp::ApplicationTools::displayWarning("Codon models are experimental in the current version... use with caution!");

                DLOG(WARNING) << "CODONS activated but the program is not fully tested under these settings!";
            }else{
                throw Exception("Unknown option");
            }

        }

        // Get transition model from substitution model
        if (!this->PAR_model_indels) {
            local_dmodel = PhylogeneticsApplicationTools::getTransitionModel(local_alphabetDistMethod,
                                                                              this->gCode.get(),
                                                                              local_sitesDistMethod,
                                                                              local_parmap);
        } else {
            unique_ptr<TransitionModel> test;
            test.reset(local_smodel);
            local_dmodel = test.release();
        }

        // Add a ASRV distribution
        if (local_dmodel->getNumberOfStates() > local_dmodel->getAlphabet()->getSize()) {
            //Markov-modulated Markov model!
            local_rDist = new ConstantRateDistribution();
        } else {
            local_rDist = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());
        }

        this->removeGaps(local_sitesDistMethod);

        if (this->PAR_distance_method.find("-ml") != std::string::npos) {

            this->infereDistanceTreeML(local_dmodel,local_rDist,local_sitesDistMethod);

        } else {

            if(!tree) {
                // Fast but rough estimate of the initial tree topology
                // (distance based without optimisation -ML)
                this->infereDistanceTreeFast(local_dmodel,local_sitesDistMethod,local_rDist);
            }

        }

        delete local_sitesDistMethod;
        delete distMethod;
        delete local_allSites;

    }

}

void CastorApplication::initTreeMethodDistance(){

    std::string token;

    this->PAR_distance_method = bpp::ApplicationTools::getStringParameter("init.distance.method",this->getParams(), "nj");

    bpp::ApplicationTools::displayResult("Initial tree reconstruction method", this->PAR_distance_method);

    token = this->PAR_distance_method.substr(0, this->PAR_distance_method.find("-"));

    this->distances = nullptr;
    this->distMethod = nullptr;

    if (token == "wpgma") {
        this->initTreeMethodDistanceWPGMA();
    } else if (token == "upgma") {
        this->initTreeMethodDistanceUPGMA();
    } else if (token == "nj") {
        this->initTreeMethodDistanceNJ();
    } else if (token == "bionj") {
        this->initTreeMethodDistanceBIONJ();
    } else if (token == "distmatrix") {
        this->initTreeMethodDistanceDistmatrix();
    } else if (token == "infere_distance_matrix") {
        this->initTreeMethodDistanceInfereDistanceMatrix();
    } else{
        throw Exception("Unknown tree reconstruction method.");
    }

    this->infereDistanceTree();

}

void CastorApplication::initTreeMethod(){

    if (initTreeOpt == "user") {

        this->initTreeMethodUser();

    } else if (initTreeOpt == "random") {

        this->initTreeMethodRandom();

    } else if (this->initTreeOpt == "distance") {

        this->initTreeMethodDistance();

    } else{
        throw Exception("Unknown init tree method.");
    }

}

void CastorApplication::resolveMultifurcations(){

    // If the tree has multifurcation, then resolve it with midpoint rooting
    TreeTemplate<Node> *ttree_=nullptr;

    ttree_ = new TreeTemplate<Node>(*this->tree);

    if (ttree_->getRootNode()->getNumberOfSons() > 2) {
        TreeTemplateTools::midRoot(*(ttree_), TreeTemplateTools::MIDROOT_VARIANCE, false);
        this->tree = ttree_;
    }

}

void CastorApplication::renameTreeNodes(){

    // Rename internal nodes with standard Vxx * where xx is a progressive number
    this->tree->setNodeName(tree->getRootId(), "root");

    UtreeBppUtils::renameInternalNodes(tree);

}

void CastorApplication::initTreeBranchLengthExponentialDistribution(std::map<std::string, std::string> &cmdArgs){

    double value = 0.0;

    value = bpp::ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, 2);

    if (value <= 0){
        throw Exception("Value for branch length must be greater than 0");
    }

    bpp:: ApplicationTools::displayResult("Random branch lengths distributed accordint to exponential distribution with mean ", value);

    double const exp_dist_mean   = value;
    double const exp_dist_lambda = 1 / exp_dist_mean;
    std::random_device rd;
    std::exponential_distribution<> rng (exp_dist_lambda);
    std::mt19937 rnd_gen (rd ());
    double r;
    int rootId = tree->getRootId();
    std::vector<int> nodeIds = tree->getNodesId();
    for(int i=0;i<nodeIds.size();i++){
        r = rng (rnd_gen);
        if(nodeIds.at(i) != rootId){
            tree->setDistanceToFather(nodeIds.at(i),r);
        }

    }

}

void CastorApplication::fixSmallBranchLength(){

    bpp:: ApplicationTools::displayMessage("fixing small branch lengths ");

    double r = MIN_BRANCH_LEN;
    int rootId = tree->getRootId();
    std::vector<int> nodeIds = tree->getNodesId();
    for(int i=0;i<nodeIds.size();i++){
        if(nodeIds.at(i) != rootId){
            tree->setDistanceToFather(nodeIds.at(i),r);
        }

    }

}

void CastorApplication::initTreeBranchLengthInput(std::map<std::string, std::string> &cmdArgs){

    bool midPointRootBrLengths = false;

    // Is the root has to be moved to the midpoint position along the branch that contains it ? If no, do nothing!
    midPointRootBrLengths = ApplicationTools::getBooleanParameter("midpoint_root_branch", cmdArgs, false,"", true, 2);

    if (midPointRootBrLengths){
        bpp::TreeTools::constrainedMidPointRooting(*this->tree);
    }
}

void CastorApplication::initTreeBranchLengthEqual(std::map<std::string, std::string> &cmdArgs){

    double value = 0.0;

    value = bpp::ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, 2);

    if (value <= 0){
        throw Exception("Value for branch length must be superior to 0");
    }

    bpp:: ApplicationTools::displayResult("Branch lengths set to", value);

    this->tree->setBranchLengths(value);

}

void CastorApplication::initTreeBranchLengthClock() {

    bpp::TreeTools::convertToClockTree(*this->tree, this->tree->getRootId(), true);

}

void CastorApplication::initTreeBranchLengthGrafen(std::map<std::string, std::string> &cmdArgs){

    std::string grafenHeight;
    double h = 0.0;
    double rho = 0.0;
    double nh = 0.0;

    grafenHeight = bpp::ApplicationTools::getStringParameter("height", cmdArgs, "input", "", true, 2);

    if (grafenHeight == "input") {
        h = bpp::TreeTools::getHeight(*this->tree, this->tree->getRootId());
    } else {
        h = bpp::TextTools::toDouble(grafenHeight);
        if (h <= 0){
            throw Exception("Height must be positive in Grafen's method.");
        }
    }

    bpp::ApplicationTools::displayResult("Total height", bpp::TextTools::toString(h));

    rho = bpp::ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, 2);

    bpp::ApplicationTools::displayResult("Grafen's rho", rho);

    bpp::TreeTools::computeBranchLengthsGrafen(*this->tree, rho);

    nh = bpp::TreeTools::getHeight(*this->tree, this->tree->getRootId());

    this->tree->scaleTree(h / nh);

}

void CastorApplication::initTreeBranchLength(){

    std::string cmdName;
    std::map<std::string, std::string> cmdArgs;

    // Setting branch lengths?
    this->initBrLenMethod = bpp::ApplicationTools::getStringParameter("init.brlen.method", this->getParams(), "Input","", true, 1);

    bpp::KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);

    if (cmdName == "Input") {

        this->initTreeBranchLengthInput(cmdArgs);

    } else if (cmdName == "Equal") {

        this->initTreeBranchLengthEqual(cmdArgs);

    } else if (cmdName == "Clock") {

        this->initTreeBranchLengthClock();

    } else if (cmdName == "Grafen") {

        this->initTreeBranchLengthGrafen(cmdArgs);

    } else if (cmdName == "exp") {

        this->initTreeBranchLengthExponentialDistribution(cmdArgs);

    } else{

        throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");

    }

    bpp::ApplicationTools::displayResult("Branch lengths", cmdName);

}

void CastorApplication::convertBppTree2Utree(){

    // Convert the bpp into utree for tree-search engine
    //auto utree = new Utree();
    this->utree = new UtreeBppUtils::Utree();

    UtreeBppUtils::convertTree_b2u(this->tree, this->utree, this->tm);

    if (this->PAR_alignment) {
        UtreeBppUtils::associateNode2Alignment(this->sequences, this->utree);
    } else {
        UtreeBppUtils::associateNode2Alignment(this->sites, this->utree);
    }

}

void CastorApplication::getTree(){

    this->tree = nullptr;
    this->initTreeOpt = bpp::ApplicationTools::getStringParameter("init.tree", this->getParams(), "user", "", false,1);

    bpp::ApplicationTools::displayResult("Initial tree", this->initTreeOpt);

    this->initTreeMethod();

    this->resolveMultifurcations();

    this->renameTreeNodes();

    // Try to write the current tree to file. This will be overwritten
    // by the optimized tree, but allow to check file existence before
    // running optimization!
    bpp::PhylogeneticsApplicationTools::writeTree(*this->tree, this->getParams());

    this->initTreeBranchLength();

    this->fixSmallBranchLength();

    DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(this->tree);

    this->convertBppTree2Utree();

    DLOG(INFO) << "Bidirectional map size: " << this->tm.size();
    DLOG(INFO) << "[Initial Utree Topology] " << this->utree->printTreeNewick(true);

    this->utree->addVirtualRootNode();

    // Once the tree has the root, then map it as well
    this->tm.insert(UtreeBppUtils::nodeassoc(this->tree->getRootId(), this->utree->rootnode->getVnode_id()));

}