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
#include <Bpp/Seq/Container/SiteContainerTools.h>
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
#include "inference_indel_rates.hpp"

#include "DistanceFactory.hpp"
#include "DistanceFactoryAngle.hpp"
#include "DistanceFactoryAlign.hpp"
#include "DistanceFactoryPrealigned.hpp"

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

void CastorApplication::getSubstitutionModel(bpp::Tree *tree){

    std::vector<std::string> keys;

    this->estimatePIPparameters = false;

    // Instantiate a substitution model and extend it with PIP
    if (this->PAR_model_indels) {

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
                        computeFrequenciesFromData = true;
                    }
                }

                this->baseModel += ",";
            }

            this->baseModel.pop_back();
            this->baseModel += ")";
            this->modelMap["model"] = this->baseModel;
        }

        // Instantiation of the canonical substitution model
        if (this->PAR_Alphabet.find("Codon") != std::string::npos) {
            this->smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabetNoGaps,this->gCode.get(), this->sites,this->modelMap, "", true,false, 0);
        }else if (this->PAR_Alphabet.find("Protein") != std::string::npos){
            this->smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabetNoGaps, this->gCode.get(), this->sites,this->modelMap, "", true, false, 0);
        } else {
            this->smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabet, this->gCode.get(), this->sites,this->modelMap, "", true, false, 0);
        }

        // If PIP, then check if lambda/mu initial values are estimated from the data
        if (this->modelMap.find("estimated") != this->modelMap.end()) {
            bpp::ApplicationTools::displayError("The use of the tag [observed] is obsolete. Use the tag [initFreqs] instead");
            exit(EXIT_FAILURE);
        } else if (this->modelMap.find("initFreqs") != this->modelMap.end()) {
            if (this->modelMap["initFreqs"] == "observed") {
                estimatePIPparameters = true;
            }
        } else if( (this->modelMap.find("lambda") == this->modelMap.end()) || (this->modelMap.find("mu") == this->modelMap.end())) {
            estimatePIPparameters = true;
        }

        if (estimatePIPparameters) {

            inference_indel_rates::infere_indel_rates_from_sequences(this->PAR_input_sequences,
                                                                     this->PAR_Alphabet,
                                                                     this->PAR_alignment,
                                                                     this->PAR_model_indels,
                                                                     this->getParams(),
                                                                     tree,
                                                                     this->lambda,
                                                                     this->mu,
                                                                     this->gCode.get(),
                                                                     this->modelMap);


            DLOG(INFO) << "[PIP model] Estimated PIP parameters from data using input sequences (lambda=" <<lambda << ",mu=" << mu << "," "I=" << lambda * mu << ")";

        } else {
            lambda = (this->modelMap.find("lambda") == this->modelMap.end()) ? 0.1 : std::stod(this->modelMap["lambda"]);
            mu = (this->modelMap.find("mu") == this->modelMap.end()) ? 0.2 : std::stod(this->modelMap["mu"]);
        }

        DLOG(INFO) << "[PIP model] Fixed PIP parameters to (lambda=" << lambda << ",mu=" << mu << "," "I="<< lambda * mu << ")";

        // Instantiate the corrisponding PIP model given the alphabet
        if (this->PAR_alignment) {
            if (this->PAR_Alphabet.find("DNA") != std::string::npos && this->PAR_Alphabet.find("Codon") == std::string::npos) {

                this->smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(this->alphabet), this->smodel, *this->sequences, lambda, mu,computeFrequenciesFromData);

            } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {

                this->smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(this->alphabet), this->smodel, *this->sequences, lambda, mu,computeFrequenciesFromData);

            } else if (this->PAR_Alphabet.find("Codon") != std::string::npos) {

                this->smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(this->alphabet), this->gCode.get(), this->smodel,*this->sequences, lambda, mu,computeFrequenciesFromData);

                bpp::ApplicationTools::displayWarning("Codon models are experimental in the current version... use with caution!");
                DLOG(WARNING) << "CODONS activated byt the program is not fully tested under these settings!";
            }
        } else {
            if (this->PAR_Alphabet.find("DNA") != std::string::npos && this->PAR_Alphabet.find("Codon") == std::string::npos) {

                this->smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(this->alphabet), this->smodel, *this->sites, lambda, mu,computeFrequenciesFromData);

            } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {

                this->smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(this->alphabet), this->smodel, *this->sites, lambda, mu,computeFrequenciesFromData);

            } else if (this->PAR_Alphabet.find("Codon") != std::string::npos) {

                this->smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(this->alphabet), this->gCode.get(), this->smodel,*this->sites, lambda, mu,computeFrequenciesFromData);

                bpp::ApplicationTools::displayWarning("Codon models are experimental in the current version... use with caution!");
                DLOG(WARNING) << "CODONS activated byt the program is not fully tested under these settings!";
            }
        }

    } else {
        bpp::SiteContainerTools::changeGapsToUnknownCharacters(*this->sites);
        this->smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabet, this->gCode.get(), this->sites,this->getParams(), "", true, false, 0);
    }

    DLOG(INFO) << "[Substitution model] Number of states: " << (int) this->smodel->getNumberOfStates();

    bpp::ApplicationTools::displayResult("Substitution model", this->smodel->getName());
    if (this->PAR_model_indels){
        bpp::ApplicationTools::displayResult("Indel parameter initial value",(estimatePIPparameters) ? "estimated" : "fixed");
    }

    ParameterList parameters = this->smodel->getParameters();
    for (size_t i = 0; i < parameters.size(); i++) {
        bpp::ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }

    for (size_t i = 0; i < this->smodel->getFrequencies().size(); i++) {
        bpp::ApplicationTools::displayResult("eq.freq(" + this->smodel->getAlphabet()->getName(i) + ")",TextTools::toString(this->smodel->getFrequencies()[i], 4));
    }

    bpp::StdStr s1;
    bpp::PhylogeneticsApplicationTools::printParameters(this->smodel, s1, 1, true);
    DLOG(INFO) << s1.str();

}

void CastorApplication::getTree(){

    this->tree = nullptr;
    this->initTreeOpt = ApplicationTools::getStringParameter("init.tree", this->getParams(), "user", "", false,1);

    ApplicationTools::displayResult("Initial tree", this->initTreeOpt);

    if (initTreeOpt == "user") {

        this->tree = PhylogeneticsApplicationTools::getTree(this->getParams());

        DLOG(INFO) << "[Input tree parser] Number of leaves" << this->tree->getNumberOfLeaves();

    } else if (initTreeOpt == "random") {

        vector<std::string> names = this->sites->getSequencesNames();
        this->tree = TreeTemplateTools::getRandomTree(names);
        this->tree->setBranchLengths(1.);

    } else if (this->initTreeOpt == "distance") {

        bpp::DistanceMatrix *distances;

        this->PAR_distance_method = ApplicationTools::getStringParameter("init.distance.method",this->getParams(), "nj");

        ApplicationTools::displayResult("Initial tree reconstruction method", this->PAR_distance_method);

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
            bpp::ApplicationTools::displayResult("Initial tree method", std::string("LZ compression"));

            try {
                this->PAR_distance_matrix = bpp::ApplicationTools::getAFilePath("init.distance.matrix.file",this->getParams(), true, true, "", false,"",0);
            } catch (bpp::Exception &e) {
                LOG(FATAL) << "Error when reading distance matrix file: " << e.message();
            }

            DLOG(INFO) << "initial tree method from LZ compression from matrix file" << PAR_distance_matrix;

            distances = InputUtils::parseDistanceMatrix(PAR_distance_matrix);

            bpp::BioNJ bionj(*distances, true, true, false);

            this->tree = bionj.getTree();

        } else if (token == "infere_distance_matrix") {

            int ALPHABET_DIM=0;
            int K=0;
            bool mldist_flag = true;
            bool mldist_gap_flag = false;
            double cutoff_dist = 1.0;
            double indel_rate = 1.0;

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

            DistanceFactoryPrographMSA::DistanceFactory *dist_factory_angle = new DistanceFactoryPrographMSA::DistanceFactoryAngle(ALPHABET_DIM, K);

            DistanceFactoryPrographMSA::DistanceMatrix dist_ml = dist_factory_angle->computePwDistances(this->sequences,
                                                                                                        ALPHABET_DIM,
                                                                                                        K,
                                                                                                        mldist_flag,
                                                                                                        mldist_gap_flag,
                                                                                                        cutoff_dist,
                                                                                                        indel_rate);

            bpp::DistanceMatrix *dist_ = new DistanceMatrix(this->sequences->getSequencesNames());

            for (int iii = 0; iii < this->sequences->getNumberOfSequences(); iii++) {
                for (int jjj = 0; jjj < this->sequences->getNumberOfSequences(); jjj++) {
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

            if (this->PAR_model_indels) {
                alphabetDistMethod = this->alphabet;
                parmap["model"] = this->modelMap["model"];
            } else {
                alphabetDistMethod = this->alphabetNoGaps;
                parmap["model"] = this->getParams()["model"];
            }



            allSites = SequenceApplicationTools::getSiteContainer(alphabetDistMethod, this->getParams());
            sitesDistMethod = SequenceApplicationTools::getSitesToAnalyse(*allSites, this->getParams());
            delete allSites;


            bpp::SubstitutionModel *smodel;
            if (this->PAR_model_indels) {

                // Instantiation of the canonical substitution model
                if (this->PAR_Alphabet.find("Codon") != std::string::npos ||
                        this->PAR_Alphabet.find("Protein") != std::string::npos) {
                    smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(this->alphabetNoGaps, this->gCode.get(),
                                                                                      sitesDistMethod, this->modelMap, "",
                                                                                      true,
                                                                                      false, 0);
                } else {
                    smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabetDistMethod,
                                                                                      this->gCode.get(), sitesDistMethod,
                                                                                      this->modelMap,
                                                                                      "", true,
                                                                                      false, 0);
                }


                if(this->modelMap.find("lambda") == this->modelMap.end() || this->modelMap.find("mu") == this->modelMap.end()){


                    //==================================================================================
                    //==================================================================================
                    //==================================================================================
                    //==================================================================================
                    //==================================================================================
                    // m@x:: new code
                    if(!this->tree){

                        double lambda_tmp = 10.0;
                        double mu_tmp = 0.1;

                        bpp::SubstitutionModel *smodel_tmp;
                        bpp::SubstitutionModel *smodel_copy = smodel->clone();

                        VectorSiteContainer *sitesDistMethod_tmp = sitesDistMethod->clone();
                        bpp::Alphabet *alphabetDistMethod_tmp = alphabetDistMethod->clone();

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

                        TransitionModel *dmodel_tmp;
                        // Get transition model from substitution model
                        if (!this->PAR_model_indels) {
                            dmodel_tmp = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod_tmp, this->gCode.get(),sitesDistMethod_tmp, parmap);
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
                            rDist_tmp = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());
                        }

                        // Remove gap characters since we are roughly estimating the initial topology
                        if (!this->PAR_model_indels) {
                            bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sitesDistMethod_tmp);
                        }

                        UnifiedDistanceEstimation distEstimation_tmp(dmodel_tmp, rDist_tmp, sitesDistMethod_tmp, 1, false);

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


                }else{
                    this->lambda = std::stod(this->modelMap["lambda"]);
                    this->mu = std::stod(this->modelMap["mu"]);
                }

                // Instatiate the corrisponding PIP model given the alphabet
                if (this->PAR_Alphabet.find("DNA") != std::string::npos &&
                        this->PAR_Alphabet.find("Codon") == std::string::npos) {
                    smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabetDistMethod), smodel,
                                         *sitesDistMethod, this->lambda, this->mu, false);
                } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
                    smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabetDistMethod), smodel,
                                        *sitesDistMethod, this->lambda, this->mu, false);
                } else if (this->PAR_Alphabet.find("Codon") != std::string::npos) {
                    smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabetDistMethod), this->gCode.get(),
                                           smodel, *sitesDistMethod,
                                           this->lambda,
                                           this->mu, false);
                    ApplicationTools::displayWarning(
                            "Codon models are experimental in the current version... use with caution!");
                    DLOG(WARNING) << "CODONS activated but the program is not fully tested under these settings!";
                }

            }

            TransitionModel *dmodel;
            // Get transition model from substitution model
            if (!this->PAR_model_indels) {
                dmodel = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod, this->gCode.get(),sitesDistMethod, parmap);
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
                rDist = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());
            }

            // Remove gap characters since we are roughly estimating the initial topology
            if (!this->PAR_model_indels) {
                bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sitesDistMethod);
            }

            UnifiedDistanceEstimation distEstimation(dmodel, rDist, sitesDistMethod, 1, false);

            if (this->PAR_distance_method.find("-ml") != std::string::npos) {

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
                ParameterList allParameters = dmodel->getParameters();
                allParameters.addParameters(rDist->getParameters());

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
            /**/

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
    this->tree->setNodeName(tree->getRootId(), "root");
    UtreeBppUtils::renameInternalNodes(tree);

    // Try to write the current tree to file. This will be overwritten by the optimized tree,
    // but allow to check file existence before running optimization!
    PhylogeneticsApplicationTools::writeTree(*this->tree, this->getParams());

    // Setting branch lengths?
    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", this->getParams(), "Input",
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
    //auto utree = new Utree();
    this->utree = new UtreeBppUtils::Utree();

    UtreeBppUtils::convertTree_b2u(this->tree, utree, this->tm);
    if (this->PAR_alignment) {
        UtreeBppUtils::associateNode2Alignment(this->sequences, utree);
    } else {
        UtreeBppUtils::associateNode2Alignment(this->sites, utree);
    }

    DLOG(INFO) << "Bidirectional map size: " << this->tm.size();
    DLOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

    utree->addVirtualRootNode();
    // Once the tree has the root, then map it as well
    //tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode));

    this->tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode->getVnode_id()));



}