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

void CastorApplication::getCLIarguments(std::map<std::string, std::string> &modelMap){

    this->PAR_model_substitution = ApplicationTools::getStringParameter("model", this->getParams(), "JC69","", true, true);

    this->PAR_alignment = ApplicationTools::getBooleanParameter("alignment", this->getParams(), false);

    // Split model string description and test if PIP is required
    std::string modelStringName;
    KeyvalTools::parseProcedure(PAR_model_substitution, modelStringName, modelMap);

    this->PAR_model_indels = modelStringName == "PIP";

}

void CastorApplication::getAlphabet(bpp::Alphabet *alphabetNoGaps,
                 std::unique_ptr<GeneticCode> &gCode,
                 bpp::Alphabet *alphabet){

    this->PAR_Alphabet = ApplicationTools::getStringParameter("alphabet", this->getParams(), "DNA", "",true, true);

    // Alphabet without gaps
    alphabetNoGaps = bpp::SequenceApplicationTools::getAlphabet(this->getParams(), "", false, false);

    // Codon alphabet ?
    this->codonAlphabet = false;

    // Alphabet used for all the computational steps (it can allows for gap extension)
    if (this->PAR_model_indels) {

        if (this->PAR_Alphabet.find("DNA") != std::string::npos) {
            alphabet = new bpp::DNA_EXTENDED();
        } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
            alphabet = new bpp::ProteicAlphabet_Extended();
        } else if (this->PAR_Alphabet.find("Codon") != std::string::npos){
            alphabet = new bpp::DNA_EXTENDED();
            alphabet = new CodonAlphabet_Extended(dynamic_cast<bpp::NucleicAlphabet *>(alphabet));
            this->codonAlphabet = true;
        }

    } else {

        if (this->PAR_Alphabet.find("DNA") != std::string::npos ) {
            alphabet = new bpp::DNA();
        } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
            alphabet = new bpp::ProteicAlphabet();
        } else if (this->PAR_Alphabet.find("Codon") != std::string::npos){
            alphabet = new bpp::DNA();
            alphabet = new CodonAlphabet(dynamic_cast<bpp::NucleicAlphabet *>(alphabet));
            this->codonAlphabet = true;
        }
    }

    // Alphabet used for codon models
    if (this->codonAlphabet) {
        std::string codeDesc = ApplicationTools::getStringParameter("genetic_code", this->getParams(), "Standard","", true, true);

        ApplicationTools::displayResult("Genetic Code", codeDesc);

        gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(dynamic_cast<bpp::CodonAlphabet *>(alphabetNoGaps)->getNucleicAlphabet(),codeDesc));

    }
}

void CastorApplication::getData(bpp::SequenceContainer *sequences,
             bpp::SiteContainer *sites,
             bpp::Alphabet *alphabet){

    this->PAR_input_sequences = ApplicationTools::getAFilePath("input.sequence.file", this->getParams(),
                                                               true, true, "", false, "", 1);

    try {

        ApplicationTools::displayBooleanResult("Aligned sequences", !this->PAR_alignment);

        if (this->PAR_alignment) {

            // If the user requires the computation of an alignment, then the input file is made of unaligned sequences
            bpp::Fasta seqReader;
            sequences = seqReader.readSequences(PAR_input_sequences, alphabet);
            ApplicationTools::displayResult("Number of sequences",TextTools::toString(sequences->getNumberOfSequences()));

        } else {

            VectorSiteContainer *allSites = SequenceApplicationTools::getSiteContainer(alphabet,this->getParams());

            sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, this->getParams(), "", true,!this->PAR_model_indels, true, 1);

            AlignmentUtils::checkAlignmentConsistency(*sites);
            ApplicationTools::displayResult("Number of sequences",TextTools::toString(sites->getNumberOfSequences()));
            ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

            delete allSites;

        }


    } catch (bpp::Exception &e) {
        LOG(FATAL) << "Error when reading sequence file due to: " << e.message();
    }

}

void CastorApplication::getTree(bpp::Alphabet *alphabet,bpp::Alphabet *alphabetNoGaps,bpp::Tree *tree,
             bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,std::map<std::string, std::string> &modelMap,
             std::unique_ptr<GeneticCode> &gCode,UtreeBppUtils::treemap &tm,tshlib::Utree *utree){

    this->PAR_initTreeOpt = ApplicationTools::getStringParameter("init.tree", this->getParams(), "user", "", false,1);

    ApplicationTools::displayResult("Initial tree", this->PAR_initTreeOpt);

    this->getInitTree(tree,sites,sequences,alphabet);

    // If the tree has multifurcation, then resolve it with midpoint rooting
    this->resolveMultifurcation(tree);

    // Rename internal nodes with standard Vxx * where xx is a progressive number
    tree->setNodeName(tree->getRootId(), "root");
    UtreeBppUtils::renameInternalNodes(tree);

    // Try to write the current tree to file. This will be overwritten by the optimized tree,
    // but allow to check file existence before running optimization!
    PhylogeneticsApplicationTools::writeTree(*tree, this->getParams());

    // Setting branch lengths?
    this->initBranchLength(tree);

    DLOG(INFO) << "[Initial Tree Topology] " << OutputUtils::TreeTools::writeTree2String(tree);

    // Convert the bpp into utree for tree-search engine
    utree = new UtreeBppUtils::Utree();

    UtreeBppUtils::convertTree_b2u(tree, utree, tm);
    if (this->PAR_alignment) {
        UtreeBppUtils::associateNode2Alignment(sequences, utree);
    } else {
        UtreeBppUtils::associateNode2Alignment(sites, utree);
    }

    DLOG(INFO) << "Bidirectional map size: " << tm.size();
    DLOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

    utree->addVirtualRootNode();
    // Once the tree has the root, then map it as well

    tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode->getVnode_id()));

}

void CastorApplication::getInitTree(bpp::Tree *tree,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet){

    if (this->PAR_initTreeOpt == "user") {

        this->getUserTree(tree);

    } else if (this->PAR_initTreeOpt == "random") {

        this->getRandomTree(tree,sites);

    } else if (this->PAR_initTreeOpt == "distance") {

        this->getDistanceTree(tree,sites,sequences,alphabet);

    } else{

        throw Exception("Unknown init tree method.");

    }

}

void CastorApplication::getUserTree(bpp::Tree *tree){

    tree = PhylogeneticsApplicationTools::getTree(this->getParams());

    DLOG(INFO) << "[Input tree parser] Number of leaves" << tree->getNumberOfLeaves();

}

void CastorApplication::getRandomTree(bpp::Tree *tree,bpp::SiteContainer *sites){

    vector<string> names = sites->getSequencesNames();

    tree = TreeTemplateTools::getRandomTree(names);

    tree->setBranchLengths(1.);

}

void CastorApplication::getDistanceTree(bpp::Tree *tree,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet){

    this->PAR_distance_method = ApplicationTools::getStringParameter("init.distance.method",
                                                                     this->getParams(), "nj");

    ApplicationTools::displayResult("Initial tree reconstruction method", PAR_distance_method);

    AgglomerativeDistanceMethod *distMethod = nullptr;
    std::string token = PAR_distance_method.substr(0, this->PAR_distance_method.find("-"));
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

        this->getDistmatrixTree(tree);

    } else if (token == "infere_distance_matrix") {

        this->getInfere_distance_matrix_Tree(tree,sequences,alphabet);

    } else{

        throw Exception("Unknown tree reconstruction method.");

    }

    //check_this_code(); //????????????????????????

}

void CastorApplication::getDistmatrixTree(bpp::Tree *tree){

    // Use a distance matrix provided by the user
    ApplicationTools::displayResult("Initial tree method", std::string("LZ compression"));

    bpp::DistanceMatrix *distances;



    try {
        this->PAR_distance_matrix = ApplicationTools::getAFilePath("init.distance.matrix.file",
                                                             this->getParams(), true, true, "", false,
                                                             "",
                                                             0);
    } catch (bpp::Exception &e) {

        LOG(FATAL) << "Error when reading distance matrix file: " << e.message();

    }

    DLOG(INFO) << "initial tree method from LZ compression from matrix file" << PAR_distance_matrix;

    distances = InputUtils::parseDistanceMatrix(PAR_distance_matrix);

    bpp::BioNJ bionj(*distances, true, true, false);

    tree = bionj.getTree();

}

void CastorApplication::getInfere_distance_matrix_Tree(bpp::Tree *tree,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet){

    int ALPHABET_DIM;
    int K;
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

    if (!sequences) {
        bpp::Fasta seqReader;
        sequences = seqReader.readSequences(this->PAR_input_sequences, alphabet);
    }

    DistanceFactoryPrographMSA::DistanceFactory *dist_factory_angle = new DistanceFactoryPrographMSA::DistanceFactoryAngle(ALPHABET_DIM, K);

    DistanceFactoryPrographMSA::DistanceMatrix dist_ml = dist_factory_angle->computePwDistances(sequences,
                                                                                                ALPHABET_DIM,
                                                                                                K,
                                                                                                mldist_flag,
                                                                                                mldist_gap_flag,
                                                                                                cutoff_dist,
                                                                                                indel_rate);

    bpp::DistanceMatrix *dist_ = new DistanceMatrix(sequences->getSequencesNames());

    for (int iii = 0; iii < sequences->getNumberOfSequences(); iii++) {
        for (int jjj = 0; jjj < sequences->getNumberOfSequences(); jjj++) {
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

}

void CastorApplication::resolveMultifurcation(bpp::Tree *tree){

    auto ttree_ = new TreeTemplate<Node>(*tree);

    if (ttree_->getRootNode()->getNumberOfSons() > 2) {
        TreeTemplateTools::midRoot(*(ttree_), TreeTemplateTools::MIDROOT_VARIANCE, false);
        tree = ttree_;
    }

}

void CastorApplication::initBranchLength(bpp::Tree *tree){

    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", this->getParams(), "Input",
                                                                  "", true, 1);
    string cmdName;
    map<string, string> cmdArgs;
    KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);

    if (cmdName == "Input") {

        // Is the root has to be moved to the midpoint position along the branch that contains it ? If no, do nothing!
        bool midPointRootBrLengths = ApplicationTools::getBooleanParameter("midpoint_root_branch", cmdArgs, false,"", true, 2);

        if (midPointRootBrLengths){
            TreeTools::constrainedMidPointRooting(*tree);
        }

    } else if (cmdName == "Equal") {

        double value = ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, 2);

        if (value <= 0){
            throw Exception("Value for branch length must be superior to 0");
        }

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

            if (h <= 0){
                throw Exception("Height must be positive in Grafen's method.");
            }

        }

        ApplicationTools::displayResult("Total height", TextTools::toString(h));

        double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, 2);

        ApplicationTools::displayResult("Grafen's rho", rho);

        TreeTools::computeBranchLengthsGrafen(*tree, rho);

        double nh = TreeTools::getHeight(*tree, tree->getRootId());

        tree->scaleTree(h / nh);

    } else{
        throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");
    }

    ApplicationTools::displayResult("Branch lengths", cmdName);

}

void CastorApplication::getSubstitutionModel(std::map<std::string, std::string> &modelMap,bpp::SubstitutionModel *smodel,
                          bpp::TransitionModel *model,std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                          bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                          bpp::Tree *tree){

    // Instantiate a substitution model and extend it with PIP
    if (this->PAR_model_indels) {

        this->getSubstitutionModelIndel(modelMap,smodel,model,gCode,alphabet,alphabetNoGaps,sites,sequences,tree);

    } else {

        bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites,this->getParams(), "", true, false, 0);

    }

    DLOG(INFO) << "[Substitution model] Number of states: " << (int) smodel->getNumberOfStates();

    ApplicationTools::displayResult("Substitution model", smodel->getName());

    if (this->PAR_model_indels){
        ApplicationTools::displayResult("Indel parameter initial value",(this->PAR_estimatePIPparameters) ? "estimated" : "fixed");
    }

    bpp::StdStr stream;
    bpp::PhylogeneticsApplicationTools::printParameters(smodel, stream, 1, true);
    DLOG(INFO) << stream.str();

}

void CastorApplication::getModelMap(std::map<std::string, std::string> &modelMap,bpp::SubstitutionModel *smodel){

    this->PAR_computeFrequenciesFromData = false;

    // If frequencies are estimated from the data, but there is no alignment, then flag it.
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
                    this->PAR_computeFrequenciesFromData = true;
                }
            }
            baseModel += ",";
        }
        baseModel.pop_back();
        baseModel += ")";
        modelMap["model"] = baseModel;
    }

}

void CastorApplication::getIndelRates(std::map<std::string, std::string> &modelMap,bpp::Tree *tree,std::unique_ptr<GeneticCode> &gCode){

    this->PAR_estimatePIPparameters = false;

    if (modelMap.find("estimated") != modelMap.end()) {

        ApplicationTools::displayError("The use of the tag [observed] is obsolete. Use the tag [initFreqs] instead");
        exit(EXIT_FAILURE);

    } else if (modelMap.find("initFreqs") != modelMap.end()) {

        if (modelMap["initFreqs"] == "observed") {
            this->PAR_estimatePIPparameters = true;
        }

    } else if( (modelMap.find("lambda") == modelMap.end()) || (modelMap.find("mu") == modelMap.end())) {

        this->PAR_estimatePIPparameters = true;

    }

    if (this->PAR_estimatePIPparameters) {

        inference_indel_rates::infere_indel_rates_from_sequences(this->PAR_input_sequences,this->PAR_Alphabet,tree,
                                                                 this->lambda,this->mu,gCode.get(),modelMap);

        DLOG(INFO) << "[PIP model] Estimated PIP parameters from data using input sequences (lambda=" <<
                   this->lambda << ",mu=" << this->mu << "," "I=" << this->lambda * this->mu << ")";

    } else {
        this->lambda = (modelMap.find("lambda") == modelMap.end()) ? 1.0 : std::stod(modelMap["lambda"]);
        this->mu = (modelMap.find("mu") == modelMap.end()) ? 0.1 : std::stod(modelMap["mu"]);
    }

    DLOG(INFO) << "[PIP model] Fixed PIP parameters to (lambda=" << this->lambda << ",mu=" << this->mu << "," "I="<< this->lambda * this->mu << ")";

}

void CastorApplication::getSubstitutionModelIndel(std::map<std::string, std::string> &modelMap,bpp::SubstitutionModel *smodel,
                                             bpp::TransitionModel *model,std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                                             bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                                             bpp::Tree *tree){


    this->getModelMap(modelMap,smodel);

    // Instantiation of the canonical substitution model
    this->getIndelRates(modelMap,tree,gCode);

    const SequenceContainer *data;
    if (this->PAR_alignment) {
        data = sequences;
    }else{
        data = sites;
    }

    // Instantiate the corrisponding PIP model given the alphabet
    if (this->PAR_Alphabet.find("DNA") != std::string::npos && this->PAR_Alphabet.find("Codon") == std::string::npos) {

        smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabet), smodel, *data, this->lambda, this->mu,
                             this->PAR_computeFrequenciesFromData);

    } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {

        smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabet), smodel, *data, this->lambda, this->mu,
                            this->PAR_computeFrequenciesFromData);

    } else if (this->PAR_Alphabet.find("Codon") != std::string::npos) {

        smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabet), gCode.get(), smodel,
                               *data, this->lambda, this->mu,
                               this->PAR_computeFrequenciesFromData);

        ApplicationTools::displayWarning("Codon models are experimental in the current version... use with caution!");

        DLOG(WARNING) << "CODONS activated byt the program is not fully tested under these settings!";

    }

}

void CastorApplication::getParameters(ParameterList &parameters,bpp::SubstitutionModel *smodel){

    parameters = smodel->getParameters();

    // print parameters
    for (size_t i = 0; i < parameters.size(); i++) {
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }

    for (size_t i = 0; i < smodel->getFrequencies().size(); i++) {
        ApplicationTools::displayResult("eq.freq(" + smodel->getAlphabet()->getName(i) + ")",TextTools::toString(smodel->getFrequencies()[i], 4));
    }

}

void CastorApplication::getASVR(bpp::DiscreteDistribution *rDist,
             bpp::SubstitutionModel *smodel){

    // Among site rate variation (ASVR)
    if (smodel->getNumberOfStates() >= 2 * smodel->getAlphabet()->getSize()) {

        // Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();

    } else {

        rDist = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());

    }

}

void CastorApplication::getMSA(progressivePIP *proPIP,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
            bpp::DiscreteDistribution *rDist,bpp::SubstitutionModel *smodel,UtreeBppUtils::treemap &tm,
            bpp::Tree *tree,UtreeBppUtils::Utree *utree){


    this->PAR_output_file_msa = ApplicationTools::getAFilePath("output.msa.file", this->getParams(), false,false, "", true, "", 1);

    this->PAR_alignment_version = ApplicationTools::getStringParameter("alignment.version", this->getParams(), "cpu", "",true, 0);

    this->PAR_alignment_sbsolutions = ApplicationTools::getIntParameter("alignment.sb_solutions", this->getParams(), 1, "", true, 0);

    this->PAR_alignment_sbtemperature = ApplicationTools::getDoubleParameter("alignment.sb_temperature", this->getParams(), 1.0, "", true, 0);

    ApplicationTools::displayMessage("\n[Computing the multi-sequence alignment]");
    ApplicationTools::displayResult("Aligner optimised for:", this->PAR_alignment_version);

    ApplicationTools::displayBooleanResult("Stochastic backtracking active", this->PAR_alignment_sbsolutions > 1);
    if (this->PAR_alignment_sbsolutions > 1) {
        ApplicationTools::displayResult("Number of stochastic solutions:",TextTools::toString(this->PAR_alignment_sbsolutions));
    }

    DLOG(INFO) << "[Alignment sequences] Starting MSA inference using Pro-PIP...";

    // Execute alignment on post-order node list
    std::vector<tshlib::VirtualNode *> ftn = utree->getPostOrderNodeList();

    enumDP3Dversion DPversion = CPU; // DP3D version

    if (PAR_alignment_version.find("cpu") != std::string::npos) {
        DPversion = CPU; // slower but uses less memory
        this->PAR_alignment_sbsolutions = 1;
    } else if (PAR_alignment_version.find("ram") != std::string::npos) {
        DPversion = RAM; // faster but uses more memory
        this->PAR_alignment_sbsolutions = 1;
    } else if (PAR_alignment_version.find("sb") != std::string::npos) {
        DPversion = SB;  // stochastic backtracking version
    } else {
        ApplicationTools::displayError("The user specified an unknown alignment.version. The execution will not continue.");
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    proPIP = new bpp::progressivePIP(utree,               // tshlib tree
                                     tree,                // bpp tree
                                     smodel,              // substitution model
                                     tm,                  // tree-map
                                     sequences,           // un-aligned input sequences
                                     rDist,               // rate-variation among site distribution
                                     this->getSeed());  // seed for random number generation

    proPIP->_initializePIP(ftn,          // list of tshlib nodes in on post-order (correct order of execution)
                           DPversion,    // version of the alignment algorithm
                           this->PAR_alignment_sbsolutions,       // number of suboptimal MSAs
                           this->PAR_alignment_sbtemperature); // to tune the greedyness of the sub-optimal solution

    proPIP->PIPnodeAlign(); // align input sequences with a DP algorithm under PIP

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    ApplicationTools::displayResult("\nAlignment elapsed time (msec):",TextTools::toString((double) duration, 4));


    // convert PIPmsa into a sites objects
    sites = PIPmsaUtils::PIPmsa2Sites(proPIP->alphabet_,
                                      *(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getseqNames()),
                                      *(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getMSA()) );

}

void CastorApplication::getOptParams(bpp::AbstractHomogeneousTreeLikelihood *tl){

    tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(Optimizators::optimizeParameters(tl,
                                                                                            tl->getParameters(),
                                                                                            this->getParams(),
                                                                                            "",
                                                                                            true,
                                                                                            true,
                                                                                            0));

    this->logL = tl->getLogLikelihood();

}

void CastorApplication::initLK(bpp::TransitionModel *model,bpp::SubstitutionModel *smodel,bpp::AbstractHomogeneousTreeLikelihood *tl,
            UtreeBppUtils::treemap &tm,bpp::DiscreteDistribution *rDist,std::unique_ptr<GeneticCode> &gCode,
            bpp::Alphabet *alphabet,bpp::SiteContainer *sites,bpp::Tree *tree,UtreeBppUtils::Utree *utree){

    // Get transition model from substitution model
    if (!this->PAR_model_indels) {
        model = bpp::PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), sites,this->getParams(), "", true, false, 0);
    } else {
        unique_ptr<TransitionModel> test;
        test.reset(smodel);
        model = test.release();
    }

    // Initialise likelihood functions
    if (!this->PAR_model_indels) {

        tl = new bpp::UnifiedTSHomogeneousTreeLikelihood(*tree, *sites, model, rDist, utree, &tm, true,this->getParams(), "", false, false,false);

    } else {

        tl = new bpp::UnifiedTSHomogeneousTreeLikelihood_PIP(*tree, *sites, model, rDist, utree, &tm, true,this->getParams(), "", false,false, false);
    }

    ApplicationTools::displayResult("Tree likelihood model", std::string("Homogeneous"));

    tl->initialize();

}

void CastorApplication::getParSanityCheck(bpp::AbstractHomogeneousTreeLikelihood *tl,
                       bpp::SiteContainer *sites,
                       std::unique_ptr<GeneticCode> &gCode){

    //Listing parameters
    string paramNameFile = ApplicationTools::getAFilePath("output.parameter_names.file", this->getParams(), false,
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
        if (this->codonAlphabet) {
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
                exit(EXIT_FAILURE);
        }
        bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites",
                                                                     this->getParams(), false, "",
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
            exit(EXIT_FAILURE);
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

}

void CastorApplication::output(bpp::Tree *tree,
            UtreeBppUtils::Utree *utree,
            bpp::SiteContainer *sites,
            bpp::AbstractHomogeneousTreeLikelihood *tl,
            UtreeBppUtils::treemap &tm,
            ParameterList &parameters,
            bpp::DiscreteDistribution *rDist){


    // Export alignment to file
    if (this->PAR_output_file_msa.find("none") == std::string::npos) {
        DLOG(INFO) << "[Alignment sequences]\t The final alignment can be found in " << PAR_output_file_msa;
        bpp::Fasta seqWriter;
        seqWriter.writeAlignment(TextUtils::appendToFilePath(this->PAR_output_file_msa, "initial"), *sites, true);
    }

    // Export final tree (if nexus is required, then our re-implementation of the the nexus writer is called)
    tree = new TreeTemplate<Node>(tl->getTree());

    this->PAR_output_tree_format = ApplicationTools::getStringParameter("output.tree.format",this->getParams(),"Newick","",true,true);

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
        ApplicationTools::displayResult("Output annotation to file", PAR_output_annotation_file);
        OutputUtils::writeTreeAnnotations2TSV(tree, PAR_output_annotation_file);

    }

    // Write parameters to screen:
    ApplicationTools::displayResult("Final Log likelihood", TextTools::toString(logL, 15));

    parameters = tl->getSubstitutionModelParameters();
    for (size_t i = 0; i < parameters.size(); i++) {
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }

    parameters = tl->getRateDistributionParameters();
    for (size_t i = 0; i < parameters.size(); i++) {
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }

    // Checking convergence:
    PhylogeneticsApplicationTools::checkEstimatedParameters(tl->getParameters());

    // Write parameters to file (according to arguments)
    OutputUtils::exportOutput(tl, sites, this->getParams());

    // Compute support measures
    this->PAR_support = ApplicationTools::getStringParameter("support", this->getParams(), "", "", true,
                                                             true);
    if (this->PAR_support == "bootstrap") {
        ApplicationTools::displayMessage("\n[Tree support measures]");

        bpp::Bootstrap(tl, *sites, rDist, utree, &tm, this->getParams(), "support.");
    }

}

/*
void CastorApplication::check_this_code(bpp::Alphabet *alphabet,bpp::Alphabet *alphabetNoGaps,bpp::Tree *tree,
                                        bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,std::map<std::string, std::string> &modelMap,
                                        std::unique_ptr<GeneticCode> &gCode,UtreeBppUtils::treemap &tm,tshlib::Utree *utree,
                                        bpp::DistanceMethod *distMethod){

    if (false) {

        // Compute bioNJ tree using the GTR model
        map<std::string, std::string> parmap;

        VectorSiteContainer *allSites;
        VectorSiteContainer *sitesDistMethod;
        bpp::Alphabet *alphabetDistMethod;

        if (this->PAR_model_indels) {
            alphabetDistMethod = alphabet;
            parmap["model"] = modelMap["model"];
        } else {
            alphabetDistMethod = alphabetNoGaps;
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


            if(modelMap.find("lambda") == modelMap.end() || modelMap.find("mu") == modelMap.end()){

                inference_indel_rates::infere_indel_rates_from_sequences(this->PAR_input_sequences,
                                                                         this->PAR_Alphabet,
                                                                         tree,
                                                                         this->lambda,
                                                                         this->mu,
                                                                         gCode.get(),
                                                                         modelMap);

            }else{
                lambda = std::stod(modelMap["lambda"]);
                mu = std::stod(modelMap["mu"]);
            }

            // Instatiate the corrisponding PIP model given the alphabet
            if (this->PAR_Alphabet.find("DNA") != std::string::npos &&
                this->PAR_Alphabet.find("Codon") == std::string::npos) {
                smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabetDistMethod), smodel,
                                     *sitesDistMethod, lambda, mu, false);
            } else if (this->PAR_Alphabet.find("Protein") != std::string::npos) {
                smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabetDistMethod), smodel,
                                    *sitesDistMethod, lambda, mu, false);
            } else if (this->PAR_Alphabet.find("Codon") != std::string::npos) {
                smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabetDistMethod), gCode.get(),
                                       smodel, *sitesDistMethod,
                                       lambda,
                                       mu, false);
                ApplicationTools::displayWarning(
                        "Codon models are experimental in the current version... use with caution!");
                DLOG(WARNING) << "CODONS activated but the program is not fully tested under these settings!";
            }

        }

        //Initialize model to compute the distance tree
        //TransitionModel *dmodel = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod, gCode.get(), sitesDistMethod, parmap);
        TransitionModel *dmodel;
        // Get transition model from substitution model
        if (!this->PAR_model_indels) {
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
            rDist = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());
        }

        // Remove gap characters since we are roughly estimating the initial topology
        if (!this->PAR_model_indels) {
            bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sitesDistMethod);
        }

        UnifiedDistanceEstimation distEstimation(dmodel, rDist, sitesDistMethod, 1, false);

        if (PAR_distance_method.find("-ml") != std::string::npos) {

            std::string PAR_optim_distance = ApplicationTools::getStringParameter(
                    "init.distance.optimization.method", this->getParams(),
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
                                                                           this->getParams(), 2);
            string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", this->getParams(),
                                                           false, false);
            auto *messenger = (mhPath == "none") ? nullptr : (mhPath == "std") ? ApplicationTools::message.get()
                                                                               : new StlOutputStream(
                            new ofstream(mhPath.c_str(), ios::out));
            ApplicationTools::displayResult("Initial tree optimization handler", mhPath);

            // Optimisation method profiler
            string prPath = ApplicationTools::getAFilePath("optimization.profiler", this->getParams(), false,
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

        //delete sitesDistMethod;
        //delete distMethod;

    } //m@x
    //----------------------------------------------------------------

}

*/
