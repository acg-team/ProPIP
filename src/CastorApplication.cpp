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

    // Init values
    PAR_model_substitution_ = "";
    PAR_alphabet_ = "";
    PAR_input_sequences_ = "";
    PAR_distance_method_ = "";
    PAR_distance_matrix_ = "";
    PAR_optim_distance_ = "";
    PAR_output_file_msa_ = "";
    PAR_alignment_version_ = "";
    PAR_output_tree_format_ = "";
    PAR_output_annotation_file_ = "";
    PAR_support_ = "";
    PAR_init_tree_opt_ = "";

    PAR_alignment_sb_solutions_ = 1;

    PAR_alignment_sb_temperature_ = 0.0;

    PAR_estimate_pip_parameters_ = false;
    PAR_compute_frequencies_from_data_ = false;
    PAR_alignment_ = false;
    PAR_model_indels_ = false;

    PAR_codon_alphabet_ = false;

    lambda_ = 1.0;
    mu_ = 0.1;
    logL_ = 0.0;

    smodel = nullptr;

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

void CastorApplication::checkInputArgument(int argc) {

    if (argc < 2) {
        this->help();
        exit(EXIT_SUCCESS);
    } else {
        this->banner();
        this->startTimer();
    }

}

void CastorApplication::getParseProcedure(std::map<std::string, std::string> &modelMap,std::string &modelStringName){

    // Split model string description and test if PIP is required
    KeyvalTools::parseProcedure(this->PAR_model_substitution_, modelStringName, modelMap);

}

bool CastorApplication::isPIP(std::string modelStringName){

    return (modelStringName == "PIP");
}

std::map<std::string, std::string> CastorApplication::getCLIarguments(){

    std::map<std::string, std::string> modelMap;

    std::string modelStringName;

    this->PAR_model_substitution_ = ApplicationTools::getStringParameter("model", this->getParams(), "JC69","", true, true);

    this->getParseProcedure(modelMap,modelStringName);

    this->PAR_init_tree_opt_ = ApplicationTools::getStringParameter("init.tree", this->getParams(), "user", "", false,1);

    this->PAR_distance_method_ = ApplicationTools::getStringParameter("init.distance.method",this->getParams(), "nj");

    this->PAR_model_indels_ = isPIP(modelStringName);

    this->PAR_alignment_ = ApplicationTools::getBooleanParameter("alignment", this->getParams(), false);

    this->PAR_alphabet_ = ApplicationTools::getStringParameter("alphabet", this->getParams(), "DNA", "",true, true);

    this->PAR_input_sequences_ = ApplicationTools::getAFilePath("input.sequence.file", this->getParams(),true, true, "", false, "", 1);

    ///////////////////////////////////////////
    // tree init method
    this->PAR_optim_distance_ = ApplicationTools::getStringParameter("init.distance.optimization.method", this->getParams(),"init");

    ApplicationTools::displayResult("Initial tree model parameters estimation method",this->PAR_optim_distance_);

    if (this->PAR_optim_distance_ != "init" &&
        this->PAR_optim_distance_ != "pairwise" &&
        this->PAR_optim_distance_ != "iterations"){

        throw Exception("Unknown parameter estimation procedure '" + this->PAR_optim_distance_ + "'.");
    }
    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // Alignment options
    this->PAR_output_file_msa_ = ApplicationTools::getAFilePath("output.msa.file", this->getParams(), false,false, "", true, "", 1);

    this->PAR_alignment_version_ = ApplicationTools::getStringParameter("alignment.version", this->getParams(), "cpu", "",true, 0);

    this->PAR_alignment_sb_solutions_ = ApplicationTools::getIntParameter("alignment.sb_solutions", this->getParams(), 1, "", true, 0);

    this->PAR_alignment_sb_temperature_ = ApplicationTools::getDoubleParameter("alignment.sb_temperature", this->getParams(), 1.0, "", true, 0);
    ///////////////////////////////////////////

    return modelMap;
}

bpp::Alphabet* CastorApplication::getAlphabetIndel(){

    bpp::Alphabet *alphabet = nullptr;

    if (this->PAR_alphabet_.find("DNA") != std::string::npos) {
        alphabet = new bpp::DNA_EXTENDED();
        this->PAR_codon_alphabet_ = false;
    } else if (this->PAR_alphabet_.find("Protein") != std::string::npos) {
        alphabet = new bpp::ProteicAlphabet_Extended();
        this->PAR_codon_alphabet_ = false;
    } else if (this->PAR_alphabet_.find("Codon") != std::string::npos){
        alphabet = new CodonAlphabet_Extended(dynamic_cast<bpp::NucleicAlphabet *>(new bpp::DNA_EXTENDED()));
        this->PAR_codon_alphabet_ = true;
    }

    return alphabet;
}

bpp::Alphabet* CastorApplication::getAlphabetNoIndel(){

    bpp::Alphabet *alphabet = nullptr;

    if (this->PAR_alphabet_.find("DNA") != std::string::npos ) {
        alphabet = new bpp::DNA();
        this->PAR_codon_alphabet_ = false;
    } else if (this->PAR_alphabet_.find("Protein") != std::string::npos) {
        alphabet = new bpp::ProteicAlphabet();
        this->PAR_codon_alphabet_ = false;
    } else if (this->PAR_alphabet_.find("Codon") != std::string::npos){
        alphabet = new CodonAlphabet(dynamic_cast<bpp::NucleicAlphabet *>(new bpp::DNA()));
        this->PAR_codon_alphabet_ = true;
    }

    return alphabet;
}

bpp::Alphabet* CastorApplication::getAlphabetNoGaps(){

    bpp::Alphabet *alphabetNoGaps;

    // Alphabet without gaps
    alphabetNoGaps = bpp::SequenceApplicationTools::getAlphabet(this->getParams(), "", false, false);

    return alphabetNoGaps;
}

bpp::Alphabet* CastorApplication::getAlphabet(){

    bpp::Alphabet *alphabet;

    // Alphabet used for all the computational steps (it can allows for gap extension)
    if (this->PAR_model_indels_) {
        alphabet = this->getAlphabetIndel();
    } else {
        alphabet = this->getAlphabetNoIndel();
    }

    return alphabet;
}

std::unique_ptr<GeneticCode> CastorApplication::getGcode(bpp::Alphabet *alphabetNoGaps){

    std::unique_ptr<GeneticCode> gCode;

    if (this->PAR_codon_alphabet_) {

        std::string codeDesc = ApplicationTools::getStringParameter("genetic_code", this->getParams(), "Standard", "",true, true);

        ApplicationTools::displayResult("Genetic Code", codeDesc);

        gCode.reset(bpp::SequenceApplicationTools::getGeneticCode(dynamic_cast<bpp::CodonAlphabet *>(alphabetNoGaps)->getNucleicAlphabet(), codeDesc));

    }else{

        return NULL;

    }

    return gCode;
}

bpp::SiteContainer* CastorApplication::getAlignedSequences(bpp::Alphabet *alphabet){

    bpp::SiteContainer* sites = nullptr;
    VectorSiteContainer *allSites = nullptr;

    allSites = SequenceApplicationTools::getSiteContainer(alphabet,this->getParams());

    sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, this->getParams(), "", true,!this->PAR_model_indels_, true, 1);

    AlignmentUtils::checkAlignmentConsistency(*sites);

    delete allSites;

    return sites;

}

bpp::SequenceContainer* CastorApplication::getUnalignedSequences(bpp::Alphabet *alphabet){

    bpp::SequenceContainer *sequences = nullptr;

    // If the user requires the computation of an alignment, then the input file is made of unaligned sequences
    bpp::Fasta seqReader;
    sequences = seqReader.readSequences(PAR_input_sequences_, alphabet);

    return sequences;
}

bpp::SequenceContainer* CastorApplication::getSequences(bpp::Alphabet *alphabet){

    bpp::SequenceContainer *sequences = nullptr;

    try {

        if (this->PAR_alignment_) {
            sequences = this->getUnalignedSequences(alphabet);
        }

    } catch (bpp::Exception &e) {

        LOG(FATAL) << "Error when reading sequence file due to: " << e.message();

    }

    return sequences;
}

bpp::SiteContainer* CastorApplication::getSites(bpp::Alphabet *alphabet){

    bpp::SiteContainer *sites = nullptr;

    try {

        if (!this->PAR_alignment_) {

            sites = this->getAlignedSequences(alphabet);

            ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
        }

    } catch (bpp::Exception &e) {

        LOG(FATAL) << "Error when reading sites file due to: " << e.message();

    }

    return sites;
}

void CastorApplication::renameInternalNodes(bpp::Tree *tree){

    tree->setNodeName(tree->getRootId(), "root");
    UtreeBppUtils::renameInternalNodes(tree);

}

bpp::Tree* CastorApplication::getTree(bpp::SubstitutionModel *smodel,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet,
                                      bpp::Alphabet *alphabetNoGaps,std::map<std::string, std::string> &modelMap,
                                      std::unique_ptr<GeneticCode> &gCode,bpp::DiscreteDistribution *rDist){

    bpp::Tree *tree = nullptr;

    this->PAR_init_tree_opt_ = ApplicationTools::getStringParameter("init.tree", this->getParams(), "user", "", false,1);

    ApplicationTools::displayResult("Initial tree", this->PAR_init_tree_opt_);

    tree = this->getInitTree(smodel,sites,sequences,alphabet,alphabetNoGaps,modelMap,gCode,rDist);

    // If the tree has multifurcation, then resolve it with midpoint rooting
    this->resolveMultifurcation(tree);

    // Rename internal nodes with standard Vxx * where xx is a progressive number
    this->renameInternalNodes(tree);

    // Try to write the current tree to file. This will be overwritten by the optimized tree,
    // but allow to check file existence before running optimization!
    PhylogeneticsApplicationTools::writeTree(*tree, this->getParams());

    // Setting branch lengths?
    this->initBranchLength(tree);

    return tree;
}

tshlib::Utree * CastorApplication::getUtree(bpp::Tree *tree,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                                            UtreeBppUtils::treemap &tm){

    tshlib::Utree *utree = nullptr;

    // Convert the bpp into utree for tree-search engine
    utree = new UtreeBppUtils::Utree();

    UtreeBppUtils::convertTree_b2u(tree, utree, tm);
    if (this->PAR_alignment_) {
        UtreeBppUtils::associateNode2Alignment(sequences, utree);
    } else {
        UtreeBppUtils::associateNode2Alignment(sites, utree);
    }

    DLOG(INFO) << "Bidirectional map size: " << tm.size();
    DLOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true);

    utree->addVirtualRootNode();

    // Once the tree has the root, then map it as well
    tm.insert(UtreeBppUtils::nodeassoc(tree->getRootId(), utree->rootnode->getVnode_id()));

    return utree;
}

bpp::Tree* CastorApplication::getInitTree(bpp::SubstitutionModel *smodel,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet,
                                          bpp::Alphabet *alphabetNoGaps,std::map<std::string, std::string> &modelMap,
                                          std::unique_ptr<GeneticCode> &gCode,bpp::DiscreteDistribution *rDist){

    bpp::Tree *tree = nullptr;

    if (this->PAR_init_tree_opt_ == "user") {

        tree = this->getUserTree();

    } else if (this->PAR_init_tree_opt_ == "random") {

        tree = this->getRandomTree(sites);

    } else if (this->PAR_init_tree_opt_ == "distance") {

        tree = this->getDistanceTree(smodel,sites,sequences,alphabet,alphabetNoGaps,modelMap,gCode,rDist);

    } else{

        throw Exception("Unknown init tree method.");

    }

    return tree;
}

bpp::Tree* CastorApplication::getUserTree(){

    bpp::Tree *tree = nullptr;

    tree = PhylogeneticsApplicationTools::getTree(this->getParams());

    DLOG(INFO) << "[Input tree parser] Number of leaves" << tree->getNumberOfLeaves();

    return tree;
}

bpp::Tree * CastorApplication::getRandomTree(bpp::SiteContainer *sites){

    bpp::Tree *tree = nullptr;

    vector<string> names = sites->getSequencesNames();

    tree = TreeTemplateTools::getRandomTree(names);

    tree->setBranchLengths(1.);

    return tree;
}

bpp::Tree* CastorApplication::getDistanceTree(bpp::TransitionModel *dmodel,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                                              bpp::Alphabet *alphabet,bpp::Alphabet *alphabetNoGaps,std::map<std::string, std::string> &modelMap,
                                              std::unique_ptr<GeneticCode> &gCode,bpp::DiscreteDistribution *rDist){

    bpp::Tree *tree = nullptr;

    ApplicationTools::displayResult("Initial tree reconstruction method", this->PAR_distance_method_);

    AgglomerativeDistanceMethod *distMethod = nullptr;

    std::string token = PAR_distance_method_.substr(0, this->PAR_distance_method_.find("-"));

    if(token == "wpgma") {
        auto *wpgma = new PGMA(true);
        distMethod = wpgma;
        tree = this->infereInitTree(dmodel,alphabet,alphabetNoGaps,sites,sequences,modelMap,gCode,distMethod,rDist);
        delete wpgma;
    }else if(token == "upgma"){
        auto *upgma = new PGMA(false);
        distMethod = upgma;
        tree = this->infereInitTree(dmodel,alphabet,alphabetNoGaps,sites,sequences,modelMap,gCode,distMethod,rDist);
        delete upgma;
    }else if(token == "nj"){
        auto *nj = new NeighborJoining();
        nj->outputPositiveLengths(true);
        distMethod = nj;
        tree = this->infereInitTree(dmodel,alphabet,alphabetNoGaps,sites,sequences,modelMap,gCode,distMethod,rDist);
        delete nj;
    }else if(token == "bionj"){
        auto *bionj = new BioNJ();
        bionj->outputPositiveLengths(true);
        distMethod = bionj;
        tree = this->infereInitTree(dmodel,alphabet,alphabetNoGaps,sites,sequences,modelMap,gCode,distMethod,rDist);
        delete bionj;
    }else if(token == "distmatrix"){
        tree = this->getUserDistmatrixTree();
    }else if(token == "infere_distance_matrix") {
        tree = this->infereDistanceMatrixTree(sequences, alphabet);
    }else{
        throw Exception("Unknown tree reconstruction method.");
    }

    return tree;
}

/*
bpp::Tree* CastorApplication::getDistanceTree(bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet,
                                              bpp::Alphabet *alphabetNoGaps,std::map<std::string, std::string> &modelMap,
                                              std::unique_ptr<GeneticCode> &gCode){

    bpp::Tree *tree = nullptr;

    this->PAR_distance_method_ = ApplicationTools::getStringParameter("init.distance.method",this->getParams(), "nj");

    ApplicationTools::displayResult("Initial tree reconstruction method", PAR_distance_method_);

    AgglomerativeDistanceMethod *distMethod = nullptr;

    std::string token = PAR_distance_method_.substr(0, this->PAR_distance_method_.find("-"));

    if(token == "wpgma") {
        auto *wpgma = new PGMA(true);
        distMethod = wpgma;
        tree = this->infereInitTree(alphabet,alphabetNoGaps,sites,sequences,modelMap,gCode,distMethod);
        delete wpgma;
    }else if(token == "upgma"){
        auto *upgma = new PGMA(false);
        distMethod = upgma;
        tree = this->infereInitTree(alphabet,alphabetNoGaps,sites,sequences,modelMap,gCode,distMethod);
        delete upgma;
    }else if(token == "nj"){
        auto *nj = new NeighborJoining();
        nj->outputPositiveLengths(true);
        distMethod = nj;
        tree = this->infereInitTree(alphabet,alphabetNoGaps,sites,sequences,modelMap,gCode,distMethod);
        delete nj;
    }else if(token == "bionj"){
        auto *bionj = new BioNJ();
        bionj->outputPositiveLengths(true);
        distMethod = bionj;
        tree = this->infereInitTree(alphabet,alphabetNoGaps,sites,sequences,modelMap,gCode,distMethod);
        delete bionj;
    }else if(token == "distmatrix"){
        tree = this->getUserDistmatrixTree();
    }else if(token == "infere_distance_matrix") {
        tree = this->infereDistanceMatrixTree(sequences, alphabet);
    }else{
        throw Exception("Unknown tree reconstruction method.");
    }

    return tree;
}
*/
bpp::Tree* CastorApplication::getUserDistmatrixTree(){

    bpp::Tree* tree;

    // Use a distance matrix provided by the user
    ApplicationTools::displayResult("Initial tree method", std::string("LZ compression"));

    bpp::DistanceMatrix *distances = nullptr;

    try {
        this->PAR_distance_method_ = ApplicationTools::getAFilePath("init.distance.matrix.file",this->getParams(), true, true, "", false,"",0);
    } catch (bpp::Exception &e) {
        LOG(FATAL) << "Error when reading distance matrix file: " << e.message();
    }

    DLOG(INFO) << "initial tree method from LZ compression from matrix file" << this->PAR_distance_matrix_;

    distances = InputUtils::parseDistanceMatrix(this->PAR_distance_matrix_);

    bpp::BioNJ bionj(*distances, true, true, false);

    tree = bionj.getTree();

    delete distances;

    return tree;
}

bpp::Tree* CastorApplication::infereDistanceMatrixTree(bpp::SequenceContainer *sequences, bpp::Alphabet *alphabet){

    bpp::Tree* tree = nullptr;
    bpp::DistanceMatrix *dist_ = nullptr;
    AgglomerativeDistanceMethod *distMethod = nullptr;

    int ALPHABET_DIM = 0;
    int K = 0;
    bool mldist_flag = true;
    bool mldist_gap_flag = false;
    double cutoff_dist = 1.0;
    double indel_rate = 1.0;

    if (this->PAR_alphabet_.find("DNA") != std::string::npos) {
        ALPHABET_DIM = 4;
        K = 6;
    } else if (this->PAR_alphabet_.find("Protein") != std::string::npos) {
        ALPHABET_DIM = 20;
        K = 2;
    } else if (this->PAR_alphabet_.find("Codon") != std::string::npos) {
        ALPHABET_DIM = 60;
        K = 2;
    }

    if (!sequences) {
        bpp::Fasta seqReader;
        sequences = seqReader.readSequences(this->PAR_input_sequences_, alphabet);
    }

    DistanceFactoryPrographMSA::DistanceFactory *dist_factory_angle = new DistanceFactoryPrographMSA::DistanceFactoryAngle(ALPHABET_DIM, K);

    DistanceFactoryPrographMSA::DistanceMatrix dist_ml = dist_factory_angle->computePwDistances(sequences,ALPHABET_DIM,K,mldist_flag,mldist_gap_flag,cutoff_dist,indel_rate);

    dist_ = new DistanceMatrix(sequences->getSequencesNames());

    for (int i = 0; i < sequences->getNumberOfSequences(); i++) {
        for (int j = 0; j < sequences->getNumberOfSequences(); j++) {
            (*dist_)(i, j) = (abs(dist_ml.distances(i, j)) < DISTCUTOFF ? 0.0 : abs(dist_ml.distances(i, j)));
        }
    }

    auto *bionj = new BioNJ(true, true, false);

    bionj->outputPositiveLengths(true);

    distMethod = bionj;

    distMethod->setDistanceMatrix(*dist_);

    distMethod->computeTree();

    tree = distMethod->getTree();

    delete dist_factory_angle;
    delete dist_;
    delete bionj;

    return tree;
}

bpp::Tree * CastorApplication::resolveMultifurcation(bpp::Tree *tree){

    auto ttree_ = new TreeTemplate<Node>(*tree);

    if (ttree_->getRootNode()->getNumberOfSons() > 2) {
        TreeTemplateTools::midRoot(*(ttree_), TreeTemplateTools::MIDROOT_VARIANCE, false);
        tree = ttree_;
    }

    return tree;
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

bpp::SubstitutionModel* CastorApplication::getSubstitutionModel(std::map<std::string, std::string> &modelMap,
                          std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                          bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                          bpp::Tree *tree){

    bpp::SubstitutionModel *smodel;

    // Instantiate a substitution model and extend it with PIP
    if (this->PAR_model_indels_) {

        smodel = this->getSubstitutionModelIndel(modelMap,gCode,alphabet,alphabetNoGaps,sites,sequences,tree);

    } else {

        smodel = this->getSubstitutionModelNoIndel(gCode,sites,alphabet);

    }

    DLOG(INFO) << "[Substitution model] Number of states: " << (int) smodel->getNumberOfStates();

    ApplicationTools::displayResult("Substitution model", smodel->getName());

    if (this->PAR_model_indels_){
        ApplicationTools::displayResult("Indel parameter initial value",(this->PAR_estimate_pip_parameters_) ? "estimated" : "fixed");
    }

    bpp::StdStr stream;
    bpp::PhylogeneticsApplicationTools::printParameters(smodel, stream, 1, true);
    DLOG(INFO) << stream.str();

    return smodel;
}

void CastorApplication::updateModelMap(std::map<std::string, std::string> &modelMap){

    this->PAR_compute_frequencies_from_data_ = false;

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
                    this->PAR_compute_frequencies_from_data_ = true;
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

    this->PAR_estimate_pip_parameters_ = false;

    if (modelMap.find("estimated") != modelMap.end()) {

        ApplicationTools::displayError("The use of the tag [observed] is obsolete. Use the tag [initFreqs] instead");
        exit(EXIT_FAILURE);

    } else if (modelMap.find("initFreqs") != modelMap.end()) {

        if (modelMap["initFreqs"] == "observed") {
            this->PAR_estimate_pip_parameters_ = true;
        }

    } else if( (modelMap.find("lambda") == modelMap.end()) || (modelMap.find("mu") == modelMap.end())) {

        this->PAR_estimate_pip_parameters_ = true;

    }

    if (this->PAR_estimate_pip_parameters_) {

        inference_indel_rates::infere_indel_rates_from_sequences(this->PAR_input_sequences_,this->PAR_alphabet_,tree,this->lambda_,this->mu_,gCode.get(),modelMap);

        DLOG(INFO) << "[PIP model] Estimated PIP parameters from data using input sequences (lambda=" <<this->lambda_ << ",mu=" << this->mu_ << "," "I=" << this->lambda_ * this->mu_ << ")";

    } else {
        this->lambda_ = (modelMap.find("lambda") == modelMap.end()) ? 1.0 : std::stod(modelMap["lambda"]);
        this->mu_ = (modelMap.find("mu") == modelMap.end()) ? 0.1 : std::stod(modelMap["mu"]);
    }

    DLOG(INFO) << "[PIP model] Fixed PIP parameters to (lambda=" << this->lambda_ << ",mu=" << this->mu_ << "," "I="<< this->lambda_ * this->mu_ << ")";

}

bpp::SubstitutionModel* CastorApplication::getSubstitutionModelNoIndel(std::unique_ptr<GeneticCode> &gCode,
                                                                       bpp::SiteContainer *sites,bpp::Alphabet *alphabet){

    bpp::SubstitutionModel *smodel = nullptr;

    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites,this->getParams(), "", true, false, 0);

    return smodel;
}

bpp::SubstitutionModel* CastorApplication::initCanonicalSubstitutionModel(bpp::Alphabet *alphabetNoGaps,std::unique_ptr<GeneticCode> &gCode,
                                                                          bpp::SiteContainer *sites,std::map<std::string, std::string> &modelMap,
                                                                          bpp::Alphabet *alphabet){

    bpp::SubstitutionModel *smodel_canonical = nullptr;

    if (this->PAR_alphabet_.find("Codon") != std::string::npos){

        smodel_canonical = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabetNoGaps, gCode.get(), sites,
                                                                               modelMap, "",true, false, 0);
    }else if (this->PAR_alphabet_.find("Protein") != std::string::npos){

        smodel_canonical = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabetNoGaps, gCode.get(), sites,
                                                                               modelMap, "",true, false, 0);
    }else if (this->PAR_alphabet_.find("DNA") != std::string::npos){

        smodel_canonical = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet,gCode.get(), sites,
                                                                               modelMap,"",true,false, 0);
    }

    return smodel_canonical;

}

bpp::SubstitutionModel* CastorApplication::initPIPsubstitutionModel(bpp::Alphabet *alphabet,bpp::SubstitutionModel *smodel_init,
                                                                    const bpp::SequenceContainer *data,std::unique_ptr<GeneticCode> &gCode){

    bpp::SubstitutionModel *smodel = nullptr;

    if (this->PAR_alphabet_.find("DNA") != std::string::npos && this->PAR_alphabet_.find("Codon") == std::string::npos) {

        smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabet), smodel_init, *data, this->lambda_, this->mu_,this->PAR_compute_frequencies_from_data_);

    } else if (this->PAR_alphabet_.find("Protein") != std::string::npos) {

        smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabet), smodel_init, *data, this->lambda_, this->mu_,this->PAR_compute_frequencies_from_data_);

    } else if (this->PAR_alphabet_.find("Codon") != std::string::npos) {

        smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabet), gCode.get(), smodel_init,*data,this->lambda_, this->mu_,this->PAR_compute_frequencies_from_data_);

        ApplicationTools::displayWarning("Codon models are experimental in the current version... use with caution!");

        DLOG(WARNING) << "CODONS activated byt the program is not fully tested under these settings!";

    }

    return smodel;
}

bpp::SubstitutionModel* CastorApplication::getSubstitutionModelIndel(std::map<std::string, std::string> &modelMap,
                                             std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                                             bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                                             bpp::Tree *tree) {

    bpp::SubstitutionModel *smodel_pip = nullptr;
    bpp::SubstitutionModel *smodel_canonical = nullptr;

    this->updateModelMap(modelMap);

    // Instantiation of the canonical substitution model
    smodel_canonical = this->initCanonicalSubstitutionModel(alphabetNoGaps,gCode,sites,modelMap,alphabet);

    // get insertion rate and deletion rate
    this->getIndelRates(modelMap,tree,gCode);

    const SequenceContainer *data;
    if (this->PAR_alignment_) {
        data = sequences;
    }else{
        data = sites;
    }

    // Instantiate the corrisponding PIP model given the alphabet
    smodel_pip = this->initPIPsubstitutionModel(alphabet,smodel_canonical,data,gCode);

    delete smodel_canonical;

    return smodel_pip;
}

bpp::SubstitutionModel* CastorApplication::getSubstitutionModelIndel(std::map<std::string, std::string> &modelMap,
                                                                     std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                                                                     bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences) {

    bpp::SubstitutionModel *smodel_pip = nullptr;
    bpp::SubstitutionModel *smodel_canonical = nullptr;

    this->updateModelMap(modelMap);

    // Instantiation of the canonical substitution model
    smodel_canonical = this->initCanonicalSubstitutionModel(alphabetNoGaps,gCode,sites,modelMap,alphabet);

    const SequenceContainer *data;
    if (this->PAR_alignment_) {
        data = sequences;
    }else{
        data = sites;
    }

    // Instantiate the corrisponding PIP model given the alphabet
    smodel_pip = this->initPIPsubstitutionModel(alphabet,smodel_canonical,data,gCode);

    delete smodel_canonical;

    return smodel_pip;
}

ParameterList CastorApplication::getParameters(bpp::SubstitutionModel *smodel){

    ParameterList parameters;

    parameters = smodel->getParameters();

    return parameters;
}

DiscreteDistribution* CastorApplication::getASVR(bpp::SubstitutionModel *smodel){

    bpp::DiscreteDistribution *rDist = nullptr;

    // Among site rate variation (ASVR)
    if (smodel->getNumberOfStates() >= 2 * smodel->getAlphabet()->getSize()) {

        // Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();

    } else {

        rDist = PhylogeneticsApplicationTools::getRateDistribution(this->getParams());

    }

    return rDist;
}

bpp::SiteContainer* CastorApplication::getMSA(bpp::SequenceContainer *sequences,bpp::DiscreteDistribution *rDist,
                                              bpp::SubstitutionModel *smodel,UtreeBppUtils::treemap &tm,bpp::Tree *tree,
                                              UtreeBppUtils::Utree *utree){

    /////////////////////////////////
    ApplicationTools::displayMessage("\n[Computing the multi-sequence alignment]");

    ApplicationTools::displayResult("Aligner optimised for:", this->PAR_alignment_version_);

    ApplicationTools::displayBooleanResult("Stochastic backtracking active", this->PAR_alignment_sb_solutions_ > 1);

    if (this->PAR_alignment_sb_solutions_ > 1) {
        ApplicationTools::displayResult("Number of stochastic solutions:",TextTools::toString(this->PAR_alignment_sb_solutions_));
    }

    DLOG(INFO) << "[Alignment sequences] Starting MSA inference using Pro-PIP...";
    /////////////////////////////////

    progressivePIP *proPIP = nullptr;

    bpp::SiteContainer * sites = nullptr;

    // Execute alignment on post-order node list
    std::vector<tshlib::VirtualNode *> ftn = utree->getPostOrderNodeList();

    enumDP3Dversion DPversion = RAM;

    if (PAR_alignment_version_.find("cpu") != std::string::npos) {
        DPversion = CPU; // slower but uses less memory
        this->PAR_alignment_sb_solutions_ = 1;
    } else if (PAR_alignment_version_.find("ram") != std::string::npos) {
        DPversion = RAM; // faster but uses more memory
        this->PAR_alignment_sb_solutions_ = 1;
    } else if (PAR_alignment_version_.find("sb") != std::string::npos) {
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
                           this->PAR_alignment_sb_solutions_,       // number of suboptimal MSAs
                           this->PAR_alignment_sb_temperature_); // to tune the greedyness of the sub-optimal solution

    proPIP->PIPnodeAlign(); // align input sequences with a DP algorithm under PIP

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    ApplicationTools::displayResult("\nAlignment elapsed time (msec):",TextTools::toString((double) duration, 4));

    // convert PIPmsa into a sites objects
    sites = PIPmsaUtils::PIPmsa2Sites(proPIP->alphabet_,*(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getseqNames()),
                                      *(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getMSA()) );


    ApplicationTools::displayResult("Alignment log likelihood", TextTools::toString(proPIP->getPIPnodeRootNode()->MSA_->getMSA()->_getScore(), 15));

    delete proPIP;

    return sites;
}

void CastorApplication::getOptParams(bpp::AbstractHomogeneousTreeLikelihood *tl){

    //==================================
    // m@x
    //tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(Optimizators::optimizeParameters(tl,tl->getParameters(),this->getParams(),"",true,true,0));
    Optimizators * opt = nullptr;
    opt = new Optimizators();
    opt->init(this->getParams(),"", true, true, 0);
    tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(opt->optimizeParameters(tl,tl->getParameters()));
    delete opt;
    //==================================

    this->logL_ = tl->getLogLikelihood();

}

bpp::AbstractHomogeneousTreeLikelihood * CastorApplication::initLK(bpp::TransitionModel *model,bpp::SubstitutionModel *smodel,
            UtreeBppUtils::treemap &tm,bpp::DiscreteDistribution *rDist,std::unique_ptr<GeneticCode> &gCode,
            bpp::Alphabet *alphabet,bpp::SiteContainer *sites,bpp::Tree *tree,UtreeBppUtils::Utree *utree){

    bpp::AbstractHomogeneousTreeLikelihood *tl = nullptr;

    if (!this->PAR_model_indels_) {

        // Get transition model from substitution model
        model = bpp::PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), sites,this->getParams(), "", true, false, 0);

        // Initialise likelihood functions
        tl = new bpp::UnifiedTSHomogeneousTreeLikelihood(*tree, *sites, model, rDist, utree, &tm, true,this->getParams(), "", false, false,false);

    } else {

        // Get transition model from substitution model
        unique_ptr<TransitionModel> test;
        test.reset(smodel);
        model = test.release();

        std::map <std::string, std::string> params = this->getParams();

        // Initialise likelihood functions
        tl = new bpp::UnifiedTSHomogeneousTreeLikelihood_PIP(*tree, *sites, model, rDist, utree, &tm, true,this->getParams(), "", false,false, false);
    }

    ApplicationTools::displayResult("Tree likelihood model", std::string("Homogeneous"));

    tl->initialize();

    return tl;
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
        if (this->PAR_codon_alphabet_) {
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
            if (f){
                exit(EXIT_FAILURE);
            }
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

void CastorApplication::getBootstrap(UtreeBppUtils::Utree *utree,bpp::AbstractHomogeneousTreeLikelihood *tl,
                                     bpp::SiteContainer *sites,UtreeBppUtils::treemap &tm,bpp::DiscreteDistribution *rDist){

    // Compute support measures
    this->PAR_support_ = ApplicationTools::getStringParameter("support", this->getParams(), "", "", true,
                                                              true);
    if (this->PAR_support_ == "bootstrap") {
        ApplicationTools::displayMessage("\n[Tree support measures]");

        bpp::Bootstrap(tl, *sites, rDist, utree, &tm, this->getParams(), "support.");
    }

}

void CastorApplication::output(bpp::Tree *tree,UtreeBppUtils::Utree *utree,bpp::SiteContainer *sites,
            bpp::AbstractHomogeneousTreeLikelihood *tl,UtreeBppUtils::treemap &tm,ParameterList &parameters,
            bpp::DiscreteDistribution *rDist){

    // Export alignment to file
    if (this->PAR_output_file_msa_.find("none") == std::string::npos) {

        DLOG(INFO) << "[Alignment sequences]\t The final alignment can be found in " << PAR_output_file_msa_;

        bpp::Fasta seqWriter;
        seqWriter.writeAlignment(TextUtils::appendToFilePath(this->PAR_output_file_msa_, "initial"), *sites, true);
    }

    // Export final tree (if nexus is required, then our re-implementation of the the nexus writer is called)
    tree = new TreeTemplate<Node>(tl->getTree());

    this->PAR_output_tree_format_ = ApplicationTools::getStringParameter("output.tree.format",this->getParams(),"Newick","",true,true);

    if (this->PAR_output_tree_format_.find("Nexus") != std::string::npos) {
        std::vector<Tree *> tmp;
        tmp.push_back(tree);
        OutputUtils::writeNexusMetaTree(tmp, this->getParams());
    } else {
        PhylogeneticsApplicationTools::writeTree(*tree, this->getParams());
    }

    // Export annotation file (tab separated values)
    this->PAR_output_annotation_file_ = ApplicationTools::getAFilePath("output.annotation.file",this->getParams(),false,false,"",true,"",1);

    if (this->PAR_output_annotation_file_.find("none") == std::string::npos) {
        ApplicationTools::displayResult("Output annotation to file", PAR_output_annotation_file_);

        OutputUtils::writeTreeAnnotations2TSV(tree, PAR_output_annotation_file_);
    }

    // Write parameters to screen:
    ApplicationTools::displayResult("Final Log likelihood", TextTools::toString(logL_, 15));

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

}

void CastorApplication::getIgnoreParams(ParameterList &parametersToIgnore,bool &ignoreBrLen,bpp::TransitionModel *dmodel,bpp::DiscreteDistribution *rDist){

    ParameterList allParameters;

    allParameters = dmodel->getParameters();

    allParameters.addParameters(rDist->getParameters());

    string paramListDesc = ApplicationTools::getStringParameter("init.distance.optimization.ignore_parameter", this->getParams(),"", "", true, false);

    ignoreBrLen = false;

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
            ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
        }
    }

}

TransitionModel* CastorApplication::getTransitionModel(bpp::SubstitutionModel *smodel,bpp::Alphabet *alphabetDistMethod,
                                                       std::unique_ptr<GeneticCode> &gCode,VectorSiteContainer *sitesDistMethod,
                                                       map<std::string, std::string> &parmap){

    TransitionModel * dmodel;

    // Get transition model from substitution model
    if (this->PAR_model_indels_) {
        unique_ptr<TransitionModel> test;
        test.reset(smodel);
        dmodel = test.release();
    } else {
        dmodel = PhylogeneticsApplicationTools::getTransitionModel(alphabetDistMethod, gCode.get(),sitesDistMethod, parmap);
    }

    return dmodel;
}


bpp::Tree * CastorApplication::initTreeDistanceOptML(bpp::AgglomerativeDistanceMethod *distMethod,bool ignoreBrLen,
                                                     bpp::TransitionModel *dmodel,bpp::ParameterList &parametersToIgnore,
                                                     bpp::DiscreteDistribution *rDist,VectorSiteContainer *sitesDistMethod){

    bpp::Tree *tree = nullptr;
    Optimizators * opt = nullptr;

    UnifiedDistanceEstimation distEstimation(dmodel, rDist, sitesDistMethod, 1, false);

    opt = new Optimizators();

    // Optimisation method verbosity
    opt->setOptVerbose(ApplicationTools::getParameter<unsigned int>("optimization.verbose",this->getParams(), 2));

    opt->setMessageHandler();

    opt->setProfiler();

    // Should I ignore some parameters?
    this->getIgnoreParams(parametersToIgnore,ignoreBrLen,dmodel,rDist);

    opt->setNumEvalMax(ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval",this->getParams(), 1000000));

    ApplicationTools::displayResult("Initial tree optimization | max # ML evaluations",TextTools::toString(opt->getNumEvalMax()));

    opt->setTolerance(ApplicationTools::getDoubleParameter("optimization.tolerance",this->getParams(), .000001));

    ApplicationTools::displayResult("Initial tree optimization | Tolerance",TextTools::toString(opt->getTolerance()));

    tree = opt->buildDistanceTreeGeneric(distEstimation, *distMethod, parametersToIgnore,!ignoreBrLen,this->PAR_optim_distance_);

    delete opt;

    return tree;
}

bpp::Tree* CastorApplication::initTreeDistanceNoOpt(bpp::AgglomerativeDistanceMethod *distMethod,bpp::TransitionModel *dmodel,
                                                    bpp::DiscreteDistribution *rDist,VectorSiteContainer *sitesDistMethod){

    bpp::Tree *tree = nullptr;
    DistanceMatrix *dm = nullptr;

    UnifiedDistanceEstimation distEstimation(dmodel, rDist, sitesDistMethod, 1, false);

    distEstimation.computeMatrix();

    dm = distEstimation.getMatrix();

    distMethod->setDistanceMatrix((*dm));

    distMethod->computeTree();

    tree = distMethod->getTree();

    //delete dm;

    return tree;
}
/*
bpp::Tree* CastorApplication::infereInitTree(bpp::Alphabet *alphabet,bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,
                                              bpp::SequenceContainer *sequences,std::map<std::string, std::string> &modelMap,
                                              std::unique_ptr<GeneticCode> &gCode,bpp::AgglomerativeDistanceMethod *distMethod){
    bpp::Tree *tree = nullptr;

    //bpp::DiscreteDistribution *rDist = nullptr;
    VectorSiteContainer *allSites = nullptr;
    VectorSiteContainer *sitesDistMethod = nullptr;
    bpp::Alphabet *alphabetDistMethod = nullptr;
    //bpp::SubstitutionModel *smodel;
    //bpp::TransitionModel *dmodel = nullptr;
    bpp::ParameterList parametersToIgnore;
    bool ignoreBrLen = false;
    map<std::string, std::string> parmap;

    if (this->PAR_model_indels_) {
        alphabetDistMethod = alphabet;
        parmap["model"] = modelMap["model"];
    } else {
        alphabetDistMethod = alphabetNoGaps;
        parmap["model"] = this->getParams()["model"];
    }

    allSites = SequenceApplicationTools::getSiteContainer(alphabetDistMethod, this->getParams());
    sitesDistMethod = SequenceApplicationTools::getSitesToAnalyse(*allSites, this->getParams());

    //smodel = this->getSubstitutionModel(modelMap,gCode,alphabet,alphabetNoGaps,sites,sequences,tree);

    //Initialize model to compute the distance tree
    //dmodel = this->getTransitionModel(smodel,alphabetDistMethod,gCode,sitesDistMethod,parmap);

    // Add a ASRV distribution
    //rDist = this->getASVR(smodel);

    // Remove gap characters since we are roughly estimating the initial topology
    if (!this->PAR_model_indels_) {
        bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sitesDistMethod);
    }

    if (this->PAR_distance_method_.find("-ml") != std::string::npos) {
        // Distance based initial tree topology inference  with optimisation
        tree = initTreeDistanceOptML(distMethod,ignoreBrLen,dmodel,parametersToIgnore,rDist,sitesDistMethod);
    }else{
        // Fast but rough estimate of the initial tree topology (distance based without optimisation -ML)
        tree = initTreeDistanceNoOpt(distMethod,dmodel,rDist,sitesDistMethod);
    }

    delete allSites;
    delete sitesDistMethod;
    delete alphabetDistMethod;
    //delete smodel;
    //delete dmodel;
    //delete rDist;

    return tree;
}
*/
bpp::Tree* CastorApplication::infereInitTree(bpp::TransitionModel *dmodel,bpp::Alphabet *alphabet,bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,
                                             bpp::SequenceContainer *sequences,std::map<std::string, std::string> &modelMap,
                                             std::unique_ptr<GeneticCode> &gCode,bpp::AgglomerativeDistanceMethod *distMethod,
                                             bpp::DiscreteDistribution *rDist){

    bpp::Tree *tree = nullptr;

    //bpp::DiscreteDistribution *rDist = nullptr;
    VectorSiteContainer *allSites = nullptr;
    VectorSiteContainer *sitesDistMethod = nullptr;
    bpp::Alphabet *alphabetDistMethod = nullptr;
    //bpp::TransitionModel *dmodel = nullptr;
    bpp::ParameterList parametersToIgnore;
    bool ignoreBrLen = false;
    map<std::string, std::string> parmap;

    if (this->PAR_model_indels_) {
        alphabetDistMethod = alphabet;
        parmap["model"] = modelMap["model"];
    } else {
        alphabetDistMethod = alphabetNoGaps;
        parmap["model"] = this->getParams()["model"];
    }

    allSites = SequenceApplicationTools::getSiteContainer(alphabetDistMethod, this->getParams());
    sitesDistMethod = SequenceApplicationTools::getSitesToAnalyse(*allSites, this->getParams());

    //smodel = this->getSubstitutionModel(modelMap,gCode,alphabet,alphabetNoGaps,sites,sequences,tree);

    //Initialize model to compute the distance tree
    //dmodel = this->getTransitionModel(smodel,alphabetDistMethod,gCode,sitesDistMethod,parmap);

    // Add a ASRV distribution
    //rDist = this->getASVR(smodel);

    // Remove gap characters since we are roughly estimating the initial topology
    if (!this->PAR_model_indels_) {
        bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sitesDistMethod);
    }

    if (this->PAR_distance_method_.find("-ml") != std::string::npos) {
        // Distance based initial tree topology inference  with optimisation
        tree = initTreeDistanceOptML(distMethod,ignoreBrLen,dmodel,parametersToIgnore,rDist,sitesDistMethod);
    }else{
        // Fast but rough estimate of the initial tree topology (distance based without optimisation -ML)
        tree = initTreeDistanceNoOpt(distMethod,dmodel,rDist,sitesDistMethod);
    }

    //delete allSites;
    //delete sitesDistMethod;
    //delete alphabetDistMethod;

    //delete smodel;
    //delete dmodel;
    //delete rDist;

    return tree;
}

void CastorApplicationUtils::printParameters(ParameterList &parameters,bpp::SubstitutionModel *smodel){

    //********************************************************//
    // print parameters
    for (size_t i = 0; i < parameters.size(); i++) {
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }

    for (size_t i = 0; i < smodel->getFrequencies().size(); i++) {
        ApplicationTools::displayResult("eq.freq(" + smodel->getAlphabet()->getName(i) + ")",TextTools::toString(smodel->getFrequencies()[i], 4));
    }
    //********************************************************//

}