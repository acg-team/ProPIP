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
 * @file CastorApplication.hpp
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
 * @see For more information visit:
 */
#ifndef CASTOR_CASTORAPPLICATION_HPP
#define CASTOR_CASTORAPPLICATION_HPP
// From the STL:
#include <string>
#include <map>
#include <Bpp/Exceptions.h>
#include <glog/logging.h>
#include <iostream>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <boost/asio/ip/host_name.hpp>
#include <Bpp/Seq/AlphabetIndex/DefaultNucleotideScore.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Distance/PGMA.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include "pPIP.hpp"
#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"
#include "CompositePIPnode.hpp"

#define MIN_BRANCH_LEN 0.000001

namespace bpp {

    class CastorApplication {

    private:

        mutable std::map<std::string, std::string> params_;

        bool timerStarted_;
        bool codonAlphabet_;
        bool estimatePIPparameters;
        bool computeFrequenciesFromData;

        double lambda;
        double mu;
        double logL;
        double temperature;

        long seed_;

        int num_sb;

        enumDP3Dversion DPversion;

        std::string appName_;
        std::string appBuild_;
        std::string appVersion_;
        std::string modelStringName;
        std::string paramNameFile;
        std::string initBrLenMethod;
        std::string baseModel;
        std::string initTreeOpt;

        std::map<std::string, std::string> modelMap;
        std::map<std::string, std::string> basemodelMap;

        bpp::DistanceMatrix *distances;
        bpp::AgglomerativeDistanceMethod *distMethod;

        bpp::SubstitutionModel *smodel;
        bpp::TransitionModel *model;

        std::unique_ptr<bpp::GeneticCode> gCode;

        pPIP *alignment;

    public:

        bool PAR_model_indels;
        bool PAR_alignment;

        double PAR_alignment_sbtemperature;

        int PAR_alignment_sbsolutions;

        std::string PAR_model_substitution;
        std::string PAR_Alphabet;
        std::string PAR_input_sequences;
        std::string PAR_support;
        std::string PAR_output_tree_format;
        std::string PAR_output_annotation_file;
        std::string PAR_distance_method;
        std::string PAR_distance_matrix;
        std::string PAR_optim_distance;
        std::string PAR_output_file_msa;
        std::string PAR_alignment_version;
        std::string codeDesc;

        bpp::Alphabet *alphabetNoGaps;
        bpp::Alphabet *alphabet;

        bpp::AbstractHomogeneousTreeLikelihood *tl;

        bpp::SequenceContainer *sequences;
        bpp::SiteContainer *sites;
        bpp::DiscreteDistribution *rDist;
        bpp::Tree *tree;

        progressivePIP *proPIP;

        UtreeBppUtils::Utree *utree;
        UtreeBppUtils::treemap tm;

    public:

        CastorApplication(int argc, char *argv[], const std::string &name, const std::string &strVersion, const std::string &build_date);

        void init();

        void startTimer();

        void done();

        std::map<std::string, std::string> &getParams() { return params_; }

        const std::string &getParam(const std::string &name) const;

        std::string &getParam(const std::string &name) { return params_[name]; }

        long getSeed() {return seed_;}

        void help();

        void banner();

        void version();

        void getCLIarguments();

        void start(int argc);

        void getAlphabet();

        void getAlphabetIndel();

        void getAlphabetNoIndel();

        void getData();

        void initTreeMethod();

        void resolveMultifurcations();

        void renameTreeNodes();

        void initTreeBranchLength();

        void initTreeBranchLengthInput(std::map<std::string, std::string> &cmdArgs);

        void initTreeBranchLengthEqual(std::map<std::string, std::string> &cmdArgs);

        void initTreeBranchLengthClock();

        void initTreeBranchLengthGrafen(std::map<std::string, std::string> &cmdArgs);

        void convertBppTree2Utree();

        void initTreeMethodUser();

        void initTreeMethodRandom();

        void initTreeMethodDistance();

        void initTreeMethodDistanceWPGMA();

        void initTreeMethodDistanceUPGMA();

        void initTreeMethodDistanceNJ();

        void initTreeMethodDistanceBIONJ();

        void initTreeMethodDistanceDistmatrix();

        void initTreeMethodDistanceInfereDistanceMatrix();

        void infereDistanceTree();

        void infereDistanceTreeFast(bpp::TransitionModel *local_dmodel,
                                    bpp::VectorSiteContainer *local_sitesDistMethod,
                                    DiscreteDistribution *local_rDist);

        void infereDistanceTreeML(bpp::TransitionModel *local_dmodel,
                                  DiscreteDistribution *local_rDist,
                                  bpp::VectorSiteContainer *local_sitesDistMethod);

        void removeGaps(bpp::VectorSiteContainer *local_sitesDistMethod);

        void addASRVdistribution(DiscreteDistribution *local_rDist,
                                 bpp::TransitionModel *local_dmodel);

        void getUsersIndelRates();

        void getTree();

        void getUnalignedSequences();

        void getAlignedSequences();

        void getSubstitutionIndelModel();

        void getSubstitutionNoIndelModel();

        void extendSubstitutionModelWithPIP(const SequenceContainer *data);

        void initCanonicalSubstitutionModel();

        void getIndelRates();

        void getBackgroundFrequencies();

        void getSubstitutionModel();

        void getASRV();

        void computeMSA();

        void initLkFun();

        void optimizeParameters();

        void parameterSanityCheck();

        void bootstrapping();

        void output();

        bpp::TransitionModel * getTransitionModelFromSubsModel(bool PAR_model_indels,
                                                               bpp::SubstitutionModel *smodel,
                                                               const Alphabet* alphabet,
                                                               const GeneticCode* gCode,
                                                               const SiteContainer* data,
                                                               std::map<std::string, std::string>& params,
                                                               const std::string& suffix,
                                                               bool suffixIsOptional,
                                                               bool verbose,
                                                               int warn);

        bpp::ParameterList getParametersList();

        double checkLkValue(bpp::ParameterList &pl);

        void resolveZeroLKValue();

        void checkStopCodon();

        void removeSaturatedSite();

    };


} //end of namespace bpp;



#endif //CASTOR_JATIAPPLICATION_HPP
