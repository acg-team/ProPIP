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

namespace bpp {
    class CastorApplication {
    private:
        std::string appName_;
        std::string appBuild_;
        std::string appVersion_;
        mutable std::map<std::string, std::string> params_;
        bool timerStarted_;
        long seed_;

    public:


        std::string PAR_model_substitution;
        std::string modelStringName;
        std::map<std::string, std::string> modelMap;
        bool PAR_model_indels;
        bool PAR_alignment;
        std::string PAR_Alphabet;
        bpp::Alphabet *alphabetNoGaps;
        std::unique_ptr<bpp::GeneticCode> gCode;
        bool codonAlphabet;
        bpp::Alphabet *alphabet;
        std::string codeDesc;
        std::string PAR_input_sequences;
        bpp::SequenceContainer *sequences;
        bpp::SiteContainer *sites;
        bpp::DiscreteDistribution *rDist;
        pPIP *alignment = nullptr;
        progressivePIP *proPIP = nullptr;
        std::string PAR_output_file_msa;
        std::string PAR_alignment_version;
        int PAR_alignment_sbsolutions;
        double PAR_alignment_sbtemperature;
        int num_sb;
        double temperature;
        bpp::AbstractHomogeneousTreeLikelihood *tl;
        bpp::SubstitutionModel *smodel;
        bpp::TransitionModel *model;
        double logL;
        std::string paramNameFile;
        std::string PAR_support;
        std::string PAR_output_tree_format;
        std::string PAR_output_annotation_file;
        bool estimatePIPparameters;
        bool computeFrequenciesFromData;
        std::string baseModel;
        std::map<std::string, std::string> basemodelMap;
        double lambda;
        double mu;
        std::string initTreeOpt;
        bpp::Tree *tree;
        std::string PAR_distance_method;
        std::string PAR_distance_matrix;
        std::string PAR_optim_distance;
        UtreeBppUtils::treemap tm;
        UtreeBppUtils::Utree *utree;
        enumDP3Dversion DPversion;
        std::string initBrLenMethod;

        CastorApplication(int argc, char *argv[], const std::string &name, const std::string &strVersion, const std::string &build_date);

    public:
        void startTimer();

        void done();

        std::map<std::string, std::string> &getParams() { return params_; }

        const std::string &getParam(const std::string &name) const {
            if (params_.find(name) == params_.end()) throw bpp::Exception("BppApplication::getParam(). Parameter '" + name + "' not found.");
            return params_[name];
        }

        std::string &getParam(const std::string &name) { return params_[name]; }

        long getSeed() {return seed_;}

        void help() {
            std::cout << appName_ << std::endl << std::endl;
            std::cout << "Usage: Castor [arguments] or [params=file.txt]" << std::endl;
            std::cout << "Documentation can be found at https://bitbucket.org/lorenzogatti89/castor/" << std::endl;
        }


        void banner() {

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

        void version() {
            std::cout << appName_ << std::endl;
            std::cout << appVersion_ << std::endl;
            std::cout << appBuild_ << std::endl;
        }


        void getCLIarguments();

        void start(int argc);

        void getAlphabet();

        void getAlphabetIndel();

        void getAlphabetNoIndel();

        void getData();

        void initTreeMethod();

        void resolveMultifurcation();

        void renameTreeNodes();

        void bpp2utree();

        void setTreeBranchLengths();

        void setTreeBranchLengthsInput(std::map<std::string, std::string> &cmdArgs);

        void setTreeBranchLengthsEqual(std::map<std::string, std::string> &cmdArgs);

        void setTreeBranchLengthsClock();

        void setTreeBranchLengthsGrafen(std::map<std::string, std::string> &cmdArgs);

        void getTree();

        void getUnalignedSequences();

        void getAlignedSequences();

        void getSubstitutionIndelModel();

        void getSubstitutionNoIndelModel();

        void initCanonicalSubstitutionModel();

        void instantiateCanonicalSubstitutionModel();

        void instantiatePIPSubstitutionModel(double lambda,double mu);

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

    };


} //end of namespace bpp;



#endif //CASTOR_JATIAPPLICATION_HPP
