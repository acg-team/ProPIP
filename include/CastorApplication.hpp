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

#include <boost/asio/ip/host_name.hpp>

#include <Bpp/Text/KeyvalTools.h>

#include <memory>

/*
* From SeqLib:
*/
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include "ExtendedAlphabet.hpp"
#include "Utils.hpp"

#include <Bpp/Phyl/Distance/DistanceMethod.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Distance/PGMA.h>

#include "DistanceFactory.hpp"
#include "DistanceFactoryAngle.hpp"
#include "DistanceFactoryAlign.hpp"
#include "DistanceFactoryPrealigned.hpp"

#include "inference_indel_rates.hpp"

#include "PIP.hpp"
#include "ExtendedAlphabet.hpp"

#include "Optimizators.hpp"
#include "UnifiedDistanceEstimation.hpp"

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

#include <Utree.hpp>
#include <TreeRearrangment.hpp>

#include "progressivePIP.hpp"
#include "FactoryPIPnode.hpp"
#include "CompositePIPnode.hpp"

#include "SupportMeasures.hpp"

#include "RHomogeneousTreeLikelihood_PIP.hpp"
#include "RHomogeneousTreeLikelihood_Generic.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_PIP.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_Generic.hpp"

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
        CastorApplication(int argc, char *argv[], const std::string &name, const std::string &strVersion, const std::string &build_date);

    public:

        std::string PAR_model_substitution;
        bool PAR_alignment;
        bool PAR_model_indels;
        std::string PAR_Alphabet;
        std::string PAR_input_sequences;
        std::string PAR_distance_method;
        std::string PAR_output_file_msa;
        std::string PAR_alignment_version;
        std::string PAR_output_tree_format;
        std::string PAR_output_annotation_file;
        std::string PAR_support;
        std::string PAR_initTreeOpt;

        int PAR_alignment_sbsolutions;
        double PAR_alignment_sbtemperature;
        bool PAR_estimatePIPparameters;
        bool PAR_computeFrequenciesFromData;

        bool codonAlphabet;

        double lambda;
        double mu;
        double logL;

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
            std::cout << "Usage: ProPIP [arguments] or [params=file.txt]" << std::endl;
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

        void checkInputArgument(int argc) {

            if (argc < 2) {
                this->help();
                exit(0);
            } else {
                this->banner();
                this->startTimer();
            }

        }

        void getCLIarguments(std::map<std::string, std::string> &modelMap);

        void getAlphabet(bpp::Alphabet *alphabetNoGaps,std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet);

        void getData(bpp::SequenceContainer *sequences,bpp::SiteContainer *sites,bpp::Alphabet *alphabet);

        void getInitTree(bpp::Tree *tree,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet);

        void initBranchLength(bpp::Tree *tree);

        void resolveMultifurcation(bpp::Tree *tree);

        void getUserTree(bpp::Tree *tree);

        void getRandomTree(bpp::Tree *tree,bpp::SiteContainer *sites);

        void getDistanceTree(bpp::Tree *tree,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet);

        void getDistmatrixTree(bpp::Tree *tree);

        void getInfere_distance_matrix_Tree(bpp::Tree *tree,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet);

        void getTree(bpp::Alphabet *alphabet,bpp::Alphabet *alphabetNoGaps,bpp::Tree *tree,
                     bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,std::map<std::string, std::string> &modelMap,
                     std::unique_ptr<GeneticCode> &gCode,UtreeBppUtils::treemap &tm,tshlib::Utree *utree);

        /*
        void check_this_code(bpp::Alphabet *alphabet,bpp::Alphabet *alphabetNoGaps,bpp::Tree *tree,
                                                bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,std::map<std::string, std::string> &modelMap,
                                                std::unique_ptr<GeneticCode> &gCode,UtreeBppUtils::treemap &tm,tshlib::Utree *utree,
                                                bpp::DistanceMethod *distMethod);
        */

        void getModelMap(std::map<std::string, std::string> &modelMap,bpp::SubstitutionModel *smodel);

        void getIndelRates(std::map<std::string, std::string> &modelMap,bpp::Tree *tree,std::unique_ptr<GeneticCode> &gCode);

        void getSubstitutionModel(std::map<std::string, std::string> &modelMap,bpp::SubstitutionModel *smodel,
                                  bpp::TransitionModel *model,std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                                  bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                                  bpp::Tree *tree);

        void getSubstitutionModelIndel(std::map<std::string, std::string> &modelMap,bpp::SubstitutionModel *smodel,
                                                          bpp::TransitionModel *model,std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                                                          bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                                                          bpp::Tree *tree);

        void getParameters(ParameterList &parameters,bpp::SubstitutionModel *smodel);

        void getASVR(bpp::DiscreteDistribution *rDist,bpp::SubstitutionModel *smodel);

        void getMSA(progressivePIP *proPIP,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::DiscreteDistribution *rDist,
                    bpp::SubstitutionModel *smodel,UtreeBppUtils::treemap &tm,bpp::Tree *tree,UtreeBppUtils::Utree *utree);

        void getOptParams(bpp::AbstractHomogeneousTreeLikelihood *tl);

        void initLK(bpp::TransitionModel *model,bpp::SubstitutionModel *smodel,bpp::AbstractHomogeneousTreeLikelihood *tl,
                    UtreeBppUtils::treemap &tm,bpp::DiscreteDistribution *rDist,std::unique_ptr<GeneticCode> &gCode,
                    bpp::Alphabet *alphabet,bpp::SiteContainer *sites,bpp::Tree *tree,UtreeBppUtils::Utree *utree);

        void getParSanityCheck(bpp::AbstractHomogeneousTreeLikelihood *tl,bpp::SiteContainer *sites,std::unique_ptr<GeneticCode> &gCode);

        void output(bpp::Tree *tree,UtreeBppUtils::Utree *utree,bpp::SiteContainer *sites,bpp::AbstractHomogeneousTreeLikelihood *tl,
                    UtreeBppUtils::treemap &tm,ParameterList &parameters,bpp::DiscreteDistribution *rDist);

    };


} //end of namespace bpp;



#endif //CASTOR_JATIAPPLICATION_HPP
