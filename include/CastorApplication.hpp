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

        bpp::SubstitutionModel *smodel;

        std::string PAR_model_substitution_;
        std::string PAR_alphabet_;
        std::string PAR_input_sequences_;
        std::string PAR_distance_method_;
        std::string PAR_distance_matrix_;
        std::string PAR_optim_distance_;
        std::string PAR_output_file_msa_;
        std::string PAR_alignment_version_;
        std::string PAR_output_tree_format_;
        std::string PAR_output_annotation_file_;
        std::string PAR_support_;
        std::string PAR_init_tree_opt_;

        int PAR_alignment_sb_solutions_;

        double PAR_alignment_sb_temperature_;

        bool PAR_estimate_pip_parameters_;
        bool PAR_compute_frequencies_from_data_;
        bool PAR_alignment_;
        bool PAR_model_indels_;

        bool PAR_codon_alphabet_;

        double lambda_;
        double mu_;
        double logL_;

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

        void checkInputArgument(int argc);

        void getParseProcedure(std::map<std::string, std::string> &modelMap,std::string &modelStringName);

        bool isPIP(std::string modelStringName);

        std::map<std::string, std::string> getCLIarguments();

        bpp::Alphabet* getAlphabetIndel();

        bpp::Alphabet* getAlphabetNoIndel();

        bpp::Alphabet* getAlphabetNoGaps();

        bpp::Alphabet* getAlphabet();

        std::unique_ptr<GeneticCode> getGcode(bpp::Alphabet *alphabetNoGaps);

        bpp::SiteContainer* getAlignedSequences(bpp::Alphabet *alphabet);

        bpp::SequenceContainer* getUnalignedSequences(bpp::Alphabet *alphabet);

        bpp::SequenceContainer* getSequences(bpp::Alphabet *alphabet);

        bpp::SiteContainer* getSites(bpp::Alphabet *alphabet);

        bpp::Tree* getInitTree(bpp::SubstitutionModel *smodel,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet,
                               bpp::Alphabet *alphabetNoGaps,std::map<std::string, std::string> &modelMap,
                               std::unique_ptr<GeneticCode> &gCode,bpp::DiscreteDistribution *rDist);

        void initBranchLength(bpp::Tree *tree);

        bpp::Tree* resolveMultifurcation(bpp::Tree *tree);

        bpp::Tree* getUserTree();

        bpp::Tree* getRandomTree(bpp::SiteContainer *sites);

        //=======================================//
        // overlapping functions
        bpp::Tree* getDistanceTree(bpp::TransitionModel *dmodel,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet,
                                   bpp::Alphabet *alphabetNoGaps,std::map<std::string, std::string> &modelMap,
                                   std::unique_ptr<GeneticCode> &gCode,bpp::DiscreteDistribution *rDist);

        //bpp::Tree* getDistanceTree(bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet,
        //                                              bpp::Alphabet *alphabetNoGaps,std::map<std::string, std::string> &modelMap,
        //                                              std::unique_ptr<GeneticCode> &gCode);
        //=======================================//

        bpp::Tree* getUserDistmatrixTree();

        bpp::Tree* infereDistanceMatrixTree(bpp::SequenceContainer *sequences, bpp::Alphabet *alphabet);

        void renameInternalNodes(bpp::Tree *tree);

        //=======================================//
        // overlapping functions
        bpp::Tree * infereInitTree(bpp::TransitionModel *dmodel,bpp::Alphabet *alphabet,bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,
                                   bpp::SequenceContainer *sequences,std::map<std::string, std::string> &modelMap,
                                   std::unique_ptr<GeneticCode> &gCode,bpp::AgglomerativeDistanceMethod *distMethod,
                                   bpp::DiscreteDistribution *rDist);

        //bpp::Tree* infereInitTree(bpp::Alphabet *alphabet,bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,
        //                                             bpp::SequenceContainer *sequences,std::map<std::string, std::string> &modelMap,
        //                                             std::unique_ptr<GeneticCode> &gCode,bpp::AgglomerativeDistanceMethod *distMethod);
        //=======================================//

        bpp::Tree* getTree(bpp::SubstitutionModel *smodel,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,bpp::Alphabet *alphabet,
                           bpp::Alphabet *alphabetNoGaps,std::map<std::string, std::string> &modelMap,
                           std::unique_ptr<GeneticCode> &gCode,bpp::DiscreteDistribution *rDist);

        tshlib::Utree* getUtree(bpp::Tree *tree,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                                UtreeBppUtils::treemap &tm);

        void updateModelMap(std::map<std::string, std::string> &modelMap);

        void getIndelRates(std::map<std::string, std::string> &modelMap,bpp::Tree *tree,std::unique_ptr<GeneticCode> &gCode);

        bpp::SubstitutionModel* getSubstitutionModel(std::map<std::string, std::string> &modelMap,
                                  std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                                  bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                                  bpp::Tree *tree);

        bpp::SubstitutionModel* initCanonicalSubstitutionModel(bpp::Alphabet *alphabetNoGaps,std::unique_ptr<GeneticCode> &gCode,
                                                               bpp::SiteContainer *sites,std::map<std::string, std::string> &modelMap,
                                                               bpp::Alphabet *alphabet);

        bpp::SubstitutionModel* initPIPsubstitutionModel(bpp::Alphabet *alphabet,bpp::SubstitutionModel *smodel_init,
                                                                            const bpp::SequenceContainer *data,std::unique_ptr<GeneticCode> &gCode);

        //=======================================//
        // overlapping functions
        bpp::SubstitutionModel* getSubstitutionModelIndel(std::map<std::string, std::string> &modelMap,
                                                                             std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                                                                             bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences,
                                                                             bpp::Tree *tree);

        bpp::SubstitutionModel* getSubstitutionModelIndel(std::map<std::string, std::string> &modelMap,
                                                                             std::unique_ptr<GeneticCode> &gCode,bpp::Alphabet *alphabet,
                                                                             bpp::Alphabet *alphabetNoGaps,bpp::SiteContainer *sites,bpp::SequenceContainer *sequences);
        //=======================================//

        bpp::SubstitutionModel* getSubstitutionModelNoIndel(std::unique_ptr<GeneticCode> &gCode,bpp::SiteContainer *sites,bpp::Alphabet *alphabet);

        ParameterList getParameters(bpp::SubstitutionModel *smodel);

        DiscreteDistribution* getASVR(bpp::SubstitutionModel *smodel);

        bpp::SiteContainer* getMSA(bpp::SequenceContainer *sequences,bpp::DiscreteDistribution *rDist,bpp::SubstitutionModel *smodel,
                                   UtreeBppUtils::treemap &tm,bpp::Tree *tree,UtreeBppUtils::Utree *utree);

        void getOptParams(bpp::AbstractHomogeneousTreeLikelihood *tl);

        bpp::AbstractHomogeneousTreeLikelihood* initLK(bpp::TransitionModel *model,bpp::SubstitutionModel *smodel,
                    UtreeBppUtils::treemap &tm,bpp::DiscreteDistribution *rDist,std::unique_ptr<GeneticCode> &gCode,
                    bpp::Alphabet *alphabet,bpp::SiteContainer *sites,bpp::Tree *tree,UtreeBppUtils::Utree *utree);

        void getParSanityCheck(bpp::AbstractHomogeneousTreeLikelihood *tl,bpp::SiteContainer *sites,std::unique_ptr<GeneticCode> &gCode);

        void getBootstrap(UtreeBppUtils::Utree *utree,bpp::AbstractHomogeneousTreeLikelihood *tl,
                                             bpp::SiteContainer *sites,UtreeBppUtils::treemap &tm,bpp::DiscreteDistribution *rDist);

        void output(bpp::Tree *tree,UtreeBppUtils::Utree *utree,bpp::SiteContainer *sites,bpp::AbstractHomogeneousTreeLikelihood *tl,
                    UtreeBppUtils::treemap &tm,ParameterList &parameters,bpp::DiscreteDistribution *rDist);


        TransitionModel* getTransitionModel(bpp::SubstitutionModel *smodel,bpp::Alphabet *alphabetDistMethod,
                                                               std::unique_ptr<GeneticCode> &gCode,VectorSiteContainer *sitesDistMethod,
                                                               map<std::string, std::string> &parmap);

        void getIgnoreParams(ParameterList &parametersToIgnore,bool &ignoreBrLen,bpp::TransitionModel *dmodel,bpp::DiscreteDistribution *rDist);

        bpp::Tree* initTreeDistanceOptML(bpp::AgglomerativeDistanceMethod *distMethod,bool ignoreBrLen,
                                                             bpp::TransitionModel *dmodel,bpp::ParameterList &parametersToIgnore,
                                                             bpp::DiscreteDistribution *rDist,VectorSiteContainer *sitesDistMethod);

        bpp::Tree* initTreeDistanceNoOpt(bpp::AgglomerativeDistanceMethod *distMethod,bpp::TransitionModel *dmodel,
                                                            bpp::DiscreteDistribution *rDist,VectorSiteContainer *sitesDistMethod);

    };


} //end of namespace bpp;

namespace CastorApplicationUtils {

    void printParameters(ParameterList &parameters,bpp::SubstitutionModel *smodel);

}

#endif //CASTOR_JATIAPPLICATION_HPP
