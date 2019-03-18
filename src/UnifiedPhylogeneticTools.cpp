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
 * @file UnifiedPhylogeneticTools.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 28 06 2018
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
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <glog/logging.h>
#include "UnifiedPhylogeneticTools.hpp"
#include "PIP.hpp"

bpp::SubstitutionModel *UnifiedPhylogeneticTools::getUnifiedSubstitutionModel(const Alphabet *alphabet,
                                                                              const GeneticCode *gCode,
                                                                              const SiteContainer *data,
                                                                              std::map<std::string, std::string> &params,
                                                                              const std::string &suffix,
                                                                              bool suffixIsOptional,
                                                                              bool verbose,
                                                                              int warn) throw (bpp::Exception){
    /*
    // Get user arguments
    std::string PAR_Alphabet = ApplicationTools::getStringParameter("alphabet", params, "DNA", "", true, true);
    bool PAR_alignment = ApplicationTools::getBooleanParameter("alignment", params, false);
    std::string PAR_model_substitution = ApplicationTools::getStringParameter("model", params, "JC69", "", true, true);

    // Split model string description and test if PIP is required
    std::string modelStringName;
    std::map<std::string, std::string> modelMap;
    KeyvalTools::parseProcedure(PAR_model_substitution, modelStringName, modelMap);
    bool PAR_model_indels = modelStringName == "PIP";

    bpp::SubstitutionModel *smodel = nullptr;
    double lambda;
    double mu;
    bool estimatePIPparameters;

    // Instantiate a substitution model and extend it with PIP
    if (PAR_model_indels) {

        bool computeFrequenciesFromData = false;

        // If frequencies are estimated from the data, but there is no alignment, then flag it.
        if (PAR_alignment) {
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
                            computeFrequenciesFromData = true;
                        }
                    }
                    baseModel += ",";
                }
                baseModel.pop_back();
                baseModel += ")";
                modelMap["model"] = baseModel;
            }

        }

        // Instantiation of the canonical substitution model
        if (PAR_Alphabet.find("Codon") != std::string::npos || PAR_Alphabet.find("Protein") != std::string::npos) {
            smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabetNoGaps, gCode.get(), sites, modelMap, "", true, false, 0);
        } else {
            smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, modelMap, "", true, false, 0);
        }

        // If PIP, then check if lambda/mu initial values are estimated from the data
        estimatePIPparameters = !(modelMap.find("initFreqs=observed") == modelMap.end());

        if (estimatePIPparameters) {

            if (PAR_alignment) {
                lambda = bpp::estimateLambdaFromData(tree, sequences, PAR_proportion);
                mu = bpp::estimateMuFromData(tree, PAR_proportion);
                DLOG(INFO) << "[PIP model] Estimated PIP parameters from data using input sequences (lambda=" << lambda << ",mu=" << mu << ")";
            } else {
                lambda = bpp::estimateLambdaFromData(tree, sites);
                mu = bpp::estimateMuFromData(tree, sites);
                DLOG(INFO) << "[PIP model] Estimated PIP parameters from data using input alignment (lambda=" << lambda << ",mu=" << mu << ")";
            }

        } else {
            lambda = (modelMap.find("lambda") == modelMap.end()) ? 0.1 : std::stod(modelMap["lambda"]);
            mu = (modelMap.find("mu") == modelMap.end()) ? 0.2 : std::stod(modelMap["mu"]);
        }

        // Instatiate the corrisponding PIP model given the alphabet
        if (PAR_Alphabet.find("DNA") != std::string::npos && PAR_Alphabet.find("Codon") == std::string::npos) {
            smodel = new PIP_Nuc(dynamic_cast<NucleicAlphabet *>(alphabet), smodel, *sequences, lambda, mu, computeFrequenciesFromData);
        } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
            smodel = new PIP_AA(dynamic_cast<ProteicAlphabet *>(alphabet), smodel, *sequences, lambda, mu, computeFrequenciesFromData);
        } else if (PAR_Alphabet.find("Codon") != std::string::npos) {
            smodel = new PIP_Codon(dynamic_cast<CodonAlphabet_Extended *>(alphabet), gCode.get(), smodel, *sequences, lambda, mu, computeFrequenciesFromData);
            ApplicationTools::displayWarning("Codon models are experimental in the current version... use with caution!");
           DLOG(WARNING) << "CODONS activated byt the program is not fully tested under these settings!";
        }

    } else {
        // if the alphabet is not extended, then the gap character is not supported
        if (!PAR_alignment) bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        smodel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, jatiapp.getParams(), "", true, false, 0);
    }
 return smodel;

*/
   return nullptr;
}
