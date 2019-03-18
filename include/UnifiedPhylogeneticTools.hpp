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
 * @file UnifiedPhylogeneticTools.hpp
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
#ifndef CASTOR_UNIFIEDPHYLOGENETICTOOLS_HPP
#define CASTOR_UNIFIEDPHYLOGENETICTOOLS_HPP

#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/Container/SiteContainer.h>

using namespace bpp;

namespace UnifiedPhylogeneticTools {

    bpp::SubstitutionModel *getUnifiedSubstitutionModel(const Alphabet* alphabet,
                                                    const GeneticCode* gCode,
                                                    const SiteContainer* data,
                                                    std::map<std::string, std::string>& params,
                                                    const std::string& suffix = "",
                                                    bool suffixIsOptional = true,
                                                    bool verbose = true,
                                                    int warn = 1) throw (bpp::Exception);

}


#endif //CASTOR_UNIFIEDPHYLOGENETICTOOLS_HPP
