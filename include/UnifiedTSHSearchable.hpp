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
 * @file UnifiedTSHSearchable.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 04 2018
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
#ifndef CASTOR_UNIFIEDTSHSEARCHABLE_HPP
#define CASTOR_UNIFIEDTSHSEARCHABLE_HPP

#include <Bpp/Numeric/Function/AbstractNumericalDerivative.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/FivePointsNumericalDerivative.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>

namespace bpp {

    class UnifiedTSHSearchable {

    protected:

        enum class LKDataClass { sub, empty };

        mutable AbstractHomogeneousTreeLikelihood *lk_;
        mutable AbstractOptimizer *optimiser_;
        std::string optMethodModel_;
        bool optNumericalDerivatives_;

        // Temporary object for storing FV components (unlimited) during tree-search
        // - allocation at the beginning of the tree-search cycle
        // - deallocation at the end of the tree-search cycle
        // - one or more VVVdouble object per each tree-rearrangement (map classVVVdouble (map idxRearrangment, map idxNode, VVVdouble))
        mutable std::map<LKDataClass,std::map<int,std::map<int, VVVdouble>>> testVectorLikelihoodData_;



    public:

        std::string OPTIMIZATION_NEWTON = "newton";
        std::string OPTIMIZATION_GRADIENT = "gradient";
        std::string OPTIMIZATION_BRENT = "Brent";
        std::string OPTIMIZATION_BFGS = "BFGS";

        UnifiedTSHSearchable() : lk_(nullptr), optimiser_(nullptr), optMethodModel_(""), optNumericalDerivatives_(false) {}

        virtual ~UnifiedTSHSearchable() {

            delete optimiser_->getFunction();

        };

        void setOptimiser(AbstractHomogeneousTreeLikelihood *lk,
                                        bool optNumericalDerivatives,
                                        std::map<std::string, std::string> &params,
                                        const string &suffix,
                                        bool suffixIsOptional,
                                        bool verbose,
                                        int warn);

        void fireBranchOptimisation(std::vector<bpp::Node *> extractionNodes);

        std::map<int,std::map<int, VVVdouble>> *getTestLikelihoodData(LKDataClass classTestLKData = LKDataClass::sub) {
            return &testVectorLikelihoodData_[classTestLKData];
        }

        virtual void addTestLikelihoodData(int idxThread){}

        virtual void removeTestLikelihoodData(int idxThread){}

    };
}


#endif //CASTOR_UNIFIEDTSHSEARCHABLE_HPP
