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
 * @file DistanceFactoryML.hpp
 * @author Adam Szalkowski
 * @author Massimo Maiolo
 * @date 06 08 2020
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

#ifndef CASTOR_DISTANCEFACTORYML_HPP
#define CASTOR_DISTANCEFACTORYML_HPP

<<<<<<< HEAD
#include "DistanceFactory.hpp"
#include <eigen3/Eigen/Dense>

#define DIST_MAX_AA 2.2
#define VAR_MAX_AA 1e3
#define VAR_MIN_AA 1e-5
#define DIST_MAX_Codon 5.2
#define VAR_MAX_Codon 5e3
#define VAR_MIN_Codon 1e-5
#define DIST_MAX_DNA 2.2
#define VAR_MAX_DNA 1e3
#define VAR_MIN_DNA 1e-5

namespace DistanceFactoryPrographMSA {

    struct distvar_t {
        distvar_t(double dist, double var) {
            this->dist = dist;
            this->var = var;
        }

        double dist;
        double var;
    };

    class DistanceFactoryML : public DistanceFactory {

    public:

        int ALPHABET_DIM_;
        typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> CountMatrix;
        CountMatrix CountMatrix_;

        DistanceFactoryML(bpp::SubstitutionModel *smodel,int ALPHABET_DIM);

        virtual ~DistanceFactoryML();

    private:

        distvar_t computeMLDist(const CountMatrix &counts,
                unsigned int gaps,
                double seqlen,
                double dist0,
                double var0,
                bool mldist_gap_flag,
                double indel_rate);

        double DIST_MAX;
        double VAR_MAX;
        double VAR_MIN;

    protected:

        distvar_t computeDistance(const CountMatrix &counts,
                                  unsigned int gaps,
                                  double seqlen,
                                  bool mldist_flag,
                                  bool mldist_gap_flag,
                                  double cutoff_dist,
                                  double indel_rate);

        bpp::SubstitutionModel *smodel_;

    };

}

=======
>>>>>>> indel_rates_inference
#endif //CASTOR_DISTANCEFACTORYML_HPP
