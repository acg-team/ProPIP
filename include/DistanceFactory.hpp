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
 * @file DistanceFactory.hpp
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
#ifndef CASTOR_DISTANCEFACTORY_HPP
#define CASTOR_DISTANCEFACTORY_HPP

#include <Bpp/Phyl/Model/SubstitutionModel.h>

#include <eigen3/Eigen/Core>
#include <map>
#include <vector>
#include <string>

#define DISTCUTOFF 0.000001

namespace DistanceFactoryPrographMSA {

    struct DistanceMatrix {
        DistanceMatrix(int dim) {
            distances = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(dim, dim);
            variances = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(dim, dim);
        }

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> distances;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> variances;

        DistanceMatrix reduce(int i) {
            int dim = this->distances.rows();
            DistanceMatrix reduced_dist(dim - 1);

            if (i == 0) {
                reduced_dist.distances = this->distances.bottomRightCorner(dim - 1, dim - 1);
                reduced_dist.variances = this->variances.bottomRightCorner(dim - 1, dim - 1);
            } else if (i == dim - 1) {
                reduced_dist.distances = this->distances.topLeftCorner(dim - 1, dim - 1);
                reduced_dist.variances = this->variances.topLeftCorner(dim - 1, dim - 1);
            } else {
                reduced_dist.distances << this->distances.topLeftCorner(i, i), this->distances.topRightCorner(i,dim - i -1),
                        this->distances.bottomLeftCorner(dim - i - 1, i), this->distances.bottomRightCorner(dim - i - 1,dim - i -1);

                reduced_dist.variances << this->variances.topLeftCorner(i, i), this->variances.topRightCorner(i,dim - i -1),
                        this->variances.bottomLeftCorner(dim - i - 1, i), this->variances.bottomRightCorner(dim - i - 1,dim - i -1);
            }

            return reduced_dist;
        }

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    };

    class DistanceFactory {
    public:
        virtual DistanceMatrix computePwDistances(bpp::SequenceContainer *sequences,
                int ALPHABET_DIM,
                int K,
                bool mldist_flag,
                bool mldist_gap_flag,
                double cutoff_dist,
                double indel_rate) = 0;

        virtual ~DistanceFactory() {};

        static DistanceFactory *getDefault(bpp::SubstitutionModel *smodel = NULL,
                bool prealigned = false,
                int ALPHABET_DIM=20,
                bool nwdist_flag=true);

    };

}

#endif //CASTOR_DISTANCEFACTORY_HPP
