/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of miniJATI
 *
 * miniJATI is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * miniJATI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with miniJATI. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file DistanceFactoryAngle.hpp
 * @author Adam Szalkowski
 * @author Massimo Maiolo
 * @date 06 08 2020
 * @version 1.0.7
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

#ifndef CASTOR_DISTANCEFACTORYANGLE_HPP
#define CASTOR_DISTANCEFACTORYANGLE_HPP

#include "DistanceFactory.hpp"

namespace DistanceFactoryPrographMSAutils {

    int power(int X, int K);
}

namespace DistanceFactoryPrographMSA {

    class DistanceFactoryAngle : public DistanceFactory{
    private:

        typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> CountMatrix;
        typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> CountMatrixSVD;

        typedef Eigen::Matrix<double, Eigen::Dynamic, 1> SingularValues;

        int ALPHABET_DIM_;
        CountMatrix CountMatrix_;
        CountMatrixSVD CountMatrixSVD_;
        SingularValues SingularValues_;

    public:
        DistanceFactoryAngle(int ALPHABET_DIM,int K){
            this->ALPHABET_DIM_=ALPHABET_DIM;
        };

        virtual ~DistanceFactoryAngle() {};

        virtual DistanceMatrix computePwDistances(bpp::SequenceContainer *sequences,
                int ALPHABET_DIM,
                int K,
                bool mldist_flag,
                bool mldist_gap_flag,
                double cutoff_dist,
                double indel_rate);


    };

}

#endif //CASTOR_DISTANCEFACTORYANGLE_HPP
