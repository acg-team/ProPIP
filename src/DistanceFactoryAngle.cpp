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
 * @file DistanceFactoryAngle.cpp
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

#include "DistanceFactory.hpp"
#include "DistanceFactoryAngle.hpp"

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>

using namespace DistanceFactoryPrographMSA;

DistanceFactoryPrographMSA::DistanceMatrix DistanceFactoryAngle::computePwDistances(bpp::SequenceContainer *sequences,
                                                                                    int ALPHABET_DIM,
                                                                                    int K,
                                                                                    bool mldist_flag,
                                                                                    bool mldist_gap_flag,
                                                                                    double cutoff_dist,
                                                                                    double indel_rate){

    int num_seq = sequences->getNumberOfSequences();

    std::vector<std::string> seq_names = sequences->getSequencesNames();

    DistanceMatrix distances(num_seq);

    int dimension = DistanceFactoryPrographMSAutils::power(ALPHABET_DIM,K);
    this->CountMatrix_= CountMatrix::Zero(num_seq,dimension);

    Eigen::Matrix<int, Eigen::Dynamic, 1> seq_len(num_seq);

    for(unsigned int i=0; i < num_seq; ++i) {

        const bpp::Sequence &seq = sequences->getSequence(seq_names.at(i));

        seq_len(i) = seq.size();

        Eigen::Array<int,Eigen::Dynamic,1> chars;
        chars.resize(K);
        chars.fill(-1);

        for(unsigned int j=0; j<seq.size(); ++j) {
            for(unsigned int k=1; k<K; ++k) {
                chars[k-1] = chars[k];
            }

            chars[K-1] = seq.getValue(j);

            if(chars[K-1] < 0 || chars[K-1] >= ALPHABET_DIM) {
                chars[K - 1] = -1;
            }

            int index = 0;
            for(unsigned int k=0; k<K; ++k) {
                if(chars[k] == -1) {
                    index = -1;
                    break;
                }

                index *= ALPHABET_DIM;
                index += chars[k];
            }

            if(index != -1) {
                this->CountMatrix_(i,index) += 1;
            }
        }
    }

    CountMatrixSVD counts2 = this->CountMatrix_.template cast<double>().transpose();

    distances.distances = counts2.colwise().norm().array().inverse().matrix().asDiagonal() * counts2.transpose() * counts2 * counts2.colwise().norm().array().inverse().matrix().asDiagonal();
    distances.distances = -1.0 * ((distances.distances.array().square() + 0.4) / 1.4).log();
    if(!(mldist_flag) && !(mldist_gap_flag)) {
        distances.distances = distances.distances.array().exp();
        distances.distances = - 0.5 * (5.0*distances.distances.array() - (45.0*distances.distances.array().square() - 20.0*distances.distances.array()).sqrt()) * distances.distances.array().inverse();
    }

    distances.variances.setConstant(1);
    distances.variances = distances.variances * seq_len.cast<double>().asDiagonal();
    distances.variances = ((distances.variances + distances.variances.transpose())/2).eval().array().inverse();
    distances.variances.array() *= distances.distances.array();

    //adjust minimum variances
    distances.variances.array() = distances.variances.array().max(1e-5 * Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>::Ones(distances.variances.rows(),distances.variances.cols()));

    return distances;
}

namespace DistanceFactoryPrographMSAutils {

    int power(int X, int K) {
        if(K==0){
            return 1;
        }else{
            return X*power(X,K-1);
        }
    }
}
