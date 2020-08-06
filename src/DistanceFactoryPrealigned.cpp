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
 * @file DistanceFactoryPrealigned.cpp
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

#include <iostream>
#include <fstream>

#include "DistanceFactory.hpp"
#include "DistanceFactoryPrealigned.hpp"

using namespace DistanceFactoryPrographMSA;

DistanceFactoryPrographMSA::DistanceMatrix DistanceFactoryPrealigned::computePwDistances(bpp::SequenceContainer *sequences,
                                                                                         int ALPHABET_DIM,
                                                                                         int K,
                                                                                         bool mldist_flag,
                                                                                         bool mldist_gap_flag,
                                                                                         double cutoff_dist,
                                                                                         double indel_rate){
    int num_seq = sequences->getNumberOfSequences();

    std::vector<std::string> seq_names = sequences->getSequencesNames();

    DistanceMatrix distances(num_seq);

    Eigen::Matrix<int, Eigen::Dynamic, 1> seq_len(num_seq);

    for(unsigned int i=0; i < num_seq; ++i) {

        const bpp::Sequence &seq1 = sequences->getSequence(seq_names.at(i));

        seq_len(i) = seq1.size();

        for(unsigned int j=i+1; j < num_seq; ++j) {

            const bpp::Sequence &seq2 = sequences->getSequence(seq_names.at(j));

            CountMatrix counts = CountMatrix::Zero(this->ALPHABET_DIM_,this->ALPHABET_DIM_);

            unsigned int gaps = 0;

            bool gap_opened1 = false;
            bool gap_opened2 = false;

            for(unsigned int k=0; k<seq1.size(); ++k) {

                if(!(seq1[k]=='-') && !(seq2[k]=='-')) {
                    int c1 = seq1.getValue(k);
                    int c2 = seq2.getValue(k);
                    if(c1 >= 0 && c1 < 20 && c2 >= 0 && c2 < 20)
                        ++counts(c1,c2);
                    gap_opened1 = false;
                    gap_opened2 = false;

                } else if (seq1[k]=='-' && seq2[k]=='-') {
                    // skip
                } else if(!(seq1[k]=='-') && !gap_opened1) {
                    ++gaps;
                    gap_opened1 = true;
                    gap_opened2 = false;
                } else if(!(seq2[k]=='-') && !gap_opened2) {
                    ++gaps;
                    gap_opened1 = false;
                    gap_opened2 = true;
                }
            }

            distvar_t distvar = this->computeDistance(counts,
                                                      gaps,
                                                      (seq1.size()+seq2.size())/2.0,
                                                      mldist_flag,
                                                      mldist_gap_flag,
                                                      cutoff_dist,
                                                      indel_rate);

            distances.distances(i, j) = distances.distances(j, i) = distvar.dist;
            distances.variances(i, j) = distances.variances(j, i) = distvar.var;
        }
    }

    return distances;
}