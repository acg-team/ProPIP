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
 * @file DistanceFactoryAlign.cpp
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
#include "DistanceFactoryAlign.hpp"

using namespace DistanceFactoryPrographMSA;

DistanceMatrix DistanceFactoryAlign::computePwDistances(bpp::SequenceContainer *sequences,
                                                        int ALPHABET_DIM,
                                                        int K,
                                                        bool mldist_flag,
                                                        bool mldist_gap_flag,
                                                        double cutoff_dist,
                                                        double indel_rate) {

    int num_seq = sequences->getNumberOfSequences();

    std::vector<std::string> seq_names = sequences->getSequencesNames();

    DistanceMatrix distances(sequences->getNumberOfSequences());
    Eigen::Matrix<int, Eigen::Dynamic, 1> seq_len(num_seq);

    for(unsigned int i = 0;i<num_seq;++i) {

        const bpp::Sequence &seq1 = sequences->getSequence(seq_names.at(i));
        std::string s_seq1 = seq1.toString();

        for(unsigned int j = i + 1;j<num_seq;++j) {

            const bpp::Sequence &seq2 = sequences->getSequence(seq_names.at(j));
            std::string s_seq2 = seq2.toString();

            distvar_t distvar = this->alignPair(s_seq1,
                                                s_seq2,
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

distvar_t DistanceFactoryAlign::alignPair(const std::string &seq1,
                                          const std::string &seq2,
                                          bool mldist_flag,
                                          bool mldist_gap_flag,
                                          double cutoff_dist,
                                          double indel_rate) {

    const int minfty = -10000;
    DynProgMatrix W((int) seq2.size() + 1, (int) seq1.size() + 1);
    DynProgMatrix X((int) seq2.size() + 1, (int) seq1.size() + 1);
    DynProgMatrix Y((int) seq2.size() + 1, (int) seq1.size() + 1);

    Eigen::Matrix<int, Eigen::Dynamic, 1> s1((int) seq1.size() + 1);
    Eigen::Matrix<int, Eigen::Dynamic, 1> s2((int) seq2.size() + 1);

    W(0, 0) = 0;

    for (unsigned int x = 1; x < seq1.size() + 1; ++x) {
        s1(x) = seq1[x - 1];
        if (s1(x) < 0) s1(x) = 20;
        X(0, x) = W(0, x) = this->gap_open + (x - 1) * this->gap_extend;
        Y(0, x) = minfty;
    }

    for (unsigned int y = 1; y < seq2.size() + 1; ++y) {
        s2(y) = seq2[y - 1];
        if (s2(y) < 0) s2(y) = 20;
        Y(y, 0) = W(y, 0) = this->gap_open + (y - 1) * this->gap_extend;
        X(y, 0) = minfty;
    }

    for (unsigned int y = 1; y < seq2.size() + 1; ++y) {
        for (unsigned int x = 1; x < seq1.size() + 1; ++x) {
            W(y, x) = W(y - 1, x - 1) + this->scoring_matrix(s2(y), s1(x));
            X(y, x) = std::max(X(y, x - 1) + this->gap_extend, W(y, x - 1) + this->gap_open);
            Y(y, x) = std::max(Y(y - 1, x) + this->gap_extend, W(y - 1, x) + this->gap_open);
            W(y, x) = std::max(std::max(X(y, x), Y(y, x)), W(y, x));
        }
    }


    CountMatrix counts = CountMatrix::Zero(this->ALPHABET_DIM_,this->ALPHABET_DIM_);

    unsigned int gaps = 0;

    bool gap_opened1 = false;
    bool gap_opened2 = false;

    // backtrack and count identical characters in matches
    for (unsigned int y = seq2.size(), x = seq1.size(); y != 0 && x != 0;) {
        if (W(y, x) == W(y - 1, x - 1) + this->scoring_matrix(s2(y), s1(x))) {
            if (s1(x) < this->ALPHABET_DIM_ && s2(y) < this->ALPHABET_DIM_)
                ++counts(s1(x), s2(y));
            gap_opened1 = false;
            gap_opened2 = false;
            --x;
            --y;

        } else if (W(y, x) == X(y, x)) {
            if (!gap_opened1)
                ++gaps;
            gap_opened1 = true;
            gap_opened2 = false;
            --x;

        } else if (W(y, x) == Y(y, x)) {
            if (!gap_opened2)
                ++gaps;
            gap_opened1 = false;
            gap_opened2 = true;
            --y;

        } else {
            LOG(FATAL) << "error while backtracking: ";
        }
    }

    distvar_t distvar = this->computeDistance(counts,
                                              gaps,
                                              (seq1.length() + seq2.length()) / 2.0,
                                              mldist_flag,
                                              mldist_gap_flag,
                                              cutoff_dist,
                                              indel_rate);

    return distvar;
}

void DistanceFactoryAlign::initMatrix_AA() {
    // BLOSUM64
    this->scoring_matrix << 4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1,
            1, 0, 0, -3, -2, 0, 0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3,
            -1, -1, -1, -2, -2, -2, -2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2,
            0, -1, -3, -4, -3, -1, -1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0,
            -1, -2, -3, -2, -1, -2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2,
            -2, -1, 1, 3, -1, 0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0,
            -2, -3, -2, -3, -1, -2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1,
            -2, -3, -2, 2, -1, -1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2,
            -1, 3, -3, -1, -1, -1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1,
            -2, -3, -2, -1, -1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1,
            1, -2, -1, -1, -1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1,
            -1, -1, -1, -2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4,
            -2, -1, -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2,
            -4, -3, -2, -1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2,
            -1, -1, -1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3,
            -2, -1, 1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2,
            0, 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2, 0,
            0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1, -1,
            -3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2,
            -2, -2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7,
            -1, 0, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1, 0, 0, -1, -2, -1,
            -1;

    this->gap_open = -10;
    this->gap_extend = -2;

}

void DistanceFactoryAlign::initMatrix_Codon() {
    // BLOSUM64 expanded to codon alphabet

    this->scoring_matrix << 6, 6, 0, 0, -2, -2, -2, -2, 3, 3, -2, -2, 1, 0, 0,
            0, 0, -4, -4, -4, -4, -1, -1, -3, -3, -3, -3, -3, -3, 0, 0, 0, 0, -2, -2,
            -2, -2, -3, -3, -3, -3, -2, -2, -3, -3, -1, -1, -1, -1, -2, -2, -2, -2, -3,
            -3, -3, -3, -3, -3, -3, -3, -1, 6, 6, 0, 0, -2, -2, -2, -2, 3, 3, -2, -2, 1,
            0, 0, 0, 0, -4, -4, -4, -4, -1, -1, -3, -3, -3, -3, -3, -3, 0, 0, 0, 0, -2,
            -2, -2, -2, -3, -3, -3, -3, -2, -2, -3, -3, -1, -1, -1, -1, -2, -2, -2, -2,
            -3, -3, -3, -3, -3, -3, -3, -3, -1, 0, 0, 4, 4, -2, -2, -2, -2, -1, -1, -1,
            -1, -2, 4, 4, 4, 4, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, 2, 2, 2,
            2, -1, -1, -1, -1, -3, -3, -2, -2, -2, -2, -2, -2, 1, 1, 1, 1, -1, -1, -1,
            -1, -4, -4, -3, -3, -4, -4, -4, -4, -1, 0, 0, 4, 4, -2, -2, -2, -2, -1, -1,
            -1, -1, -2, 4, 4, 4, 4, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, 2,
            2, 2, 2, -1, -1, -1, -1, -3, -3, -2, -2, -2, -2, -2, -2, 1, 1, 1, 1, -1, -1,
            -1, -1, -4, -4, -3, -3, -4, -4, -4, -4, -1, -2, -2, -2, -2, 4, 4, 4, 4, -2,
            -2, -1, -1, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1,
            -1, -2, -2, -2, -1, 1, 1, 1, 1, 1, 1, 0, 0, 4, 4, -1, -1, -2, -2, -2, -2, 1,
            1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, -2, -2, 4, 4, 4, 4, -2, -2, -1,
            -1, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -2,
            -2, -2, -1, 1, 1, 1, 1, 1, 1, 0, 0, 4, 4, -1, -1, -2, -2, -2, -2, 1, 1, 1,
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, -2, -2, 4, 4, 4, 4, -2, -2, -1, -1,
            -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -2, -2,
            -2, -1, 1, 1, 1, 1, 1, 1, 0, 0, 4, 4, -1, -1, -2, -2, -2, -2, 1, 1, 1, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, -2, -2, -2, -2, 4, 4, 4, 4, -2, -2, -1, -1, -3, -2,
            -2, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -2, -2, -2, -1, 1,
            1, 1, 1, 1, 1, 0, 0, 4, 4, -1, -1, -2, -2, -2, -2, 1, 1, 1, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 3, 3, -1, -1, -2, -2, -2, -2, 7, 7, -2, -2, 2, -1, -1, -1,
            -1, -3, -3, -3, -3, 2, 2, -1, -1, -2, -2, -2, -2, -1, -1, -1, -1, -2, -2,
            -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -2, -2, -2, -2, -3,
            -3, -2, -2, -3, -3, -3, -3, -1, 3, 3, -1, -1, -2, -2, -2, -2, 7, 7, -2, -2,
            2, -1, -1, -1, -1, -3, -3, -3, -3, 2, 2, -1, -1, -2, -2, -2, -2, -1, -1, -1,
            -1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -2, -2,
            -2, -2, -3, -3, -2, -2, -3, -3, -3, -3, -1, -2, -2, -1, -1, -1, -1, -1, -1,
            -2, -2, 9, 9, -2, -1, -1, -1, -1, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
            -3, -3, -1, -1, -1, -1, -1, -1, -1, -1, -3, -3, -3, -3, -1, -1, -3, -3, -1,
            -1, -1, -1, 0, 0, 0, 0, -3, -3, -4, -4, -3, -3, -3, -3, -2, -2, -2, -1, -1,
            -1, -1, -1, -1, -2, -2, 9, 9, -2, -1, -1, -1, -1, -3, -3, -3, -3, -3, -3,
            -3, -3, -3, -3, -3, -3, -1, -1, -1, -1, -1, -1, -1, -1, -3, -3, -3, -3, -1,
            -1, -3, -3, -1, -1, -1, -1, 0, 0, 0, 0, -3, -3, -4, -4, -3, -3, -3, -3, -2,
            1, 1, -2, -2, -3, -3, -3, -3, 2, 2, -2, -2, 11, -2, -2, -2, -2, -4, -4, -4,
            -4, -2, -2, -2, -2, -3, -3, -3, -3, -3, -3, -3, -1, -2, -2, -2, -2, -4, -4,
            -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -4, -4, -3, -3, -2,
            -2, -2, -2, -2, 0, 0, 4, 4, -2, -2, -2, -2, -1, -1, -1, -1, -2, 4, 4, 4, 4,
            -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, 2, 2, 2, 2, -1, -1, -1, -1,
            -3, -3, -2, -2, -2, -2, -2, -2, 1, 1, 1, 1, -1, -1, -1, -1, -4, -4, -3, -3,
            -4, -4, -4, -4, -1, 0, 0, 4, 4, -2, -2, -2, -2, -1, -1, -1, -1, -2, 4, 4, 4,
            4, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, 2, 2, 2, 2, -1, -1, -1,
            -1, -3, -3, -2, -2, -2, -2, -2, -2, 1, 1, 1, 1, -1, -1, -1, -1, -4, -4, -3,
            -3, -4, -4, -4, -4, -1, 0, 0, 4, 4, -2, -2, -2, -2, -1, -1, -1, -1, -2, 4,
            4, 4, 4, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, 2, 2, 2, 2, -1, -1,
            -1, -1, -3, -3, -2, -2, -2, -2, -2, -2, 1, 1, 1, 1, -1, -1, -1, -1, -4, -4,
            -3, -3, -4, -4, -4, -4, -1, 0, 0, 4, 4, -2, -2, -2, -2, -1, -1, -1, -1, -2,
            4, 4, 4, 4, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, 2, 2, 2, 2, -1,
            -1, -1, -1, -3, -3, -2, -2, -2, -2, -2, -2, 1, 1, 1, 1, -1, -1, -1, -1, -4,
            -4, -3, -3, -4, -4, -4, -4, -1, -4, -4, -3, -3, -1, -1, -1, -1, -3, -3, -3,
            -3, -4, -3, -3, -3, -3, 7, 7, 7, 7, -2, -2, -1, -1, -2, -2, -2, -2, -3, -3,
            -3, -2, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -2, -2, -2, -2, -2, -1,
            -1, -1, -1, -1, -1, -1, -1, -2, -2, -2, -2, -2, -4, -4, -3, -3, -1, -1, -1,
            -1, -3, -3, -3, -3, -4, -3, -3, -3, -3, 7, 7, 7, 7, -2, -2, -1, -1, -2, -2,
            -2, -2, -3, -3, -3, -2, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -2, -2,
            -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -2, -2, -2, -2, -2, -4, -4, -3,
            -3, -1, -1, -1, -1, -3, -3, -3, -3, -4, -3, -3, -3, -3, 7, 7, 7, 7, -2, -2,
            -1, -1, -2, -2, -2, -2, -3, -3, -3, -2, -1, -1, -1, -1, -2, -2, -1, -1, -1,
            -1, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -2, -2, -2, -2,
            -2, -4, -4, -3, -3, -1, -1, -1, -1, -3, -3, -3, -3, -4, -3, -3, -3, -3, 7,
            7, 7, 7, -2, -2, -1, -1, -2, -2, -2, -2, -3, -3, -3, -2, -1, -1, -1, -1, -2,
            -2, -1, -1, -1, -1, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1,
            -2, -2, -2, -2, -2, -1, -1, -3, -3, -1, -1, -1, -1, 2, 2, -3, -3, -2, -3,
            -3, -3, -3, -2, -2, -2, -2, 8, 8, 0, 0, 0, 0, 0, 0, -3, -3, -3, -2, -2, -2,
            -2, -2, 1, 1, -1, -1, -1, -1, 0, 0, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1,
            0, 0, -2, -2, -2, -2, -1, -1, -1, -3, -3, -1, -1, -1, -1, 2, 2, -3, -3, -2,
            -3, -3, -3, -3, -2, -2, -2, -2, 8, 8, 0, 0, 0, 0, 0, 0, -3, -3, -3, -2, -2,
            -2, -2, -2, 1, 1, -1, -1, -1, -1, 0, 0, -3, -3, -3, -3, -2, -2, -2, -2, -1,
            -1, 0, 0, -2, -2, -2, -2, -1, -3, -3, -2, -2, 0, 0, 0, 0, -1, -1, -3, -3,
            -2, -2, -2, -2, -2, -1, -1, -1, -1, 0, 0, 5, 5, 1, 1, 1, 1, -3, -3, -3, 0,
            -1, -1, -1, -1, 0, 0, 1, 1, 0, 0, 1, 1, -2, -2, -2, -2, -1, -1, -1, -1, 0,
            0, 2, 2, -2, -2, -2, -2, -1, -3, -3, -2, -2, 0, 0, 0, 0, -1, -1, -3, -3, -2,
            -2, -2, -2, -2, -1, -1, -1, -1, 0, 0, 5, 5, 1, 1, 1, 1, -3, -3, -3, 0, -1,
            -1, -1, -1, 0, 0, 1, 1, 0, 0, 1, 1, -2, -2, -2, -2, -1, -1, -1, -1, 0, 0, 2,
            2, -2, -2, -2, -2, -1, -3, -3, -2, -2, -1, -1, -1, -1, -2, -2, -3, -3, -3,
            -2, -2, -2, -2, -2, -2, -2, -2, 0, 0, 1, 1, 5, 5, 5, 5, -3, -3, -3, -1, -1,
            -1, -1, -1, 0, 0, 2, 2, -1, -1, 5, 5, -3, -3, -3, -3, -1, -1, -1, -1, -2,
            -2, 0, 0, -2, -2, -2, -2, -1, -3, -3, -2, -2, -1, -1, -1, -1, -2, -2, -3,
            -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, 0, 0, 1, 1, 5, 5, 5, 5, -3, -3, -3,
            -1, -1, -1, -1, -1, 0, 0, 2, 2, -1, -1, 5, 5, -3, -3, -3, -3, -1, -1, -1,
            -1, -2, -2, 0, 0, -2, -2, -2, -2, -1, -3, -3, -2, -2, -1, -1, -1, -1, -2,
            -2, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, 0, 0, 1, 1, 5, 5, 5, 5, -3,
            -3, -3, -1, -1, -1, -1, -1, 0, 0, 2, 2, -1, -1, 5, 5, -3, -3, -3, -3, -1,
            -1, -1, -1, -2, -2, 0, 0, -2, -2, -2, -2, -1, -3, -3, -2, -2, -1, -1, -1,
            -1, -2, -2, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, 0, 0, 1, 1, 5, 5, 5,
            5, -3, -3, -3, -1, -1, -1, -1, -1, 0, 0, 2, 2, -1, -1, 5, 5, -3, -3, -3, -3,
            -1, -1, -1, -1, -2, -2, 0, 0, -2, -2, -2, -2, -1, 0, 0, 2, 2, -2, -2, -2,
            -2, -1, -1, -1, -1, -3, 2, 2, 2, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
            -3, -3, 4, 4, 4, 1, -1, -1, -1, -1, -3, -3, -3, -3, -2, -2, -3, -3, 3, 3, 3,
            3, -1, -1, -1, -1, -3, -3, -3, -3, -4, -4, -4, -4, -1, 0, 0, 2, 2, -2, -2,
            -2, -2, -1, -1, -1, -1, -3, 2, 2, 2, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3,
            -3, -3, -3, 4, 4, 4, 1, -1, -1, -1, -1, -3, -3, -3, -3, -2, -2, -3, -3, 3,
            3, 3, 3, -1, -1, -1, -1, -3, -3, -3, -3, -4, -4, -4, -4, -1, 0, 0, 2, 2, -2,
            -2, -2, -2, -1, -1, -1, -1, -3, 2, 2, 2, 2, -3, -3, -3, -3, -3, -3, -3, -3,
            -3, -3, -3, -3, 4, 4, 4, 1, -1, -1, -1, -1, -3, -3, -3, -3, -2, -2, -3, -3,
            3, 3, 3, 3, -1, -1, -1, -1, -3, -3, -3, -3, -4, -4, -4, -4, -1, 0, 0, 2, 2,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, 2, 2, 2, 2, -2, -2, -2, -2, -2, -2, 0,
            0, -1, -1, -1, -1, 1, 1, 1, 5, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1,
            -1, 1, 1, 1, 1, -1, -1, -1, -1, -3, -3, -2, -2, -3, -3, -3, -3, -1, -2, -2,
            -1, -1, 1, 1, 1, 1, -2, -2, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -2,
            -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 5, 5, 5, 5, 0, 0, -1, -1, 1, 1,
            -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -2, -2, -2, -2, 0, -2, -2,
            -1, -1, 1, 1, 1, 1, -2, -2, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -2,
            -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 5, 5, 5, 5, 0, 0, -1, -1, 1, 1,
            -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -2, -2, -2, -2, 0, -2, -2,
            -1, -1, 1, 1, 1, 1, -2, -2, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -2,
            -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 5, 5, 5, 5, 0, 0, -1, -1, 1, 1,
            -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -2, -2, -2, -2, 0, -2, -2,
            -1, -1, 1, 1, 1, 1, -2, -2, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -2,
            -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 5, 5, 5, 5, 0, 0, -1, -1, 1, 1,
            -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -2, -2, -2, -2, 0, -3, -3,
            -3, -3, 1, 1, 1, 1, -2, -2, -3, -3, -4, -3, -3, -3, -3, -2, -2, -2, -2, 1,
            1, 0, 0, 0, 0, 0, 0, -3, -3, -3, -2, 0, 0, 0, 0, 6, 6, 0, 0, 1, 1, 0, 0, -3,
            -3, -3, -3, -2, -2, -2, -2, 1, 1, 0, 0, 0, 0, 0, 0, -1, -3, -3, -3, -3, 1,
            1, 1, 1, -2, -2, -3, -3, -4, -3, -3, -3, -3, -2, -2, -2, -2, 1, 1, 0, 0, 0,
            0, 0, 0, -3, -3, -3, -2, 0, 0, 0, 0, 6, 6, 0, 0, 1, 1, 0, 0, -3, -3, -3, -3,
            -2, -2, -2, -2, 1, 1, 0, 0, 0, 0, 0, 0, -1, -3, -3, -2, -2, 0, 0, 0, 0, -2,
            -2, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, 1, 1, 2, 2, 2, 2,
            -3, -3, -3, -1, -1, -1, -1, -1, 0, 0, 5, 5, 0, 0, 2, 2, -2, -2, -2, -2, -1,
            -1, -1, -1, -1, -1, 1, 1, -2, -2, -2, -2, -1, -3, -3, -2, -2, 0, 0, 0, 0,
            -2, -2, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, 1, 1, 2, 2, 2,
            2, -3, -3, -3, -1, -1, -1, -1, -1, 0, 0, 5, 5, 0, 0, 2, 2, -2, -2, -2, -2,
            -1, -1, -1, -1, -1, -1, 1, 1, -2, -2, -2, -2, -1, -2, -2, -2, -2, 4, 4, 4,
            4, -2, -2, -1, -1, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1,
            -1, -1, -2, -2, -2, -1, 1, 1, 1, 1, 1, 1, 0, 0, 4, 4, -1, -1, -2, -2, -2,
            -2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, -2, -2, 4, 4, 4, 4, -2,
            -2, -1, -1, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1,
            -1, -2, -2, -2, -1, 1, 1, 1, 1, 1, 1, 0, 0, 4, 4, -1, -1, -2, -2, -2, -2, 1,
            1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, -3, -2, -2, -1, -1, -1, -1, -2, -2,
            -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, 0, 0, 1, 1, 5, 5, 5, 5, -3, -3,
            -3, -1, -1, -1, -1, -1, 0, 0, 2, 2, -1, -1, 5, 5, -3, -3, -3, -3, -1, -1,
            -1, -1, -2, -2, 0, 0, -2, -2, -2, -2, -1, -3, -3, -2, -2, -1, -1, -1, -1,
            -2, -2, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, 0, 0, 1, 1, 5, 5, 5, 5,
            -3, -3, -3, -1, -1, -1, -1, -1, 0, 0, 2, 2, -1, -1, 5, 5, -3, -3, -3, -3,
            -1, -1, -1, -1, -2, -2, 0, 0, -2, -2, -2, -2, -1, -1, -1, 1, 1, -2, -2, -2,
            -2, -1, -1, -1, -1, -3, 1, 1, 1, 1, -2, -2, -2, -2, -3, -3, -2, -2, -3, -3,
            -3, -3, 3, 3, 3, 1, 0, 0, 0, 0, -3, -3, -2, -2, -2, -2, -3, -3, 4, 4, 4, 4,
            0, 0, 0, 0, -3, -3, -2, -2, -3, -3, -3, -3, -1, -1, -1, 1, 1, -2, -2, -2,
            -2, -1, -1, -1, -1, -3, 1, 1, 1, 1, -2, -2, -2, -2, -3, -3, -2, -2, -3, -3,
            -3, -3, 3, 3, 3, 1, 0, 0, 0, 0, -3, -3, -2, -2, -2, -2, -3, -3, 4, 4, 4, 4,
            0, 0, 0, 0, -3, -3, -2, -2, -3, -3, -3, -3, -1, -1, -1, 1, 1, -2, -2, -2,
            -2, -1, -1, -1, -1, -3, 1, 1, 1, 1, -2, -2, -2, -2, -3, -3, -2, -2, -3, -3,
            -3, -3, 3, 3, 3, 1, 0, 0, 0, 0, -3, -3, -2, -2, -2, -2, -3, -3, 4, 4, 4, 4,
            0, 0, 0, 0, -3, -3, -2, -2, -3, -3, -3, -3, -1, -1, -1, 1, 1, -2, -2, -2,
            -2, -1, -1, -1, -1, -3, 1, 1, 1, 1, -2, -2, -2, -2, -3, -3, -2, -2, -3, -3,
            -3, -3, 3, 3, 3, 1, 0, 0, 0, 0, -3, -3, -2, -2, -2, -2, -3, -3, 4, 4, 4, 4,
            0, 0, 0, 0, -3, -3, -2, -2, -3, -3, -3, -3, -1, -2, -2, -1, -1, 1, 1, 1, 1,
            -2, -2, 0, 0, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, -2, -2, -1, -1, 1, 1, -1, -1, 0, 0, 0,
            0, 4, 4, 4, 4, -2, -2, -1, -1, 0, 0, 0, 0, 0, -2, -2, -1, -1, 1, 1, 1, 1,
            -2, -2, 0, 0, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, -2, -2, -1, -1, 1, 1, -1, -1, 0, 0, 0,
            0, 4, 4, 4, 4, -2, -2, -1, -1, 0, 0, 0, 0, 0, -2, -2, -1, -1, 1, 1, 1, 1,
            -2, -2, 0, 0, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, -2, -2, -1, -1, 1, 1, -1, -1, 0, 0, 0,
            0, 4, 4, 4, 4, -2, -2, -1, -1, 0, 0, 0, 0, 0, -2, -2, -1, -1, 1, 1, 1, 1,
            -2, -2, 0, 0, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, -2, -2, -1, -1, 1, 1, -1, -1, 0, 0, 0,
            0, 4, 4, 4, 4, -2, -2, -1, -1, 0, 0, 0, 0, 0, -3, -3, -4, -4, 0, 0, 0, 0,
            -3, -3, -3, -3, -4, -4, -4, -4, -4, -1, -1, -1, -1, -1, -1, 0, 0, -2, -2,
            -2, -2, -3, -3, -3, -3, -1, -1, -1, -1, 1, 1, -1, -1, 0, 0, -2, -2, -3, -3,
            -3, -3, -2, -2, -2, -2, 6, 6, 2, 2, -1, -1, -1, -1, -1, -3, -3, -4, -4, 0,
            0, 0, 0, -3, -3, -3, -3, -4, -4, -4, -4, -4, -1, -1, -1, -1, -1, -1, 0, 0,
            -2, -2, -2, -2, -3, -3, -3, -3, -1, -1, -1, -1, 1, 1, -1, -1, 0, 0, -2, -2,
            -3, -3, -3, -3, -2, -2, -2, -2, 6, 6, 2, 2, -1, -1, -1, -1, -1, -3, -3, -3,
            -3, 0, 0, 0, 0, -2, -2, -4, -4, -3, -3, -3, -3, -3, -1, -1, -1, -1, 0, 0, 2,
            2, 0, 0, 0, 0, -3, -3, -3, -2, -1, -1, -1, -1, 0, 0, 1, 1, 0, 0, 0, 0, -2,
            -2, -2, -2, -1, -1, -1, -1, 2, 2, 5, 5, -2, -2, -2, -2, -1, -3, -3, -3, -3,
            0, 0, 0, 0, -2, -2, -4, -4, -3, -3, -3, -3, -3, -1, -1, -1, -1, 0, 0, 2, 2,
            0, 0, 0, 0, -3, -3, -3, -2, -1, -1, -1, -1, 0, 0, 1, 1, 0, 0, 0, 0, -2, -2,
            -2, -2, -1, -1, -1, -1, 2, 2, 5, 5, -2, -2, -2, -2, -1, -3, -3, -4, -4, 0,
            0, 0, 0, -3, -3, -3, -3, -2, -4, -4, -4, -4, -2, -2, -2, -2, -2, -2, -2, -2,
            -2, -2, -2, -2, -4, -4, -4, -3, -2, -2, -2, -2, 0, 0, -2, -2, 0, 0, -2, -2,
            -3, -3, -3, -3, 0, 0, 0, 0, -1, -1, -2, -2, 6, 6, 6, 6, -1, -3, -3, -4, -4,
            0, 0, 0, 0, -3, -3, -3, -3, -2, -4, -4, -4, -4, -2, -2, -2, -2, -2, -2, -2,
            -2, -2, -2, -2, -2, -4, -4, -4, -3, -2, -2, -2, -2, 0, 0, -2, -2, 0, 0, -2,
            -2, -3, -3, -3, -3, 0, 0, 0, 0, -1, -1, -2, -2, 6, 6, 6, 6, -1, -3, -3, -4,
            -4, 0, 0, 0, 0, -3, -3, -3, -3, -2, -4, -4, -4, -4, -2, -2, -2, -2, -2, -2,
            -2, -2, -2, -2, -2, -2, -4, -4, -4, -3, -2, -2, -2, -2, 0, 0, -2, -2, 0, 0,
            -2, -2, -3, -3, -3, -3, 0, 0, 0, 0, -1, -1, -2, -2, 6, 6, 6, 6, -1, -3, -3,
            -4, -4, 0, 0, 0, 0, -3, -3, -3, -3, -2, -4, -4, -4, -4, -2, -2, -2, -2, -2,
            -2, -2, -2, -2, -2, -2, -2, -4, -4, -4, -3, -2, -2, -2, -2, 0, 0, -2, -2, 0,
            0, -2, -2, -3, -3, -3, -3, 0, 0, 0, 0, -1, -1, -2, -2, 6, 6, 6, 6, -1, -1,
            -1, -1, -1, 0, 0, 0, 0, -1, -1, -2, -2, -2, -1, -1, -1, -1, -2, -2, -2, -2,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1,
            0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1,
            -1;

    this->gap_open = -10;
    this->gap_extend = -2;
}

void DistanceFactoryAlign::initMatrix_DNA() {
    // FIXME only very simple model for DNA alignment

    this->scoring_matrix << 1, -1, -2, -2, 0, -1, 1, -2, -2, 0, -2, -2, 1, -1, 0,
            -2, -2, -1, 1, 0, 0, 0, 0, 0, 0;

    this->gap_open = -5;
    this->gap_extend = -2;
}