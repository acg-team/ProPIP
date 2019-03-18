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
 * @file PIPmsa.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 07 08 2018
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

/*
* From SeqLib:
*/
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <glog/logging.h>
#include <string>

#include "PIPmsa.hpp"
#include "progressivePIP.hpp"

using namespace bpp;

void PIPmsa::_setSeqNameLeaf(std::string &seqName) {

    // leaf has only one sequence
    seqNames_.push_back(seqName);

}

void PIPmsa::_setMSAleaf(const bpp::Sequence *sequence) {

    // convert from Sequence to string
    std::string sequenceString = sequence->toString();

    /* convert a string into a vector of single char strings */
    msa_.resize(sequenceString.size());
    for (int i = 0; i < sequenceString.size(); i++) {

        const char ch = sequenceString.at(i);

        // check sequence content: if 'X', ' ' or '-' exit with error
        if (ch == 'X' || ch == ' ' || ch == '-') {
            LOG(FATAL) << "\nERROR sequence contains 'X' or ' ' or '-'";
        }

        // convert character to string column vector
        MSAcolumn_t msa_col(1, ch);

        // assign MSA column to MSA
        msa_.at(i) = msa_col;
    }

}

void PIPmsa::_setSeqNameNode(std::vector<std::string> &seqNamesL,
                             std::vector<std::string> &seqNamesR) {

    // add first all the sequences name of the left subtree
    for (int i = 0; i < seqNamesL.size(); i++) {
        seqNames_.push_back(seqNamesL.at(i));
    }

    // then add all the sequences name of the right subtree
    for (int i = 0; i < seqNamesR.size(); i++) {
        seqNames_.push_back(seqNamesR.at(i));
    }

}

void PIPmsa::_setFVleaf(int numCatg,
                        const bpp::Alphabet *alphabet) {

    // get the number of compressed sites
    int lenComprSeqs = rev_map_compressed_seqs_.size();

    // get the size of the alphabet (extended)
    int lenAlphabet = alphabet->getSize();

    // resize fv data([site][catg][states])
    fv_data_.resize(lenComprSeqs);
    for (int i = 0; i < lenComprSeqs; i++) {
        fv_data_[i].resize(numCatg);
    }

    int idx;
    // go through all the sites
    for (int i = 0; i < lenComprSeqs; i++) {

        // get the index in the compressed map
        idx = rev_map_compressed_seqs_.at(i);
        MSAcolumn_t s = this->msa_.at(idx);

        // allocate fv column to the ext. alphabet size
        bpp::ColMatrix<double> fv;
        fv.resize(lenAlphabet, 1); // ColMatrix as Nx1 matrix
        bpp::MatrixTools::fill(fv, 0.0); // all zeros

        // check if the sequence contains a "forbidden" char
        if (s[0] == 'X' || s[0] == ' ' || s[0] == '-') {
            LOG(FATAL) << "\nERROR sequence contains either 'X' or ' ' or '-'";
        }

        // get the char position in the alphabet
        idx = alphabet->charToInt(&s[0]);

        // set to 1 the indicator array at the position of the observed char
        fv(idx, 0) = 1.0;

        // assign the indicator array to all the gamma categories
        for (int catg = 0; catg < numCatg; catg++) {
            fv_data_.at(i).at(catg) = fv;
        }

    }

}

void PIPmsa::_setFVsigmaLeaf(int numCatg,
                             const bpp::ColMatrix<double> &pi) {

    // get the size of the compressed sequence
    int lenComprSeqs = fv_data_.size();

    // resize the array ([site][numCatg])
    fv_sigma_.resize(lenComprSeqs);

    double fv0;
    // go through all the sites
    for (int site = 0; site < lenComprSeqs; site++) {

        fv_sigma_.at(site).resize(numCatg);

        // go through all the gamma categories
        for (int catg = 0; catg < numCatg; catg++) {

            // compute fv_sigma = fv dot pi
            fv0 = MatrixBppUtils::dotProd(fv_data_.at(site).at(catg), pi);

            fv_sigma_.at(site).at(catg) = fv0;
        }
    }

}

void PIPmsa::_setFVemptyLeaf(int numCatg,
                             const bpp::Alphabet *alphabet) {

    // get the size of the compressed sequence
    int lenAlphabet = alphabet->getSize();

    // indicator array (all zeros except the observed character)
    bpp::ColMatrix<double> fv;
    fv.resize(lenAlphabet, 1);
    bpp::MatrixTools::fill(fv, 0.0); // all zeros

    // get the gap position in the alphabet
    std::string ch(1, GAP_CHAR);
    int gapIndex = alphabet->charToInt(ch);

    fv(gapIndex, 0) = 1.0; // set gap position to 1

    // for all the gamma categories an array of fv values
    fv_empty_data_.resize(numCatg);

    // assign the indicator array to all gamma categories
    for (int catg = 0; catg < numCatg; catg++) {
        fv_empty_data_.at(catg) = fv;
    }

}

void PIPmsa::_setFVemptyNode(int numCatg,
                             PIPmsa *childL,
                             PIPmsa *childR,
                             std::vector<bpp::RowMatrix<double> > &PrL,
                             std::vector<bpp::RowMatrix<double> > &PrR) {

    fv_empty_data_.resize(numCatg);

    // array of lk (for each gamma rate) of a single column full of gaps
    for (int catg = 0; catg < numCatg; catg++) {

        // PrfvL = Pr_L * fv_L
        bpp::ColMatrix<double> PrfvL;
        bpp::MatrixTools::mult(PrL.at(catg), childL->fv_empty_data_.at(catg), PrfvL);

        // PrfvR = Pr_R * fv_R
        bpp::ColMatrix<double> PrfvR;
        bpp::MatrixTools::mult(PrR.at(catg), childR->fv_empty_data_.at(catg), PrfvR);

        // fv = PrfvL * PrfvR
        bpp::ColMatrix<double> fv;
        bpp::MatrixTools::hadamardMult(PrfvL, PrfvR, fv);

        fv_empty_data_.at(catg) = fv;
    }

}

void PIPmsa::_setFVsigmaEmptyLeaf(int numCatg) {

    // allocate memory ([numCatg] x 1)
    fv_empty_sigma_.resize(numCatg);

    for (int catg = 0; catg < numCatg; catg++) {
        // fv_empty_sigma = fv dot pi
        // fv_empty_sigma is always 0 at the leaves
        fv_empty_sigma_.at(catg) = 0.0;
    }

}

void PIPmsa::_setFVsigmaEmptyNode(int numCatg,
                                  PIPmsa *childL,
                                  PIPmsa *childR,
                                  double bL,
                                  double bR,
                                  const std::vector<double> &mu) {

    double zetaL;
    double zetaR;

    // resize to the number of categories
    fv_empty_sigma_.resize(numCatg);

    for (int catg = 0; catg < numCatg; catg++) {

        zetaL = exp(-mu.at(catg) * bL); // pure survival probability on the left child
        zetaR = exp(-mu.at(catg) * bR); // pure survival probability on the right child

        // fv_empty_sigma = dot(fv_empty,pi)
        // which corresponds to
        // fv_empty_sigma = not_survival_L * not_survival_R +
        //                  not_survival_L * survival_R * not_survival_below_R +
        //                  survival_L * not_survival_below_L * not_survival_R +
        //                  survival_L * not_survival_below_L * not_survival_R * not_survival_below_R
        fv_empty_sigma_.at(catg) = \
                            (1 - zetaL) * (1 - zetaR) + \
                            (1 - zetaL) * zetaR * childR->fv_empty_sigma_.at(catg) + \
                            zetaL * childL->fv_empty_sigma_.at(catg) * (1 - zetaR) + \
                            zetaL * childL->fv_empty_sigma_.at(catg) * zetaR * childR->fv_empty_sigma_.at(catg);

    }

}

void PIPmsa::_compressFv(std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed) {

    // compress an array of fv values

    int comprMSAlen = rev_map_compressed_seqs_.size();

    int id_map;

    fv_data_.resize(comprMSAlen);

    for (int i = 0; i < comprMSAlen; i++) {
        id_map = rev_map_compressed_seqs_.at(i);

        fv_data_.at(i) = fv_data_not_compressed.at(id_map);
    }

}

void PIPmsa::_compressFvSigma(std::vector<std::vector<double>> &fv_sigma_not_compressed) {

    // compress an array of fv_sigma values

    int comprMSAlen = rev_map_compressed_seqs_.size();

    int id_map;

    fv_sigma_.resize(comprMSAlen);

    for (int i = 0; i < comprMSAlen; i++) {
        id_map = rev_map_compressed_seqs_.at(i);

        fv_sigma_.at(i) = fv_sigma_not_compressed.at(id_map);
    }

}

void PIPmsa::_compressMSA(const bpp::Alphabet *alphabet) {

    auto sequences = new bpp::VectorSequenceContainer(alphabet);

    std::vector<std::string> seqs = PIPmsaUtils::siteContainer2sequenceVector(msa_);

    for (int i = 0; i < seqs.size(); i++) {
        sequences->addSequence(*(new bpp::BasicSequence(seqNames_.at(i),
                                                        seqs.at(i),
                                                        alphabet)), true);
    }

    auto siteContainer = new bpp::VectorSiteContainer(*sequences);
    auto siteContCompr = bpp::PatternTools::shrinkSiteSet(*siteContainer);
    auto map_seqs = bpp::PatternTools::getIndexes(*siteContainer, *siteContCompr);

    map_compressed_seqs_ = map_seqs;

    std::vector<int> rev_map_seqs = PIPmsaUtils::reverseMap(map_seqs);

    rev_map_compressed_seqs_ = rev_map_seqs;

}

void PIPmsa::_build_MSA(MSA_t &msaL, MSA_t &msaR) {

    // convert traceback path into an MSA

    // get dimension of the left/right MSA column
    int lenColL = msaL.at(0).size();
    int lenColR = msaR.at(0).size();

    int idx_i = 0;
    int idx_j = 0;
    for (int j = 0; j < traceback_path_.size(); j++) {

        if (traceback_path_.at(j) == (int) MATCH_STATE) {

            // in MATCH case concatenate left_column (from seq1) with right_column (from seq2)
            msa_.push_back(msaL.at(idx_i) + msaR.at(idx_j));
            idx_i++;
            idx_j++;

        } else if (traceback_path_.at(j) == (int) GAP_X_STATE) {

            // in GAPX case concatenate left_column (from seq1) with a column full of gaps (right)
            std::string gapCol(lenColR, GAP_CHAR);
            msa_.push_back(msaL.at(idx_i) + gapCol);
            idx_i++;

        } else if (traceback_path_.at(j) == GAP_Y_STATE) {

            // in GAPY case concatenate a column (left) full of gaps with right_column (from seq2)
            std::string gapCol(lenColL, GAP_CHAR);
            msa_.push_back(gapCol + msaR.at(idx_j));
            idx_j++;

        } else {
            LOG(FATAL) << "\nSomething went wrong during the traceback in function "
                          "pPIP::_build_MSA. Check call stack below.";
        }
    }

}

void PIPmsa::_computeLkEmptyLeaf(const bpp::progressivePIP *pPIP,
                                 std::vector<double> &iotasNode,
                                 std::vector<double> &betasNode){

    // compute the lk of an empty column at the leaf

    // allocate memory ([numCatg] x 1)
    lk_empty_.resize(pPIP->numCatg_);

    // only 1 column
    for (int catg=0; catg<pPIP->numCatg_; catg++) {
        // compute the lk of an empty column at the leaf
        lk_empty_.at(catg) = pPIP->rDist_->getProbability((size_t) catg) * \
                             iotasNode.at(catg) * (1 - betasNode.at(catg));
    }

}

void PIPmsa::_computeLkLeaf(const bpp::progressivePIP *pPIP,
                            std::vector<double> &iotasNode,
                            std::vector<double> &betasNode){

    // compute the lk at the leaf

    // get the number of gamma categories
    int numCatg = pPIP->numCatg_;

    // get the size of the compressed sequences
    int msaLen = _getCompressedMSAlength();

    // allocate memory ([site])
    log_lk_down_.resize(msaLen);

    // compute the marginal lk over all the gamma categories
    for(int site=0;site<msaLen;site++){
        // init to 0.0
        log_lk_down_.at(site) = 0.0;
        for (int catg=0; catg<numCatg; catg++) {
            log_lk_down_.at(site) += pPIP->rDist_->getProbability((size_t) catg) * \
                                     iotasNode.at(catg) * betasNode.at(catg) * \
                                     fv_sigma_.at(site).at(catg);
        }
        // compute the log lk
        log_lk_down_.at(site) = log(log_lk_down_.at(site));
    }

}

void PIPmsa::_compressLK(std::vector<double> &lk_down_not_compressed){

    // compress an array of lk values

    int comprMSAlen = _getCompressedMSAlength();

    int id_map;

    log_lk_down_.resize(comprMSAlen);

    for(int i=0;i<comprMSAlen;i++){
        id_map = rev_map_compressed_seqs_.at(i);
        log_lk_down_.at(i)=lk_down_not_compressed.at(id_map);
    }

}

void PIPmsa::_setTracebackPathLeaf() {

    // get the MSA size
    int MSAlen = msa_.size();

    // resize traceback path
    traceback_path_.resize(MSAlen);

    for (int i = 0; i < MSAlen; i++) {
        // assign MATCH STATE to all the sites
        traceback_path_.at(i) = (int) MATCH_STATE;
    }

}


//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
std::vector<std::string> PIPmsaUtils::siteContainer2sequenceVector(std::vector<bpp::MSAcolumn_t> &MSA) {

    std::vector<std::string> seqs;

    int len = MSA.size();
    int nseq = MSA.at(0).size();

    seqs.resize(nseq);
    for (int i = 0; i < nseq; i++) {
        std::string s;
        s.resize(len);
        for (int j = 0; j < len; j++) {
            s.at(j) = MSA.at(j).at(i);
        }
        seqs[i] = s;
    }

    return seqs;

};

std::vector<int> PIPmsaUtils::reverseMap(std::vector<int> &m) {

    std::vector<int> rev_m;

    for (int i = 0; i < m.size(); i++) {
        if ((m.at(i) + 1) > rev_m.size()) {
            if (m.at(i) - rev_m.size() > 0) {
                LOG(FATAL) << "\nERROR in reverse_map";
            }
            rev_m.push_back(i);
        }
    }

    return rev_m;
}

bpp::SiteContainer *PIPmsaUtils::PIPmsa2Sites(const bpp::Alphabet *alphabet,
                                              std::vector<std::string> &seqNames,
                                              std::vector<std::string> &MSA) {

    auto sequences = new bpp::VectorSequenceContainer(alphabet);

    int msaLen = MSA.size();

    int numLeaves = seqNames.size();
    for (int j = 0; j < numLeaves; j++) {
        std::string seqname = seqNames.at(j);
        std::string seqdata;
        seqdata.resize(msaLen);
        for (int i = 0; i < msaLen; i++) {
            seqdata.at(i) = MSA.at(i).at(j);
        }
        sequences->addSequence(*(new bpp::BasicSequence(seqname, seqdata, alphabet)), true);
    }

    return new bpp::VectorSiteContainer(*sequences);
}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
