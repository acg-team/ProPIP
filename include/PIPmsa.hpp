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
 * @file PIPmsa.hpp
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

#ifndef MINIJATI_PIPMSA_HPP
#define MINIJATI_PIPMSA_HPP

#include <string>
#include <vector>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Node.h>
#include "progressivePIP.hpp"

namespace bpp {

    //*******************************************************************************************

    typedef std::string MSAcolumn_t; // MSA column type
    typedef std::vector<MSAcolumn_t> MSA_t; // MSA as vector of columns

    //*******************************************************************************************
    class PIPmsa {

    private:

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        double score_; // lk of the MSA

        std::vector<std::string> seqNames_; // vector of strings (sequence names)

        MSA_t msa_; // vector of columns (strings). Each column is a site

        std::vector<int> traceback_path_;   // MSA traceback_path as vector of int where
        // 1: match state
        // 2: gapx state
        // 3: gapy state
        // 4: stop state (error or un-reachable entry)

        std::vector<int> map_compressed_seqs_; // map (columns position) between the compressed sequence and
        // the original sequence

        std::vector<int> rev_map_compressed_seqs_; // reverse map (columns position) between the compressed
        // sequence and the original sequence

        std::vector<std::vector<bpp::ColMatrix<double>>> fv_data_; // fv values (columns): [site][gamma cat.][state]

        std::vector<std::vector<double> > fv_sigma_; // fv_sigma = fv dot pi : [site][gamma cat.]

        std::vector<bpp::ColMatrix<double> > fv_empty_data_; // fv column of an empty column [gamma cat.][state]

        std::vector<double> fv_empty_sigma_; // fv_empty_sigma = fv_empty dot pi [gamma cat.]

        std::vector<double> log_lk_down_; // array of lk [site] (marginalized over gamma cat.)

        std::vector<double> lk_empty_; // array of lk of an empty column (for each gamma category)

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPmsa() { score_ = -std::numeric_limits<double>::infinity(); }; // constructor

        ~PIPmsa(){};

        void _setSeqNameLeaf(std::string &seqName); // assign the sequence name at a leaf

        void _setMSAleaf(const bpp::Sequence *sequence); // assign the MSA at a leaf (which corresponds
        // to the input sequence)

        // set sequences name at an internal node
        void _setSeqNameNode(
                std::vector<std::string> &seqNamesL,  // array of sequence names (strings) coming from the left subtree
                std::vector<std::string> &seqNamesR); // array of sequence names (strings) coming from the right subtree

        int _getMSAlength() { return msa_.size(); } // return the MSA length

        int _getCompressedMSAlength() { return rev_map_compressed_seqs_.size(); }; // return the compressed MSA length

        double _getScore() { return score_; } // get the score (lk)

        MSA_t *_getMSA() { return &msa_; } // get the MSA

        std::vector<std::string> *_getseqNames() { return &seqNames_; }; // get the sequences name

        void _setTracebackPathLeaf(); // initialized the traceback path (match state for all the sites)

        std::vector<std::string> getSequenceNames() { return seqNames_; }; // return the sequences name

        void _setFVleaf(int numCatg,
                        const bpp::Alphabet *alphabet);

        void _setFVsigmaLeaf(int numCatg,
                             const bpp::ColMatrix<double> &pi);

        void _setFVemptyNode(int numCatg,
                             PIPmsa *childL,
                             PIPmsa *childR,
                             std::vector<bpp::RowMatrix<double> > &PrL,
                             std::vector<bpp::RowMatrix<double> > &PrR);

        void _setFVsigmaEmptyLeaf(int numCatg);

        void _setFVsigmaEmptyNode(int numCatg,
                                  PIPmsa *childL,
                                  PIPmsa *childR,
                                  double bL,
                                  double bR,
                                  const std::vector<double> &mu);

        void _setFVemptyLeaf(int numCatg,
                             const bpp::Alphabet *alphabet);

        void _computeLkEmptyLeaf(const bpp::progressivePIP *pPIP,
                                 std::vector<double> &iotasNode,
                                 std::vector<double> &betasNode);

        void _computeLkLeaf(const bpp::progressivePIP *pPIP,
                                    std::vector<double> &iotasNode,
                                    std::vector<double> &betasNode);

        void _compressFv(std::vector<std::vector<bpp::ColMatrix<double> > > &fv_data_not_compressed); // compress an array of fv values

        void _compressFvSigma(std::vector<std::vector<double>> &fv_sigma_not_compressed); // compress an array of fv_sigma values

        void _compressMSA(const bpp::Alphabet *alphabet); // compress its MSA

        void _compressLK(std::vector<double> &lk_down_not_compressed); // store the compressed arry of lk

        void _build_MSA(MSA_t &msaL, MSA_t &msaR); // build an new MSA using the left/right children MSA and
                                                   // the traceback path

    };

}

//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
namespace PIPmsaUtils {

    std::vector<std::string> siteContainer2sequenceVector(std::vector<bpp::MSAcolumn_t> &MSA);

    std::vector<int> reverseMap(std::vector<int> &m); // flip a vector (left to right)

    bpp::SiteContainer *PIPmsa2Sites(const bpp::Alphabet *alphabet,
                                     std::vector<std::string> &seqNames,
                                     std::vector<std::string> &MSA);

}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************

#endif //MINIJATI_PIPMSA_HPP
