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
 * @file progressivePIP.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 19 02 2018
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

#ifndef MINIJATI_PROGRESSIVEPIP_HPP
#define MINIJATI_PROGRESSIVEPIP_HPP

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Node.h>
#include <random>
#include <Utree.hpp>
#include <glog/logging.h>

#include <Bpp/Phyl/BipartitionList.h>

#include "Utilities.hpp"
#include "Utils.hpp"

#define SMALL_DOUBLE 1e-8
#define GAP_CHAR '-'
#define ERR_STATE (-999)
#define DBL_EPSILON std::numeric_limits<double>::min()




#define STFT 1
//#define COMPRESS_TR 1
//#define EXTRA_TR 1
//#define TR_GAP 1

#ifdef COMPRESS_TR
    #define STOP_STATE  0b00000000
    #define MATCH_STATE 0b00000001
    #define GAP_X_STATE 0b00000010
    #define GAP_Y_STATE 0b00000011
#else

    #ifdef EXTRA_TR

        #define MATCH_STATE 1
        #define GAP_X_STATE 2
        #define GAP_Y_STATE 3
        #define STOP_STATE 4

        #define MX_STATE 5
        #define MY_STATE 6
        #define MXY_STATE 7
        #define XM_STATE 8
        #define XY_STATE 9
        #define XMY_STATE 10
        #define YM_STATE 11
        #define YX_STATE 12
        #define YMX_STATE 13

    #else
        #define MATCH_STATE 1
        #define GAP_X_STATE 2
        #define GAP_Y_STATE 3
        #define STOP_STATE 4
    #endif

#endif








#define LEFT 0
#define RIGHT 1

namespace bpp {

    //*******************************
    // DP3D versions
    // CPU: more CPU less memory
    // RAM: less CPU more memory
    // SB: stochastic backtracing version
    enum enumDP3Dversion{CPU,RAM,SB};
    //*******************************

    class CompositePIPnode; // forward declaration
    class PIPnode; // forward declaration

    //*******************************************************************************************
    class progressivePIP {

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

        tshlib::Utree *utree_; // tshlib:: tree

        bpp::TreeTemplate<bpp::Node> *tree_; // bpp::tree

        long seed_; //jatiapp seed for the random numbers generation

        mutable UtreeBppUtils::treemap treemap_; // bpp::Node * <-> tshlib::VirtualNode *

        double lambda0_; // original lambda (no Gamma distribution)

        double mu0_; // original mu (no Gamma distribution)

        PIPnode *PIPnodeRoot_; // root PIPnode (root of the tree of PIPnodes)

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

        void _setLambda(double lambda); // set lambda for different rates

        void _setMu(double mu); // set mu for different rates

        void _setPi(const Vdouble &pi); // set Pi

        void _computeTauRec_(PIPnode *pipnode); // recursive computation of the total tree length and
                                                // of the total left/right subtree length

        //***************************************************************************************

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        CompositePIPnode *compositePIPaligner_; // array of PIPnodes

        unsigned long numNodes_; // total number of nodes

        unsigned long numCatg_; // number of discrete categories (gamma distribution)

        int num_sb_; // number of sub-optimal solutions saved at each DP instance

        double temperature_; // to tune the greedyness of the sub-optimal solution

        std::vector<double> lambda_; // vector[rate] of lambda rate with Gamma distribution

        std::vector<double> mu_; // vector[rate] of mu rate with Gamma distribution

        bpp::DiscreteDistribution *rDist_; // distribution for rate variation among sites

        bpp::SubstitutionModel *substModel_; // substitution model

        bpp::SequenceContainer *sequences_; // un-aligned input sequences

        const bpp::Alphabet *alphabet_; // extended alphabet (alphabet U {'-'})

        long alphabetSize_; // original alphabet size

        long extendedAlphabetSize_; // extended alphabet size

        bpp::ColMatrix<double> pi_; // steady state base frequencies

        double tau_; // total tree length

        std::vector<double> nu_; // normalizing Poisson intensity for each gamma category

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        // constructor
        progressivePIP(tshlib::Utree *utree,              // tshlib:: tree
                       bpp::Tree *tree,                   // bpp::tree
                       bpp::SubstitutionModel *smodel,    // extended substitution model
                       UtreeBppUtils::treemap &inTreeMap, // bpp::Node * <-> tshlib::VirtualNode *
                       bpp::SequenceContainer *sequences, // un-aligned input sequences
                       bpp::DiscreteDistribution *rDist,  // distribution for rate variation among sites
                       long seed);                         // seed for the random numbers generation


        ~progressivePIP(){}; // destructor

        void _initializePIP(std::vector<tshlib::VirtualNode *> &list_vnode_to_root, // list of nodes
                            enumDP3Dversion DPversion, // DP3D version
                            int num_sb, // number of sub. optimal solutions (MSAs)
                            double temperature);  // to tune the greedyness of the sub-optimal solution

        void PIPnodeAlign();

        long getSeed() const { return seed_; }; // return seed for the random number generation

        void _setAllIotas(); // compute all the insertion probabilities (iota function)

        void _setAllBetas(); // compute all the survival probabilities (beta function)

        void _setAllAlphas(); // alpha(v) = sum_from_v_to_root ( iota * beta * zeta )
                              // zeta = exp(- mu *b ) is the "pure" survival probability

        void _setAllEtas();

        bpp::Node *getBPProotNode(){ return tree_->getRootNode(); }; // get the root of the tree

        PIPnode *getPIPnodeRootNode(){ return PIPnodeRoot_; }; // get the PIPnodeRoot of the PIPnode tree

        const Alphabet *getAlphabet() const { return alphabet_; }; // get the alphabet

        void _buildPIPnodeTree(); // build a binary tree of PIPnode that dictates the alignmanet order

        void _computeTau_(); // compute the total tree length and the length of the left/right subtree
                             // of the tree rooted at a given PIPnode

        void _computeNu(); // compute the normalizing Poisson intensity (expected MSA length)

        void _computeLengthPathToRoot(); // compute the distance from any node to the root

    };

}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************
namespace progressivePIPutils {

    // sum of logs
    double add_lns(double a_ln,double b_ln);

}
//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************

#endif //MINIJATI_PROGRESSIVEPIP_HPP
