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
 * @file pPIP.hpp
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
 * @see For more information visit: https://bitbucket.org/acg-team/minijati/wiki/Home
 */
#ifndef MINIJATI_PPIP_HPP
#define MINIJATI_PPIP_HPP

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Node.h>
#include <random>
#include <Utree.hpp>
#include <Bpp/Phyl/Model/SubstitutionModel.h>

#include "Utils.hpp"
#include "Utilities.hpp"

#define SMALL_DOUBLE 1e-8
#define LEFT 0
#define RIGHT 1

namespace bpp {

    class pPIP {

    public:

        typedef std::string MSAcolumn_t;           // MSA column type
        typedef std::vector<MSAcolumn_t> MSA_t;    // MSA as vector of columns
        typedef std::vector<MSA_t> MSAensemble_t;  // for SB: set (vector) of MSAs
        //typedef std::string TracebackPath_t;       // traceback path type
        //typedef std::vector<TracebackPath_t> TracebackEnsemble_t; //for SB: set (vector) of traceback paths

        pPIP(tshlib::Utree *utree,              // tshlib:: tree
             bpp::Tree *tree,                   // bpp::tree
             bpp::SubstitutionModel *smodel,    // extended substitution model
             UtreeBppUtils::treemap &inTreeMap, // bpp::Node * <-> tshlib::VirtualNode *
             bpp::SequenceContainer *sequences, // un-aligned input sequences
             bpp::DiscreteDistribution *rDist,  // distribution for rate variation among sites
             long seed,                         // seed for the random numbers generation
             int num_sb);

        ~pPIP(){};

        void PIPAligner(std::vector<tshlib::VirtualNode *> &list_vnode_to_root,
                        bool local,
                        bool flag_map,
                        bool flag_fv);

        std::vector< std::string > getMSA(bpp::Node *node,int idx);

        std::vector<double> getScore(bpp::Node *node);

        std::vector< std::string > getSeqnames(bpp::Node *node);

        bpp::Node *getRootNode();

        const Alphabet *getAlphabet() const;

        void setSubstModel(bpp::SubstitutionModel *smodel);

        void setTree(const Tree *tree);

    protected:

    private:

        tshlib::Utree *utree_;                   // tshlib:: tree
        bpp::TreeTemplate<bpp::Node> *tree_;     // bpp::tree
        bpp::SubstitutionModel *substModel_;
        long seed_;                              //jatiapp seed for the random numbers generation

    private:

        //***************************************************************************************
        // extended substitution model
        mutable UtreeBppUtils::treemap treemap_; // bpp::Node * <-> tshlib::VirtualNode *
        bpp::SequenceContainer *sequences_;      // un-aligned input sequences
        bpp::DiscreteDistribution *rDist_;       // distribution for rate variation among sites

        std::vector<std::vector<double>> iotasNode_; //map of nodeIDs and vector of iotas (1 for each rate (Gamma,...) category
        std::vector<std::vector<double>> betasNode_; //map of nodeIDs and vector of betas (1 for each rate (Gamma,...) category
        std::vector<std::vector<bpp::RowMatrix<double> > >prNode_; // map of NodeIDs of Pr = exp(branchLength * rate * Q), rate under Gamma distribution

        std::vector< std::vector<std::string> > seqNames_;          // vector[nodeId] of sequence names (MSAs seq. names at each internal node) node
        std::vector< std::vector<MSA_t> > MSA_;                                   // vector[nodeId] MSA at each node
        std::vector< std::vector<vector<int> > > traceback_path_;              // vector[nodeId] of traceback paths (1 at each internal node)
        std::vector< std::vector< vector< vector<int> > > > traceback_map_;

        std::vector< std::vector< std::vector<int> > > subMSAidx_;

        std::vector< std::vector<double> > score_;                                // vector[nodeId] of likelihood score
        double lambda0_;                                            // original lambda (no Gamma distribution)
        double mu0_;                                                // original mu (no Gamma distribution)
        std::vector<double> lambda_;                               // vector[rate] of lambda rate with Gamma distribution
        std::vector<double> mu_;                                   // vector[rate] of mu rate with Gamma distribution
        std::vector<double> nu_;                                   // vector[rate] of nu (normalizing constant) with Gamma distribution
        double tau_;                                               // total tree length

        int num_sb_; // number of stochastic backtracking solutions

        std::vector<double> taus_;

        std::vector< std::vector<double> > nus_;

        std::vector< vector<double> > log_lk_down_;                      //each node a vector of lk
        std::vector< vector<double> > log_lk_empty_down_;                //each node a vector of lk_empty (for each gamma category)

        std::vector< std::vector< vector< vector< bpp::ColMatrix<double> > > > > fv_data_; // [node][site][catg][fv]
        std::vector< vector< bpp::ColMatrix<double> > > fv_empty_data_; // [node][catg][fv]
        std::vector< std::vector< vector<int> > > map_compressed_seqs_; // [node][idx]
        std::vector< std::vector< vector<int> > > rev_map_compressed_seqs_; // [node][idx]

        std::vector< std::vector< std::vector< std::vector<double> > > > fv_sigma_; // [node][site][catg]
        std::vector< std::vector<double> > fv_empty_sigma_; // [node][catg]

        bpp::ColMatrix<double> pi_;                                // steady state base frequencies

        const bpp::Alphabet *alphabet_;                            // extended alphabet (alphabet U {'-'})

        long alphabetSize_;                                        // original alphabet size

        long extendedAlphabetSize_;                                // extended alphabet size
        //***************************************************************************************
        void _reserve(std::vector<tshlib::VirtualNode *> &nodeList);

        std::vector<double> _computeNu(int nodeID);

        void _computeNus(std::vector<tshlib::VirtualNode *> &nodeList);

        void _setTree(const Tree *tree);

        void _setLambda(double lambda);

        void _setMu(double mu);

        void _setPi(const Vdouble &pi);

        double _computeTauRecursive(tshlib::VirtualNode *vnode);

        double _computeTau(tshlib::VirtualNode *vnode);

        void _computeTaus(std::vector<tshlib::VirtualNode *> &nodeList);

        void _setAllIotas(bpp::Node *node,bool local_root);

        void _setAllBetas(bpp::Node *node,bool local_root);

        void _getPrFromSubstutionModel(std::vector<tshlib::VirtualNode *> &listNodes);

        void _setFVleaf(bpp::Node *node);

        void _setFVemptyLeaf(bpp::Node *node);

        void _setFVsigmaEmptyLeaf(bpp::Node *node);

        void _setFVsigmaLeaf(bpp::Node *node);

        void _set_lk_leaf(bpp::Node *node);

        void _set_lk_empty_leaf(bpp::Node *node);

        void _compressMSA(bpp::Node *node,int idx_sb);

        void _compress_lk_components(bpp::Node *node,
                                     std::vector<double> &lk_down_not_compressed,
                                     std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed,
                                     int idx_sb);

        void _compress_Fv(bpp::Node *node,
                          std::vector<std::vector<double>> &fv_sigma_not_compressed,
                          std::vector<vector<bpp::ColMatrix<double> > > &fv_data_not_compressed,
                          int idx_sb);

        bool is_inside(int x0,
                       int y0,
                       int xf,
                       int yf,
                       int xt,
                       int yt);

        void set_indices_M(int &up_corner_i,
                           int &up_corner_j,
                           int &bot_corner_i,
                           int &bot_corner_j,
                           int level,
                           int h,
                           int w);

        void set_indices_X(int &up_corner_i,
                           int &up_corner_j,
                           int &bot_corner_i,
                           int &bot_corner_j,
                           int level,
                           int h,
                           int w);

        void set_indices_Y(int &up_corner_i,
                           int &up_corner_j,
                           int &bot_corner_i,
                           int &bot_corner_j,
                           int level,
                           int h,
                           int w);

        signed long get_indices_M(int nx,
                                  int ny,
                                  int up_corner_i,
                                  int up_corner_j,
                                  int bot_corner_i,
                                  int bot_corner_j,
                                  int m,
                                  int h,
                                  int w);

        signed long get_indices_X(int nx,
                                  int ny,
                                  int up_corner_i,
                                  int up_corner_j,
                                  int bot_corner_i,
                                  int bot_corner_j,
                                  int m,
                                  int h,
                                  int w);

        signed long get_indices_Y(int nx,
                                  int ny,
                                  int up_corner_i,
                                  int up_corner_j,
                                  int bot_corner_i,
                                  int bot_corner_j,
                                  int m,
                                  int h,
                                  int w);

        void set_indices_T(int &up_corner_i,
                           int &up_corner_j,
                           int &bot_corner_i,
                           int &bot_corner_j,
                           int level,
                           int h,
                           int w);

        void reset_corner(int &up_corner_i,
                          int &up_corner_j,
                          int &bot_corner_i,
                          int &bot_corner_j,
                          int h,
                          int w);

        int get_indices_T(int nx,
                          int ny,
                          int up_corner_i,
                          int up_corner_j,
                          int bot_corner_i,
                          int bot_corner_j,
                          int m,
                          int h,
                          int w);

        bool _index_of_max(double m,
                           double x,
                           double y,
                           double epsilon,
                           std::default_random_engine &generator,
                           std::uniform_real_distribution<double> &distribution,
                           bool flag_RAM,
                           int &index,
                           double &val);

        double max_of_three(double a,
                            double b,
                            double c,
                            double epsilon,
                            bool flag_RAM);

        bool checkboundary(int up_corner_i,
                           int up_corner_j,
                           int bot_corner_i,
                           int bot_corner_j,
                           int h,
                           int w);

        std::string createGapCol(int len);

        void _build_MSA(bpp::Node *node, int idx_sb);

        void _setMSAsequenceNames(bpp::Node *node);

        void _setMSAsequenceNames(bpp::Node *node,
                                  std::string seqname);

        void _setMSAleaves(bpp::Node *node,
                           const std::string &sequence);

        void _setTracebackPathleaves(bpp::Node *node);

        void _setSubMSAindexLeaves(bpp::Node *node);

        bpp::ColMatrix<double> fv_observed(MSAcolumn_t &s, int &idx);

        bpp::ColMatrix<double> computeFVrec(bpp::Node *node,
                                            MSAcolumn_t &s,
                                            int &idx,
                                            int catg);

        void allgaps(bpp::Node *node,MSAcolumn_t &s, int &idx,bool &flag);

        double compute_lk_gap_down(bpp::Node *node,MSAcolumn_t &s,int catg);

        std::vector<double> computeLK_GapColumn_local(bpp::Node *node,
                                                      MSAcolumn_t &sL,
                                                      MSAcolumn_t &sR,
                                                      bool flag_RAM);

        std::vector<double> computeLK_GapColumn_local(int nodeID,
                                                      int sonLeftID,
                                                      int sonRightID,
                                                      std::vector< bpp::ColMatrix<double> > &fvL,
                                                      std::vector< bpp::ColMatrix<double> > &fvR,
                                                      std::vector< bpp::ColMatrix<double> > &Fv_gap);

        std::vector<double> computeLK_GapColumn_local(int nodeID,
                                                            int sonLeftID,
                                                            int sonRightID,
                                                            std::vector< bpp::ColMatrix<double> > &fvL,
                                                            std::vector< bpp::ColMatrix<double> > &fvR,
                                                            std::vector< bpp::ColMatrix<double> > &Fv_gap,
                                                            std::vector<double> &fv_empty_sigma_,
                                                            std::vector<double> &lk_empty_down_L,
                                                            std::vector<double> &lk_empty_down_R);

        void _compute_lk_empty_down_rec(bpp::Node *node,std::vector<double> &lk);

        std::vector<double> _compute_lk_empty_down(bpp::Node *node);

        double _compute_lk_down_rec(bpp::Node *node,int idx,double lk,int position);

        std::vector<double> _compute_lk_down(bpp::Node *node,int idx_sb);

        double _compute_lk_down(bpp::Node *node, MSAcolumn_t &s, int catg);

        double computeLK_MXY_local(double log_phi_gamma,
                                   double valM,
                                   double valX,
                                   double valY,
                                   double log_pr);

        double computeLK_M_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 std::string &sL,
                                 std::string &sR,
                                 int m,
                                 std::map<MSAcolumn_t, double> &lkM,
                                 std::vector< std::vector<double> > &lkM_pattern,
                                 bool flag_map,
                                 bool flag_RAM,
                                 bool flag_pattern);

        double computeLK_M_local(int nodeID,
                                 int sonLeftID,
                                 int sonRightID,
                                 std::vector< bpp::ColMatrix<double> > &fvL,
                                 std::vector< bpp::ColMatrix<double> > &fvR,
                                 std::vector< bpp::ColMatrix<double> > &Fv_M_ij,
                                 std::vector<double> &Fv_sigma_M_ij);

        double computeLK_X_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn_t &sL,
                                 MSAcolumn_t &col_gap_R,
                                 int m,
                                 std::map<MSAcolumn_t, double> &lkX,
                                 std::vector< std::vector<double> > &lkX_pattern,
                                 bool flag_map,
                                 bool flag_RAM,
                                 bool flag_pattern);

        double computeLK_X_local(int nodeID,
                                 int sonLeftID,
                                 int sonRightID,
                                 std::vector< bpp::ColMatrix<double> > &fvL,
                                 std::vector< bpp::ColMatrix<double> > &fvR,
                                 std::vector< bpp::ColMatrix<double> > &Fv_X_ij,
                                 std::vector<double> &Fv_sigma_X_ij);

        double computeLK_Y_local(double NU,
                                 double valM,
                                 double valX,
                                 double valY,
                                 bpp::Node *node,
                                 MSAcolumn_t &col_gap_L,
                                 MSAcolumn_t &sR,
                                 int m,
                                 std::map<MSAcolumn_t, double> &lkY,
                                 std::vector< std::vector<double> > &lkY_pattern,
                                 bool flag_map,
                                 bool flag_RAM,
                                 bool flag_pattern);


        double computeLK_Y_local(int nodeID,
                                 int sonLeftID,
                                 int sonRightID,
                                 std::vector< bpp::ColMatrix<double> > &fvL,
                                 std::vector< bpp::ColMatrix<double> > &fvR,
                                 std::vector< bpp::ColMatrix<double> > &Fv_Y_ij,
                                 std::vector<double> &Fv_sigma_Y_ij);

        void DP3D_PIP(bpp::Node *node, bool local,bool flag_map);

        void extract_N_best(bpp::Node *node,int totalSubMSAs);

        void startingLevelSB(std::vector< vector< vector<double> > > &Log3DM,
                             std::vector< vector< vector<double> > > &Log3DX,
                             std::vector< vector< vector<double> > > &Log3DY,
                             double epsilon,
                             std::default_random_engine &generator,
                             std::uniform_real_distribution<double> &distribution,
                             int d,
                             int h,
                             int w,
                             int &lev,
                             double &val,
                             int &state);

        //void DP3D_PIP_no_gamma(bpp::Node *node, bool local,bool flag_map);

//        void DP3D_PIP_RAM(bpp::Node *node,
//                          bool local,
//                          bool flag_map,
//                          bool flag_pattern);

        void DP3D_PIP_RAM_FAST(bpp::Node *node);

        void DP3D_PIP_RAM_FAST_SB_cross(bpp::Node *node);

        void DP3D_PIP_RAM_FAST_SB(bpp::Node *node,int msa_idx_L, int msa_idx_R, int &position);

//        void DP3D_PIP_SB(bpp::Node *node,
//                         UtreeBppUtils::treemap *tm,
//                         double gamma_rate,
//                         bool local,
//                         double temperature,
//                         int num_SB);

    };

}

namespace pPIPUtils {

    // convert MSA PIP into sites container
    bpp::SiteContainer *pPIPmsa2Sites(bpp::pPIP *progressivePIP,int idx_sb);

    // sum of logs
    double add_lns(double a_ln,double b_ln);

    void max_val_in_column(std::vector<std::vector<std::vector<double>>> &M,int depth, int height, int width, double &val, int &level);

    std::vector<std::string> siteContainer_2_sequence_vector(std::vector<bpp::pPIP::MSAcolumn_t> &MSA);

    std::vector<int> reverse_map(std::vector<int> &m);

}

#endif //MINIJATI_PPIP_HPP
