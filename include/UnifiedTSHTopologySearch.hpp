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
 * @file TSHTopologySearch.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 05 02 2018
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
 * @see For more information visit: https://bitbucket.org/lorenzogatti89/castor/wiki/Home
 */
#ifndef CASTOR_TSHTOPOLOGYSEARCH_HPP
#define CASTOR_TSHTOPOLOGYSEARCH_HPP

#include <chrono>
#include <Bpp/Phyl/TopologySearch.h>
#include <Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h>
#include <Utree.hpp>
#include <TreeRearrangment.hpp>

#include "Utils.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_Generic.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_PIP.hpp"

namespace tshlib {

    class TreeSearch {

    private:
        bpp::AbstractHomogeneousTreeLikelihood *likelihoodFunc;
        double tshinitScore;
        double tshcycleScore;
        TreeSearchHeuristics tshSearchHeuristic;
        TreeRearrangmentOperations tshRearrangementCoverage;
        StartingNodeHeuristics tshStartingNodeMethod;
        double toleranceValue;
        int maxTSCycles;
        std::string scoringMethod;
        std::vector<int> performed_moves;
        int performed_cycles;
        int search_startingnodes;
        mutable tshlib::Utree *utree_;
        bool tshStatus;
        int threads_num;

    public:


        TreeSearch() {
            likelihoodFunc = nullptr;
            tshinitScore = -std::numeric_limits<double>::infinity();
            tshcycleScore = -std::numeric_limits<double>::infinity();
            tshSearchHeuristic = TreeSearchHeuristics::nosearch;
            tshRearrangementCoverage = TreeRearrangmentOperations::classic_Mixed;
            tshStartingNodeMethod = StartingNodeHeuristics::greedy;
            performed_cycles = 0;
            search_startingnodes = 0;
            scoringMethod = "";
            tshStatus = false;
            utree_ = nullptr;
            threads_num = 1;
        };

        ~TreeSearch() = default;

        void setUtree(tshlib::Utree *inUTree) {
            utree_ = inUTree;
        }

        Utree* getUtree() {
            return utree_ ?: nullptr;
        }

        void setMaxNumThreads(int numThreads){
            TreeSearch::threads_num = numThreads;
        }

        bool isTreeSearchSuccessful() const {
            return tshStatus;
        }

        void setTreeSearchStatus(bool tshStatus) {
            TreeSearch::tshStatus = tshStatus;
        }


        void setLikelihoodFunc(bpp::AbstractHomogeneousTreeLikelihood *in_likelihoodFunc) {
            likelihoodFunc = in_likelihoodFunc;
            Utree *lkUtree;
            if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
                lkUtree = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->getUtreeTopology();
            } else {
                lkUtree = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->getUtreeTopology();
            }

            setUtree(lkUtree);
        }

        void setInitialLikelihoodValue(double in_initialLikelihoodValue) {
            tshinitScore = in_initialLikelihoodValue;
        }

        void setTreeSearchStrategy(TreeSearchHeuristics in_tshStrategy, TreeRearrangmentOperations in_tshOperations) {
            tshSearchHeuristic = in_tshStrategy;
            tshRearrangementCoverage = in_tshOperations;
        }


        void setTolerance(double inTolerance) {
            toleranceValue = inTolerance;
        }

        void setMaxCycles(int inCycles){
            maxTSCycles = inCycles;
        }

        void setStartingNodeHeuristic(StartingNodeHeuristics in_tshStartingNodeMethod, int in_search_startingnodes){
            tshStartingNodeMethod = in_tshStartingNodeMethod;
            setStartingNodes(in_search_startingnodes);
        }

        StartingNodeHeuristics getStartingNodeHeuristic() const {
            return tshStartingNodeMethod;
        }


        bpp::AbstractHomogeneousTreeLikelihood *getLikelihoodFunc() const {
            return likelihoodFunc;
        }

        double getInitialLikelihoodValue() const {
            return tshinitScore;
        }

        TreeSearchHeuristics getTshStrategy() const {
            return tshSearchHeuristic;
        }

        TreeRearrangmentOperations getTshOperations() const {
            return tshRearrangementCoverage;
        }


        const std::string getScoringMethod() const {
            return scoringMethod;
        }

        void setScoringMethod(const std::string &inScoringMethod) {
            scoringMethod = inScoringMethod;
        }

        int getStartingNodes() const {
            return search_startingnodes;
        }

        void setStartingNodes(int in_search_startingnodes) {

            if (getUtree()) {
                if (utree_->listVNodes.size() < in_search_startingnodes) {

                    search_startingnodes = (int) utree_->listVNodes.size();

                   DLOG(WARNING) << "[TreeSearch::setStartingNodes] User requested too many initial seed nodes [" << in_search_startingnodes
                                 << "] to define candidate topology. Reset value to max number of nodes in the tree = "
                                 << search_startingnodes;

                } else {
                    search_startingnodes = in_search_startingnodes;

                }
            }else{
                LOG(ERROR) << "[TreeSearch::setStartingNodes] Utree has not been set for the current tree-search object. Call setUtree() first.";
            }
        }


        std::string getRearrangmentCoverageDescription() const {
            std::string rtToken;
            switch (tshRearrangementCoverage) {
                case tshlib::TreeRearrangmentOperations::classic_NNI:
                    rtToken = "NNI-like";
                    break;
                case tshlib::TreeRearrangmentOperations::classic_SPR:
                    rtToken = "SPR-like";
                    break;
                case tshlib::TreeRearrangmentOperations::classic_TBR:
                    rtToken = "TBR-like";
                    break;
                case tshlib::TreeRearrangmentOperations::classic_Mixed:
                    rtToken = "Complete";
                    break;
            }
            return rtToken;
        }

        std::string getTreeSearchStrategyDescription() const {
            std::string rtToken;
            switch (tshSearchHeuristic) {
                case tshlib::TreeSearchHeuristics::swap:
                    rtToken = "Swap";
                    break;
                case tshlib::TreeSearchHeuristics::phyml:
                    rtToken = "PhyML";
                    break;
                case tshlib::TreeSearchHeuristics::mixed:
                    rtToken = "Mixed(Swap+phyML)";
                    break;
                case tshlib::TreeSearchHeuristics::nosearch:
                    rtToken = "no-search";
                    break;
            }
            return rtToken;
        }

        std::string getStartingNodeHeuristicDescription() const {
            std::string rtToken;
            switch (tshStartingNodeMethod) {
                case tshlib::StartingNodeHeuristics::particle_swarm:
                    rtToken = "Particle Swarm";
                    break;
                case tshlib::StartingNodeHeuristics::hillclimbing:
                    rtToken = "Hillclimbing";
                    break;
                case tshlib::StartingNodeHeuristics::greedy:
                    rtToken = "Greedy";
                    break;
                case tshlib::StartingNodeHeuristics::undef:
                    rtToken = "undef";
                    break;
            }
            return rtToken;
        }


        double executeTreeSearch();

    protected:

        /*!
         * @brief This method execute the cycle of tree search according to the settings selected by the user
         * @return the score of the tree topology found during the tree search
         */
        double iterate();

        /*!
         * @brief This method defines the candidate rearrangements to be performed on the topology
         * @return A list of candidate topologies (in the form of rearrangement operations)
         */
        tshlib::TreeRearrangment *defineMoves();

        /*!
         * @brief testCandidateMoves method tests all the candidate tree topologies in the rearrangement list (or set), and it saves the score for each of them
         * @param candidateMoves pointer to the list of candidate rearrangements on a fixed topology
         */
        void testMoves(tshlib::TreeRearrangment *candidateMoves);

        /*!
         * @brief This method preallocates the required amount of memory to test each candidate tree-rearrangement.
         *        Memory can have concurrent threads in reading and writing mode.
         * @param candidateMoves stack of the moves to test
         */
        void allocateTemporaryLikelihoodData(int numThreads = 1);

        void deallocateTemporaryLikelihoodData(int numThreads = 1);


        std::string debugStackTraceMove(Move *move, Utree *_utree,
                                               std::vector < int > listNodesInvolved,
                                               std::vector < int > updatedList,
                                               double initLK = 0,
                                               double moveLK = 0);
    };
}


#endif //CASTOR_TSHTOPOLOGYSEARCH_HPP
