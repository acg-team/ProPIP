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
 * @file TSHTopologySearch.cpp
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
#include <random>
#include <glog/logging.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_INTELTBB
#include <tbb/parallel_for.h>
#endif

#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Io/Newick.h>

#include "UnifiedTSHTopologySearch.hpp"


using namespace bpp;

tshlib::TreeRearrangment *tshlib::TreeSearch::defineMoves() {

    // Initialise a new rearrangement list
    auto _move__set = new tshlib::TreeRearrangment;
    _move__set->setTree(utree_);
    _move__set->setTreeCoverage(tshRearrangementCoverage);
    _move__set->setStrategy(tshSearchHeuristic);
    _move__set->initialize();

    std::vector<tshlib::VirtualNode *> _seed__nodes((int) search_startingnodes);

    // instantiating the method to define the number of initial nodes on which computing the candidate moves
    switch (tshStartingNodeMethod) {

        case tshlib::StartingNodeHeuristics::greedy:

            // Generate candidate list of possible moves given the tree topology and the rearrangement operation type
            for (auto &node:utree_->listVNodes) {

                // Print node description with neighbors
                //DVLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

                _move__set->setSourceNode(node->getVnode_id());

                // Get all the target nodes with distance == radius from the source node
                // excluding the starting node (false)
                // do not duplicate moves in the list
                _move__set->defineMoves(false, false);

            }
            break;

        case tshlib::StartingNodeHeuristics::hillclimbing:

            RandomTools::getSample(utree_->listVNodes, _seed__nodes);

            // Generate candidate list of possible moves given the node topology and the rearrangement operation type
            for (auto &node:_seed__nodes) {

                // Print node description with neighbors
                //DVLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

                _move__set->setSourceNode(node->getVnode_id());

                // Get all the target nodes with distance == radius from the source node
                // excluding the starting node (false)
                // do not duplicate moves in the list
                _move__set->defineMoves(false, false);

            }
            break;

        case tshlib::StartingNodeHeuristics::particle_swarm:
        case tshlib::StartingNodeHeuristics::undef:
            LOG(FATAL) << "[TSH Optimisation] Strarting node heuristic not defined. Execution aborted!";
            break;
    }

    // Print the list of moves for the current P node (source node)
    //rearrangmentList.printMoves();

    return _move__set;
}


double tshlib::TreeSearch::executeTreeSearch() {
    double newScore = 0;
    setTreeSearchStatus(false);
    // Store initial likelihood value
    setInitialLikelihoodValue(likelihoodFunc->getLogLikelihood());
    //tshinitScore =  likelihoodFunc->getLogLikelihood();
    tshcycleScore = tshinitScore;

    // Remove root on the tree topology
    utree_->removeVirtualRootNode();

    DLOG(INFO) << "[TSH Optimisation] Initial topology: " << utree_->printTreeNewick(true);
    DLOG(INFO) << "[TSH Optimisation] Initial likelihood: " << TextTools::toString(tshinitScore, 15);
    DLOG(INFO) << "[TSH Optimisation] algorithm: " << getStartingNodeHeuristicDescription()
               << ", moves: " << getRearrangmentCoverageDescription() << ", "
               << "starting nodes: " + std::to_string(search_startingnodes);

    // execute tree search according to settings
    newScore = iterate();

    // Add root on the tree topology
    utree_->addVirtualRootNode();

    return newScore;
}


void tshlib::TreeSearch::testMoves(tshlib::TreeRearrangment *candidateMoves) {

    bool status = false;

    // Allocate memory
    allocateTemporaryLikelihoodData(threads_num);

    // Count moves performed per each cycle
    std::vector<int> _count__moves_cycle(candidateMoves->getNumberOfMoves(), 0);

    // Threads should run the following tasks in parallel
    int thread_id, i;
   
#ifdef USE_OPENMP
    // Set the number of threads to use
    omp_set_num_threads(threads_num);
    #pragma omp parallel for private(thread_id, i)
    for (i = 0; i < candidateMoves->getNumberOfMoves(); i++) {
#elif USE_INTELTBB
//  tbb::task_scheduler_init init(threads_num);
//  tbb::parallel_for( tbb::blocked_range<int>(0, candidateMoves->getNumberOfMoves()), [&](tbb::blocked_range<int> r)
//  {
//        for (int i=r.begin(); i<r.end(); ++i){
#else
	for (i = 0; i < candidateMoves->getNumberOfMoves(); i++) {
#endif   
        // ------------------------------------
        // Get the number of the thread
#ifdef USE_OPENMP
        thread_id = omp_get_thread_num();
#elif USE_INTELTBB
        //thread_id = tbb::this_tbb_thread::get_id();
#else 
        thread_id = 0;
#endif
        // ------------------------------------
        // Get move reference
        tshlib::Move *currentMove = candidateMoves->getMove(i);

        // ------------------------------------
        // Copy tshlib::utree structure locally (each thread should have the same reference copy)
        auto _thread__topology = new tshlib::Utree((*utree_));

        // ------------------------------------
        // Initialise local variables
        double moveLogLK = 0;

        // ------------------------------------
        // Prepare the list of nodes involved in the move (Required here!)
        VirtualNode *_node__source_id = _thread__topology->getNode(currentMove->getSourceNode());
        VirtualNode *_node__target_id = _thread__topology->getNode(currentMove->getTargetNode());

        std::vector<int> listNodesWithinPath = _thread__topology->computePathBetweenNodes(_node__source_id, _node__target_id);
        std::vector<int> updatedNodesWithinPath = candidateMoves->updatePathBetweenNodes(i, listNodesWithinPath);

        // ------------------------------------
        // Log move details
        DVLOG(1) << "[TSH Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1 << " ] "
                 << "Move [" << _thread__topology->getNode(currentMove->getSourceNode())->getNodeName() << "->"
                 << _thread__topology->getNode(currentMove->getTargetNode())->getNodeName() << "]\t"
                 << currentMove->getClass() << "\t(" << currentMove->getRadius() << ")" << "\tdirection: "
                 << currentMove->getDirection();

        // ------------------------------------
        // Apply the move
        status = candidateMoves->applyMove(i, (*_thread__topology));

        if (status) {

            DVLOG(1) << "[TSH Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1
                     << " ] [A] MOVE#" << currentMove->getUID() << " | " << _thread__topology->printTreeNewick(true);
            DVLOG(1) << "[TSH Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1
                     << " ] [R] MOVE#" << currentMove->getUID() << " | " << utree_->printTreeNewick(true, false);

            // ------------------------------------
            // Print root reachability from every node (includes rotations)
            // inputTree->_testReachingPseudoRoot();

            // ------------------------------------
            // Print tree on file
            //inputTree->saveTreeOnFile("../data/test.txt");

            // ------------------------------------
            // Add root node to the list of nodes involved in the tree-rearrangement
            updatedNodesWithinPath.push_back(_thread__topology->rootnode->getVnode_id());
            listNodesWithinPath.push_back(_thread__topology->rootnode->getVnode_id());

            // ------------------------------------
            // Recompute the likelihood
            if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
                moveLogLK = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(
                        updatedNodesWithinPath, (*_thread__topology), thread_id);
            } else {
                moveLogLK = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(
                        updatedNodesWithinPath, (*_thread__topology));
            }


            LOG_IF(FATAL, std::isinf(moveLogLK)) << "llk[Move] value is -inf for [MOVE " << candidateMoves->getMove(i)->getUID() << "]" <<
                                                 debugStackTraceMove(candidateMoves->getMove(i), _thread__topology,
                                                                     listNodesWithinPath,
                                                                     updatedNodesWithinPath,
                                                                     tshinitScore,
                                                                     moveLogLK);

            // ------------------------------------
            // Store likelihood of the move
            candidateMoves->getMove(i)->setScore(moveLogLK);

            // ------------------------------------
            // Display status of the rearrangement (deprecated)
            //candidateMoves->displayRearrangmentStatus(i, true);

            // ------------------------------------
            // Revert the move, and return to the original tree
            //candidateMoves->revertMove(i);
            //DVLOG(1) << "[TSH Cycle - Topology] [R] MOVE#" << candidateMoves->getMove(i)->getUID() << " | " << utree_->printTreeNewick(true);
//
//                // ------------------------------------
//                // Recompute the likelihood after reverting the move
//                double moveLogLK_return = -std::numeric_limits<double>::infinity();
//
//                if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
//                    moveLogLK_return = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(
//                            listNodesWithinPath);
//                } else {
//                    moveLogLK_return = dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->updateLikelihoodOnTreeRearrangement(
//                            listNodesWithinPath);
//                }

            // ------------------------------------
            // Display status of the rearrangement (deprecated)
            //candidateMoves->displayRearrangmentStatus(i, true);

//                LOG_IF(FATAL, std::isinf(moveLogLK_return))
//                << "llk[Return] value is -inf for [MOVE " << candidateMoves->getMove(i)->getUID() << "]" <<
//                debugStackTraceMove(candidateMoves->getMove(i),
//                                    utree_,
//                                    listNodesWithinPath,
//                                    updatedNodesWithinPath,
//                                    tshinitScore,
//                                    moveLogLK,
//                                    moveLogLK_return);
//
//                LOG_IF(FATAL, !ComparisonUtils::areLogicallyEqual(moveLogLK_return, tshinitScore))
//                << "Error in evaluating likelihood [MOVE " << candidateMoves->getMove(i)->getUID() << "]" <<
//                debugStackTraceMove(candidateMoves->getMove(i),
//                                    utree_,
//                                    listNodesWithinPath,
//                                    updatedNodesWithinPath,
//                                    tshinitScore,
//                                    moveLogLK,
//                                    moveLogLK_return);



        }

        // ------------------------------------
        // Count moves performed
        _count__moves_cycle[i] = 1;
        ApplicationTools::displayGauge(VectorTools::sum(_count__moves_cycle), candidateMoves->getNumberOfMoves(), '>',
                                       std::string("node " + _thread__topology->getNode(currentMove->getSourceNode())->getNodeName()));


        // ------------------------------------
        // Delete candidate topology
        delete _thread__topology;
    }
#ifdef USE_INTELTBB
//	});
#endif

    // ------------------------------------
    // Store count of moves tested for this cycle
    performed_moves.push_back(VectorTools::sum(_count__moves_cycle));

    // ------------------------------------
    // Deallocate memory
    deallocateTemporaryLikelihoodData(threads_num);

    // ------------------------------------
    ApplicationTools::displayMessage("");

}

void tshlib::TreeSearch::allocateTemporaryLikelihoodData(int numThreads) {

    for (int i = 0; i < numThreads; i++) {
        if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
            dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->addTestLikelihoodData(i);
        } else {
            dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->addTestLikelihoodData(i);
        }
    }

}

void tshlib::TreeSearch::deallocateTemporaryLikelihoodData(int numThreads) {

    for (int i = 0; i < numThreads; i++) {
        if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
            dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->removeTestLikelihoodData(i);
        } else {
            dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->removeTestLikelihoodData(i);
        }
    }

}


double tshlib::TreeSearch::iterate() {


    // Condition handler
    double c = std::abs(tshinitScore);

    while (toleranceValue < c) {

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // ------------------------------------
        // 1. Define moves according to tree-search criteria
        tshlib::TreeRearrangment *_move__set = defineMoves();

        std::ostringstream _task;
        _task << "Tree-search cycle #" << std::setfill('0') << std::setw(3) << performed_cycles + 1
              << "| Testing " << _move__set->getNumberOfMoves() << " tree rearrangements.";
        ApplicationTools::displayMessage(_task.str());

        DVLOG(1) << "[TSH Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1 << " ] " <<
                 "Initial Log-likelihood = " << TextTools::toString(tshcycleScore, 15);
        DVLOG(1) << "[TSH Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1 << " ] " <<
                 "Testing " << _move__set->getNumberOfMoves() << " moves.";

        // ------------------------------------
        // 2. Test and record likelihood of each and every candidate move
        testMoves(_move__set);

        // ------------------------------------
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        DVLOG(1) << "[TSH Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1 << " ] "
                 << "Moves tested: " << performed_moves[performed_cycles] << std::endl;
        DVLOG(1) << "[TSH Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1 << " ] "
                 << "Elapsed time: " << duration << " microseconds" << std::endl;
        DVLOG(1) << "[TSH Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1 << " ] "
                 << "*** " << (double) duration / performed_moves[performed_cycles] << " microseconds/move *** " << std::endl;

        // ------------------------------------
        // 3. Select the best move in the list and store it
        Move *bestMove = _move__set->selectBestMove(tshcycleScore);


        // ------------------------------------
        // Print winning move
        if (bestMove) {

            setTreeSearchStatus(true);

            DVLOG(1) << "[TSH Cycle " << std::setfill('0') << std::setw(3) << performed_cycles + 1 << " ] "
                     << "Found best move | "
                     << "lk = " << std::setprecision(12) << bestMove->getScore() << " | "
                     << bestMove->getClass() << "." << std::setfill('0') << std::setw(3) << bestMove->getUID()
                     << " [" << utree_->getNode(bestMove->getSourceNode())->getNodeName()
                     << " -> " << utree_->getNode(bestMove->getTargetNode())->getNodeName() << "]";


            _task.str(std::string());
            _task << "- Found new ML tree   | ["
                  << bestMove->getClass() << "." << std::setfill('0') << std::setw(3) << bestMove->getUID()
                  << "] New likelihood: " << std::setprecision(12) << bestMove->getScore();
            ApplicationTools::displayMessage(_task.str());


            std::vector<int> listNodesWithinPath, updatedNodesWithinPath;
            listNodesWithinPath = utree_->computePathBetweenNodes(utree_->getNode(bestMove->getSourceNode()),
                                                                  utree_->getNode(bestMove->getTargetNode()));
            updatedNodesWithinPath = _move__set->updatePathBetweenNodes(bestMove->getUID(), listNodesWithinPath);
            updatedNodesWithinPath.push_back(utree_->rootnode->getVnode_id());

            // ------------------------------------
            // Commit final move on the topology
            _move__set->commitMove(bestMove->getUID(), (*utree_));

            DVLOG(1) << "utree after commit " << utree_->printTreeNewick(true);

            ApplicationTools::displayTask("- Optimising " + TextTools::toString(updatedNodesWithinPath.size()) + " branches");

            if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)) {
                dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(likelihoodFunc)->topologyChangeSuccessful(updatedNodesWithinPath);
            } else {
                dynamic_cast<UnifiedTSHomogeneousTreeLikelihood *>(likelihoodFunc)->topologyChangeSuccessful(updatedNodesWithinPath);
            }

            ApplicationTools::displayTaskDone();
            ApplicationTools::displayMessage("");

            DVLOG(1) << "bpp after commit " << OutputUtils::TreeTools::writeTree2String(likelihoodFunc->getTree().clone());

            tshcycleScore = -likelihoodFunc->getValue();

            // Update tolerance
            c = std::abs(tshinitScore) - std::abs(tshcycleScore);

            // Update likelihood
            tshinitScore = tshcycleScore;

            // ------------------------------------
            // Clean memory
            delete _move__set;

        } else {

            DLOG(INFO) << "[TSH Cycle] No further likelihood improvement after " << performed_cycles << " cycles and "
                       << VectorTools::sum(performed_moves) << " performed moves. Exit loop.";

            // ------------------------------------
            // Clean memory
            delete _move__set;

            break;
        }

        // Update number of cycles
        performed_cycles = performed_cycles + 1;

        if (performed_cycles == maxTSCycles) {
            DLOG(INFO) << "[TSH Cycle] Reached max number of tree-search cycles after " << performed_cycles << " cycles";
            break;
        }

    }

    return -tshcycleScore;
}


std::string tshlib::TreeSearch::debugStackTraceMove(Move *move, Utree *_utree,
                                                    std::vector<int> listNodesInvolved,
                                                    std::vector<int> updatedList,
                                                    double initLK,
                                                    double moveLK) {

    std::ostringstream nodepath_lni;
    for (auto &node:listNodesInvolved) {
        nodepath_lni << _utree->getNode(node)->getNodeName() << ">";
    }
    std::ostringstream nodepath_uni;
    for (auto &node:updatedList) {
        nodepath_uni << _utree->getNode(node)->getNodeName() << ">";
    }


    std::ostringstream stm;
    stm << std::endl << "*** Stack trace [MOVE " << move->getClass() << "."
        << move->getUID() << "] ("
        << _utree->getNode(move->getSourceNode())->getNodeName()
        << " -> " << _utree->getNode(move->getTargetNode())->getNodeName() << ") - "
        << move->getDirection() << " ***"
        << std::endl;

    stm << "    @        [ReferenNodeList]  " << nodepath_lni.str() << std::endl;
    stm << "    @        [UpdatedNodeList]  " << nodepath_uni.str() << std::endl;
    stm << "    @        [LLK.Initial]      " << TextTools::toString(initLK, 25) << std::endl;
    stm << "    @        [LLK.Move]         " << TextTools::toString(moveLK, 25) << std::endl;


    return stm.str();

}
