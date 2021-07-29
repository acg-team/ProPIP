/***************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 * Copyright (C) 2015-2019 by Lorenzo Gatti & Massimo Maiolo
 ***************************************************************************
 * This file is part of Tree Search Heuristics Library (TSH-LIB)
 *
 * TSH-LIB is a C++ Library whose purpose is generate alternative tree
 * toplogies applying NNI, SPR, and F-VFNI rearrangements on an initial tree
 *
 * TSH-LIB is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * TSH-LIB is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public
 * License along with TSH-LIB. If not, see <http://www.gnu.org/licenses/>.
 ***************************************************************************/

/**
 * @file Move.hpp
 * @author Lorenzo Gatti
 * @date 11 06 2018
 * @version 3.0.2
 * @maintainer Lorenzo Gatti
 * @email lg@lorenzogatti.me
 * @status Development
 *
 * @brief
 * @details
 * @pre
 * @bug
 * @warning
 *
 * @see For more information visit: https://bitbucket.org/lorenzogatti89/tshlib
 */
#ifndef TSHLIB_MOVE_HPP
#define TSHLIB_MOVE_HPP

#include "Utilities.hpp"
#include <Utree.hpp>

namespace tshlib {
    class Move {

    private:

    protected:
        int moveTargetNode_;      // Pointer to the target node found during the node search
        int moveSourceNode_;      // Pointer to the source node

    public:
        int moveUID_;                           // Move UID - Useful in case of parallel independent executions
        int moveRadius_;                        // Move Radius
        MoveDirections moveDirection_;          // Move Direction for applying a rotation to the VirtualNode pointers
        double moveScore_;                      // Likelihood of the move if applied
        std::string moveClassDescription_;      // String indicating the move class (i.e. NNI,SPR,TBR) - Needed for mixed tree-search strategies
        TreeSearchHeuristics moveStrategy_;     // Store the strategy used to generate this candidate.
        bool is_duplicate;
        std::vector<int> node2Opt;

        /*!
         * @brief Standard constructor
         */
        Move();

        /*!
        * @brief Standard deconstructor
        */
        ~Move();

        Move(const Move &inMove) {

            moveUID_ = inMove.moveUID_;
            moveRadius_ = inMove.moveRadius_;
            moveScore_ = inMove.moveScore_;
            moveClassDescription_ = inMove.moveClassDescription_;
            moveDirection_ = inMove.moveDirection_;
            moveTargetNode_ = inMove.moveTargetNode_;
            moveSourceNode_ = inMove.moveSourceNode_;
            is_duplicate = inMove.is_duplicate;
        }

        Move &operator=(const Move &inMove) {
            moveUID_ = inMove.moveUID_;
            moveRadius_ = inMove.moveRadius_;
            moveScore_ = inMove.moveScore_;
            moveClassDescription_ = inMove.moveClassDescription_;
            moveDirection_ = inMove.moveDirection_;
            moveTargetNode_ = inMove.moveTargetNode_;
            moveSourceNode_ = inMove.moveSourceNode_;

        };

        void initMove();

        std::string getDirection() const;

        void setClass(TreeSearchHeuristics tsStrategy, bool _location__overpseudoroot);

        /*!
         * @brief Reset the protected move_targetnode field
         */
        void deleteTargetNode();

        /*!
         * @brief Returns the target node pointer
         * @return VirtualNode pointer of the target node
         */
        int getTargetNode() { return moveTargetNode_; };
/*!
        * @brief Returns the source node pointer
        * @return VirtualNode pointer of the source node
        */
        int getSourceNode() { return moveSourceNode_; };

        /*!
         * @brief Set the protected move_targetnode field
         * @param target_node PhyTree Pointer to the target node
         */
        void setTargetNode(int in_target_node) { moveTargetNode_ = in_target_node; };

        void setSourceNode(int in_source_node) { moveSourceNode_ = in_source_node; };

        void setRadius(int in_radius) { Move::moveRadius_ = in_radius; };

        int getRadius() const { return moveRadius_; }

        void setDirection(MoveDirections in_direction) { Move::moveDirection_ = in_direction; };

        int getUID() const { return moveUID_; }

        void setUID(int in_moveUID_) { Move::moveUID_ = in_moveUID_; }

        double getScore() const { return moveScore_; }

        void setScore(double in_moveScore_) { Move::moveScore_ = in_moveScore_; }

        MoveDirections getMoveDirection() const { return moveDirection_; }

        TreeSearchHeuristics getMoveStrategy() const { return moveStrategy_; }

    };


}

#endif //TSHLIB_MOVE_HPP
