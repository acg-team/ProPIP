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
        //int moveStepParentNode_;  // Pointer to the node that acts as step parent during the SPR move
        //int moveStepChildNode_;   // Pointer to the node that acts as step parent during the SPR move

        //bool onPseudoRoot_;                // If the SPR node is infact a pseudoroot


    public:
        int moveUID_;                           // Move UID - Useful in case of parallel independent executions
        //std::string moveName_;                  // Move Name - Unused
        int moveRadius_;                        // Move Radius
        MoveDirections moveDirection_;          // Move Direction for applying a rotation to the VirtualNode pointers
        double moveScore_;                      // Likelihood of the move if applied
        //bool moveApplied_;                      // Indicator is set to true if the move is applied to the tree
        std::string moveClassDescription_;      // String indicating the move class (i.e. NNI,SPR,TBR) - Needed for mixed tree-search strategies
        //MoveType moveType_;                     // Integer indicating the move class (i.e. NNI,SPR,TBR) - Needed for mixed tree-search strategies
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
            //moveApplied_ = inMove.moveApplied_;
            moveClassDescription_ = inMove.moveClassDescription_;
            moveDirection_ = inMove.moveDirection_;
            moveTargetNode_ = inMove.moveTargetNode_;
            moveSourceNode_ = inMove.moveSourceNode_;
            //moveStepParentNode_ = inMove.moveStepParentNode_;
            is_duplicate = inMove.is_duplicate;
        }

        Move &operator=(const Move &inMove) {
            moveUID_ = inMove.moveUID_;
            //moveName_ = "copy_" + inMove.moveName_;
            moveRadius_ = inMove.moveRadius_;
            moveScore_ = inMove.moveScore_;
            //moveApplied_ = inMove.moveApplied_;
            moveClassDescription_ = inMove.moveClassDescription_;
            //moveType_ = inMove.moveType_;
            moveDirection_ = inMove.moveDirection_;
            moveTargetNode_ = inMove.moveTargetNode_;
            moveSourceNode_ = inMove.moveSourceNode_;
            //moveStepParentNode_ = inMove.moveStepParentNode_;

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
        //VirtualNode *getTargetNode() {
        //    return moveTargetNode_;
        //};
        int getTargetNode() { return moveTargetNode_; };
/*!
        * @brief Returns the source node pointer
        * @return VirtualNode pointer of the source node
        */
        //VirtualNode *getSourceNode() { return moveSourceNode_; };
        int getSourceNode() { return moveSourceNode_; };

        /*!
         * @brief Set the protected move_targetnode field
         * @param target_node PhyTree Pointer to the target node
         */
        //void setTargetNode(VirtualNode *in_target_node) { moveTargetNode_ = in_target_node; };
        void setTargetNode(int in_target_node) { moveTargetNode_ = in_target_node; };

        //void setSourceNode(VirtualNode *in_source_node) { moveSourceNode_ = in_source_node; };
        void setSourceNode(int in_source_node) { moveSourceNode_ = in_source_node; };

        //void setStepParentNode(VirtualNode *in_stepParent) { moveStepParentNode_ = in_stepParent; };
        //void setStepParentNode(int in_stepParent) { moveStepParentNode_ = in_stepParent; };

        //VirtualNode *getStepParentNode() { return moveStepParentNode_; };
        //int getStepParentNode() { return moveStepParentNode_; };

        //void setStepChildNode(VirtualNode *in_stepChild) { moveStepChildNode_ = in_stepChild; };
        //void setStepChildNode(int in_stepChild) { moveStepChildNode_ = in_stepChild; };

        //VirtualNode *getStepChildNode() { return moveStepChildNode_; };
        //int getStepChildNode() { return moveStepChildNode_; };

        //MoveType getType() const { return moveType_; };

        void setRadius(int in_radius) { Move::moveRadius_ = in_radius; };

        int getRadius() const { return moveRadius_; }

        void setDirection(MoveDirections in_direction) { Move::moveDirection_ = in_direction; };

        int getUID() const { return moveUID_; }

        void setUID(int in_moveUID_) { Move::moveUID_ = in_moveUID_; }

        //const std::string &getName() const { return moveName_; }

        //void setName(const std::string &in_moveName_) { Move::moveName_ = in_moveName_; }

        double getScore() const { return moveScore_; }

        void setScore(double in_moveScore_) { Move::moveScore_ = in_moveScore_; }

        //bool isOverPseudoRoot_() const { return onPseudoRoot_; }

        //void setOnPseudoRoot_(bool in_onPseudoRoot_) { Move::onPseudoRoot_ = in_onPseudoRoot_; }

        MoveDirections getMoveDirection() const { return moveDirection_; }

        TreeSearchHeuristics getMoveStrategy() const { return moveStrategy_; }

    };


}

#endif //TSHLIB_MOVE_HPP
