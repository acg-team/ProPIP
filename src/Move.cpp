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
 * @file Move.cpp
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


#include <numeric>
#include <limits>
#include <iomanip>
#include <iterator>
#include <chrono>
#include <algorithm>
#include <Utilities.hpp>


#include "Move.hpp"


namespace tshlib {
    Move::~Move() = default;

    Move::Move() = default;

    void Move::initMove() {
        moveUID_ = 0;
        //moveName_ = "undef";
        //moveType_ = MoveType::undef;
        moveClassDescription_ = "undef";
        moveDirection_ = MoveDirections::undef;

        moveScore_ = -std::numeric_limits<double>::infinity();

       // moveApplied_ = false;
    }

    std::string Move::getDirection() const {
        std::string rtToken;
        switch (moveDirection_) {
            case MoveDirections::down_left:
                rtToken = "left";
                break;
            case MoveDirections::up:
                rtToken = "up";
                break;
            case MoveDirections::up_left:
                rtToken = "up-left";
                break;
            case MoveDirections::up_right:
                rtToken = "up-right";
                break;
            case MoveDirections::down_right :
                rtToken = "right";
                break;
            case MoveDirections::undef :
                rtToken = "undef";
                break;

        }

        return rtToken;
    }

    void Move::deleteTargetNode() {

        moveTargetNode_ = 0;

    }

    /*
    void Move::setClass(TreeSearchHeuristics tsStrategy, bool _location__overpseudoroot) {

        // Set tree search strategy associated to this move
        moveStrategy_ = tsStrategy;

        // Get radius of the current move
        int radius = getRadius();

        if (radius == 3) {
            switch (tsStrategy) {
                case TreeSearchHeuristics::swap:
                case TreeSearchHeuristics::phyml:
                case TreeSearchHeuristics::mixed:
                    moveType_ = MoveType::NNI;
                    break;
                case TreeSearchHeuristics::nosearch:
                    moveType_ = MoveType::undef;
                    break;
            }
        } else if (radius == 4) {
            switch (tsStrategy) {
                case TreeSearchHeuristics::swap:
                    moveType_ = MoveType::FNNI;
                    break;
                case TreeSearchHeuristics::phyml:
                    moveType_ = MoveType::SPR;
                    break;
                case TreeSearchHeuristics::mixed:
                    moveType_ = MoveType::SPR;
                    break;
                case TreeSearchHeuristics::nosearch:
                    moveType_ = MoveType::undef;
                    break;
            }
        } else if (radius > 4) {
            switch (tsStrategy) {
                case TreeSearchHeuristics::swap:
                    moveType_ = MoveType::VFNNI;
                    break;
                case TreeSearchHeuristics::phyml:
                    moveType_ = MoveType::SPR;
                    break;
                case TreeSearchHeuristics::mixed:
                    moveType_ = MoveType::SPR;
                    break;
                case TreeSearchHeuristics::nosearch:
                    moveType_ = MoveType::undef;
                    break;
            }
        } else {

            moveType_ = MoveType::undef;
        }

        // If either the source or the target node define are on the pseudoroot, and the tree-search strategy is [PhyML], then the
        // movetype is handled as TBR
        //if ((moveTargetNode_->isPseudoRootNode() || moveSourceNode_->isPseudoRootNode()) && tsStrategy == TreeSearchHeuristics::phyml)
        if (_location__overpseudoroot && tsStrategy == TreeSearchHeuristics::phyml)
            moveType_ = MoveType::TBR;


        moveClassDescription_ = getClass();
    }
    */

    void Move::setClass(TreeSearchHeuristics tsStrategy, bool _location__overpseudoroot) {

        // Set tree search strategy associated to this move
        moveStrategy_ = tsStrategy;

        // Get radius of the current move
        int radius = getRadius();

        if (radius == 3) {
            //moveType_ = MoveType::NNI;
            moveClassDescription_ = "NNI";
        } else if(radius >=4) {
            //moveType_ = MoveType::SPR;
            moveClassDescription_ = "SPR";
        }else{
            exit(EXIT_FAILURE);
        }

        //moveClassDescription_ = getClass();
    }

}