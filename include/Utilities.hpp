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
 * @file Utilities.hpp
 * @author Lorenzo Gatti
 * @date 10 11 2017
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

#ifndef TSHLIB_UTILITIES_H
#define TSHLIB_UTILITIES_H
namespace tshlib {
    enum class NodeRotation {
        undef = 0, clockwise = 1, counterclockwise = 2
    };

    enum class NodePosition {
        left, right, up, undef
    };

    enum class MoveDirections {
        down_left, down_right, up, up_right, up_left, undef
    };

    enum class MoveType {
        NNI = 1, SPR = 2, TBR = 3, FNNI = 4, VFNNI = 5, undef = 0
    };

    enum class TreeSearchStopCondition {
        iterations, convergence
    };

    enum class TreeSearchHeuristics {
        swap, phyml, mixed, nosearch
    };

    enum class StartingNodeHeuristics {
        particle_swarm, hillclimbing, greedy, undef
    };

    enum class TreeRearrangmentOperations {
        classic_NNI, classic_SPR, classic_TBR, classic_Mixed
    };

}
#endif //TSHLIB_UTILITIES_H
