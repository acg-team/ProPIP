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
 * @file TreeSearchEngine.cpp
 * @author Lorenzo Gatti
 * @date 11 10 2017
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

#include <string>
#include "TreeSearch.hpp"

namespace tshlib {

// constructor of TreeSearchHL,
    TreeSearch::TreeSearch(std::string name) {
        this->test_field = name;
    }

//copy constructor for making a new copy of a TreeSearchHL
    TreeSearch::TreeSearch(const TreeSearch &copy_from) {

    }

//copy assignment for assigning a value from one TreeSearchHL to another
    TreeSearch &TreeSearch::operator=(const TreeSearch &copy_from) {
    }

// destructor, just an example
    TreeSearch::~TreeSearch() {
        //delete[] this->test_field ;
    }

    void TreeSearch::setString(std::string input = "") {
        this->test_field = input;
    }

    int TreeSearch::getLength() {
        return (int) this->test_field.length();
    }
}