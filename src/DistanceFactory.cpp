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
 * @file DistanceFactory.cpp
 * @author Adam Szalkowski
 * @author Massimo Maiolo
 * @date 06 08 2020
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


#include "DistanceFactory.hpp"
#include "DistanceFactoryAlign.hpp"
#include "DistanceFactoryAngle.hpp"
#include "DistanceFactoryPrealigned.hpp"


using namespace DistanceFactoryPrographMSA;

DistanceFactory *DistanceFactory::getDefault(bpp::SubstitutionModel *smodel, bool prealigned, int ALPHABET_DIM, bool nwdist_flag) {
    DistanceFactory *dist_factory = NULL;
    if(!prealigned) {
        if(nwdist_flag) {
            dist_factory = new DistanceFactoryAlign(smodel,ALPHABET_DIM);
        } else {

            int K;
            if(ALPHABET_DIM==4){
                K=6;
            }else if(ALPHABET_DIM==20){
                K=2;
            }else{
                K=2;
            }

            dist_factory = new DistanceFactoryAngle(ALPHABET_DIM,K);
        }
    } else {
        dist_factory = new DistanceFactoryPrealigned(smodel,ALPHABET_DIM);
    }
    return dist_factory;
}
