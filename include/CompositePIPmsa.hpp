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
 * @file CompositePIPmsa.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 07 08 2018
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

#ifndef MINIJATI_COMPOSITEPIPMSA_HPP
#define MINIJATI_COMPOSITEPIPMSA_HPP

#include <string>
#include <vector>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include "PIPmsa.hpp"

namespace bpp {

    //*******************************************************************************************
    class iPIPmsa {

    public:

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        iPIPmsa() {}; // constructor

        virtual ~iPIPmsa() {}; // destructor

        virtual void add(PIPmsa *msa) {}; // virtual function to add a new PIPmsa

        virtual PIPmsa *getMSA(int idx=0) {};

    };

    //*******************************************************************************************
    class PIPmsaSingle : public iPIPmsa {

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        PIPmsa *pipmsa; // PIPmsaSingle has only one MSA per PIPnode

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPmsaSingle() {}; // constructor

        void add(PIPmsa *x){ pipmsa = x; }; // add (associate) an MSA

        PIPmsa *getMSA(int idx=0) {return pipmsa;}

        ~PIPmsaSingle() { // destructor
            delete pipmsa;
        }
    };

    //*******************************************************************************************
    class PIPmsaComp : public iPIPmsa {

    private:

        //***************************************************************************************
        // PRIVATE FIELDS
        //***************************************************************************************

        //***************************************************************************************
        // PRIVATE METHODS
        //***************************************************************************************

    public:

        //***************************************************************************************
        // PUBLIC FIELDS
        //***************************************************************************************

        std::vector<PIPmsa *> pipmsa; // PIPmsaComp has one or more MSAs per PIPnode

        //***************************************************************************************
        // PUBLIC METHODS
        //***************************************************************************************

        PIPmsaComp(int size) { pipmsa.resize(size); } // constructor

        void add(PIPmsa *x){ pipmsa.push_back(x); }; // add a PIPmsa to the vector (add a new alignment object)

        PIPmsa *getMSA( int idx=0) {return pipmsa.at(idx);}

        ~PIPmsaComp() { // destructor
            for (std::vector<PIPmsa *>::const_iterator iter = pipmsa.begin(); iter != pipmsa.end(); ++iter) {
                delete *iter;
            }
        }
    };
    //*******************************************************************************************

}

#endif //MINIJATI_COMPOSITEPIPMSA_HPP
