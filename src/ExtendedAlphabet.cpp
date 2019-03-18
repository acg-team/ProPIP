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
 * @file ExtendedAlphabet.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 11 01 2018
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

#include <Bpp/Text/TextTools.h>
#include <Bpp/Utils/MapTools.h>

using namespace bpp;

// From STL:

using namespace std;

#include "ExtendedAlphabet.hpp"


DNA_EXTENDED::DNA_EXTENDED(bool exclamationMarkCountsAsGap) {
    // Alphabet content definition
    // all unresolved bases use n°15
    registerState(new NucleicAlphabetState(0, "A", 1, "Adenine"));
    registerState(new NucleicAlphabetState(1, "C", 2, "Cytosine"));
    registerState(new NucleicAlphabetState(2, "G", 4, "Guanine"));
    registerState(new NucleicAlphabetState(3, "T", 8, "Thymine"));
    registerState(new NucleicAlphabetState(4, "-", 16, "Gap"));
    registerState(new NucleicAlphabetState(5, "M", 3, "Adenine or Cytosine"));
    registerState(new NucleicAlphabetState(6, "R", 5, "Purine (Adenine or Guanine)"));
    registerState(new NucleicAlphabetState(7, "W", 9, "Adenine or Thymine"));
    registerState(new NucleicAlphabetState(8, "S", 6, "Cytosine or Guanine"));
    registerState(new NucleicAlphabetState(9, "Y", 10, "Pyrimidine (Cytosine or Thymine)"));
    registerState(new NucleicAlphabetState(10, "K", 12, "Guanine or Thymine"));
    registerState(new NucleicAlphabetState(11, "V", 7, "Adenine or Cytosine or Guanine"));
    registerState(new NucleicAlphabetState(12, "H", 11, "Adenine or Cytosine or Thymine"));
    registerState(new NucleicAlphabetState(13, "D", 13, "Adenine or Guanine or Thymine"));
    registerState(new NucleicAlphabetState(14, "B", 14, "Cytosine or Guanine or Thymine"));
    registerState(new NucleicAlphabetState(15, "N", 15, "Unresolved base"));
    registerState(new NucleicAlphabetState(15, "X", 15, "Unresolved base"));
    registerState(new NucleicAlphabetState(15, "O", 15, "Unresolved base"));
    registerState(new NucleicAlphabetState(15, "0", 15, "Unresolved base"));
    registerState(new NucleicAlphabetState(15, "?", 15, "Unresolved base"));

    if (exclamationMarkCountsAsGap)
        registerState(new NucleicAlphabetState(-1, "!", 0, "Frameshift"));
    else
        registerState(new NucleicAlphabetState(15, "!", 15, "Unresolved base"));

}

std::vector<int> DNA_EXTENDED::getAlias(int state) const throw(BadIntException) {
    if (!isIntInAlphabet(state))
        throw BadIntException(state, "DNA_EXTENDED::getAlias(int): Specified base unknown.");
    vector<int> v;
    const NucleicAlphabetState &st = getState(state);
    if (state == -1)
        v.push_back(-1);
    if (st.getBinaryCode() & 1)
        v.push_back(0);
    if (st.getBinaryCode() & 2)
        v.push_back(1);
    if (st.getBinaryCode() & 4)
        v.push_back(2);
    if (st.getBinaryCode() & 8)
        v.push_back(3);
    if (st.getBinaryCode() & 16)
        v.push_back(4);


    return v;
}

std::vector<std::string> DNA_EXTENDED::getAlias(const std::string &state) const throw(BadCharException) {
    string locstate = TextTools::toUpper(state);
    if (!isCharInAlphabet(locstate)) throw BadCharException(locstate, "DNA_EXTENDED::getAlias(int): Specified base unknown.");
    vector<int> vi = this->getAlias(this->charToInt(state));
    vector<string> v;
    for (unsigned int i = 0; i < vi.size(); i++)
        v.push_back(this->intToChar(vi[i]));
    return v;
}


int DNA_EXTENDED::getGeneric(const std::vector<int> &states) const throw(BadIntException) {
    int v = 0;
    for (size_t i = 0; i < states.size(); ++i) {
        if (!isIntInAlphabet(states[i]))
            throw BadIntException(states[i], "DNA_EXTENDED::getGeneric(const vector<int>& states): Specified base unknown.");
        v |= getState(states[i]).getBinaryCode();
    }
    return getStateByBinCode(v).getNum();
}


std::string DNA_EXTENDED::getGeneric(const std::vector<std::string> &states) const throw(BadCharException) {
    vector<int> vi;
    for (unsigned int i = 0; i < states.size(); ++i) {
        if (!isCharInAlphabet(states[i]))
            throw BadCharException(states[i], "DNA_EXTENDED::getGeneric(const vector<string>& states): Specified base unknown.");
        vi.push_back(this->charToInt(states[i]));
    }
    return intToChar(getGeneric(vi));
}


ProteicAlphabet_Extended::ProteicAlphabet_Extended() {
    // Alphabet content definition
    resize(0);
    remap();
    // Alphabet content definition
    registerState(new ProteicAlphabetState(0, "A", "ALA", "Alanine"));
    registerState(new ProteicAlphabetState(1, "R", "ARG", "Arginine"));
    registerState(new ProteicAlphabetState(2, "N", "ASN", "Asparagine"));
    registerState(new ProteicAlphabetState(3, "D", "ASP", "Asparatic Acid"));
    registerState(new ProteicAlphabetState(4, "C", "CYS", "Cysteine"));
    registerState(new ProteicAlphabetState(5, "Q", "GLN", "Glutamine"));
    registerState(new ProteicAlphabetState(6, "E", "GLU", "Glutamic acid"));
    registerState(new ProteicAlphabetState(7, "G", "GLY", "Glycine"));
    registerState(new ProteicAlphabetState(8, "H", "HIS", "Histidine"));
    registerState(new ProteicAlphabetState(9, "I", "ILE", "Isoleucine"));
    registerState(new ProteicAlphabetState(10, "L", "LEU", "Leucine"));
    registerState(new ProteicAlphabetState(11, "K", "LYS", "Lysine"));
    registerState(new ProteicAlphabetState(12, "M", "MET", "Methionine"));
    registerState(new ProteicAlphabetState(13, "F", "PHE", "Phenylalanine"));
    registerState(new ProteicAlphabetState(14, "P", "PRO", "Proline"));
    registerState(new ProteicAlphabetState(15, "S", "SER", "Serine"));
    registerState(new ProteicAlphabetState(16, "T", "THR", "Threonine"));
    registerState(new ProteicAlphabetState(17, "W", "TRP", "Tryptophan"));
    registerState(new ProteicAlphabetState(18, "Y", "TYR", "Tyrosine"));
    registerState(new ProteicAlphabetState(19, "V", "VAL", "Valine"));
    registerState(new ProteicAlphabetState(20, "-", "GAP", "Gap"));
    registerState(new ProteicAlphabetState(21, "B", "B", "N or D"));
    registerState(new ProteicAlphabetState(22, "Z", "Z", "Q or E"));
    registerState(new ProteicAlphabetState(23, "X", "X", "Unresolved amino acid"));
    registerState(new ProteicAlphabetState(23, "O", "O", "Unresolved amino acid"));
    registerState(new ProteicAlphabetState(23, "0", "0", "Unresolved amino acid"));
    registerState(new ProteicAlphabetState(23, "?", "?", "Unresolved amino acid"));
    registerState(new ProteicAlphabetState(-1, "*", "STOP", "Stop"));

}


string ProteicAlphabet_Extended::getAbbr(const string &aa) const throw(AlphabetException) {
    string AA = TextTools::toUpper(aa);
    return getState(aa).getAbbreviation();
}


string ProteicAlphabet_Extended::getAbbr(int aa) const throw(AlphabetException) {
    return getState(aa).getAbbreviation();
}


vector<int> ProteicAlphabet_Extended::getAlias(int state) const throw(BadIntException) {
    if (!isIntInAlphabet(state))
        throw BadIntException(state, "ProteicAlphabet_Extended::getAlias(int): Specified base unknown.");
    vector<int> v;
    if (state == 21) {          // N or D
        v.resize(2);
        v[0] = 2;
        v[1] = 3;
    } else if (state == 22) {   // Q or E
        v.resize(2);
        v[0] = 5;
        v[1] = 6;
    } else if (state == 23) {   // all!
        v.resize(21);
        for (size_t i = 0; i < 21; i++) {
            v[i] = static_cast<int>(i);
        }
    } else {
        v.resize(1);
        v[0] = state;
    }
    return v;
}


vector<string> ProteicAlphabet_Extended::getAlias(const string &state) const throw(BadCharException) {
    string locstate = TextTools::toUpper(state);
    if (!isCharInAlphabet(locstate))
        throw BadCharException(locstate, "ProteicAlphabet_Extended::getAlias(int): Specified base unknown.");
    vector<string> v;

    if (locstate == "B") {          // N or D
        v.resize(2);
        v[0] = "N";
        v[1] = "D";
    } else if (locstate == "Z") {   // Q or E
        v.resize(2);
        v[0] = "Q";
        v[1] = "E";
    } else if (locstate == "X"
               || locstate == "O"
               || locstate == "0"
               || locstate == "?") {   // all!
        v.resize(21);
        for (int i = 0; i < 21; i++) {
            v[static_cast<size_t>(i)] = getState(i).getLetter();
        }
    } else {
        v.resize(1);
        v[0] = locstate;
    }
    return v;
}


int ProteicAlphabet_Extended::getGeneric(const vector<int> &states) const throw(BadIntException) {
    map<int, int> m;
    for (unsigned int i = 0; i < states.size(); ++i) {
        vector<int> tmp_s = this->getAlias(states[i]); // get the states for generic characters
        for (unsigned int j = 0; j < tmp_s.size(); ++j) {
            m[tmp_s[j]]++; // add each state to the list
        }
    }
    vector<int> ve = MapTools::getKeys(m);

    string key;
    for (unsigned int i = 0; i < ve.size(); ++i) {
        if (!isIntInAlphabet(ve[i]))
            throw BadIntException(ve[i], "ProteicAlphabet_Extended::getGeneric(const vector<int>): Specified base unknown.");
        key += "_" + TextTools::toString(ve[i]);
    }
    map<string, int> g;
    g["_2_3"] = 21;
    g["_5_6"] = 22;
    int v;
    map<string, int>::iterator it = g.find(key);
    if (ve.size() == 1) {
        v = ve[0];
    } else if (it != g.end()) {
        v = it->second;
    } else {
        v = 23;
    }
    return v;
}


string ProteicAlphabet_Extended::getGeneric(const vector<string> &states) const throw(BadCharException) {
    map<string, int> m;
    for (unsigned int i = 0; i < states.size(); ++i) {
        vector<string> tmp_s = this->getAlias(states[i]); // get the states for generic characters
        for (unsigned int j = 0; j < tmp_s.size(); ++j) {
            m[tmp_s[j]]++; // add each state to the list
        }
    }
    vector<string> ve = MapTools::getKeys(m);

    string key;
    for (unsigned int i = 0; i < ve.size(); ++i) {
        if (!isCharInAlphabet(ve[i]))
            throw BadCharException(ve[i], "ProteicAlphabet_Extended::getAlias(const vector<string>): Specified base unknown.");
        key += TextTools::toString(ve[i]);
    }
    map<string, string> g;
    g["DN"] = "B";
    g["EQ"] = "Z";
    string v;
    map<string, string>::iterator it = g.find(key);
    if (ve.size() == 1) {
        v = ve[0];
    } else if (it != g.end()) {
        v = it->second;
    } else {
        v = "?";
    }
    return v;
}


void CodonAlphabet_Extended::build_() {
    vector<AlphabetState *> states(66);

    //states[0] = new AlphabetState(-1, "---", "gap");

    size_t i = 0;
    for (int i1 = 0; i1 < 4; ++i1)
        for (int i2 = 0; i2 < 4; ++i2)
            for (int i3 = 0; i3 < 4; ++i3) {
                string s = nAlph_->intToChar(i1) + nAlph_->intToChar(i2) + nAlph_->intToChar(i3);
                states[i] = new AlphabetState(static_cast<int>(i), s, s);
                i++;
            }

    states[64] = new AlphabetState(static_cast<int>(64), "---", "gap");
    states[65] = new AlphabetState(static_cast<int>(65), "NNN", "Unresolved");

    //Now register all states once for all:
    for (i = 0; i < states.size(); ++i)
        registerState(states[i]);

}

int CodonAlphabet_Extended::getGCinCodon(int codon) const {
    int i = 0;
    int j = getFirstPosition(codon);
    if (j == 1 || j == 2)
        i++;
    j = getSecondPosition(codon);
    if (j == 1 || j == 2)
        i++;
    j = getThirdPosition(codon);
    if (j == 1 || j == 2)
        i++;

    return i;
}

bool CodonAlphabet_Extended::containsUnresolved(const std::string &state) const {
    if (state.length() != 3)
        throw BadCharException(state, "CodonAlphabet_Extended::containsUnresolved", this);

    for (size_t i = 0; i < 3; i++) {
        if (nAlph_->isUnresolved(state.substr(i, 1))) {
            return true;
        }
    }
    return false;
}


bool CodonAlphabet_Extended::containsGap(const std::string &state) const {
    if (state.length() != 3)
        throw BadCharException(state, "CodonAlphabet_Extended::containsGap", this);

    for (size_t i = 0; i < 3; i++) {
        if (nAlph_->isGap(state.substr(i, 1)))
            return true;
    }

    return false;
}


Sequence *CodonAlphabet_Extended::translate(const Sequence &sequence, size_t pos) const {
    vector<int> content;

    size_t s = sequence.size();
    size_t i = pos;

    while (i + 3 <= s) {
        content.push_back(getWord(sequence, i));
        i += 3;
    }

    return new BasicSequence(sequence.getName(), content, this);
}


Sequence *CodonAlphabet_Extended::reverse(const Sequence &sequence) const {
    Sequence *pseq = new BasicSequence(sequence.getName(), "", getNAlphabet(0));

    size_t s = sequence.size();
    for (size_t i = 0; i < s; i++) {
        pseq->append(getPositions(sequence[i]));
    }

    return pseq;
}


std::vector<int> CodonAlphabet_Extended::getAlias(int state) const throw(BadIntException) {
    if (!isIntInAlphabet(state))
        throw BadIntException(state, "WordAlphabet::getAlias(int): Specified base unknown.");
    vector<int> v;

    if (state == 65) {
        v.resize(65);
        for (size_t i = 0; i < 65; ++i) {
            v[i] = static_cast<int>(i);
        }
    } else {
        v.resize(1);
        v[0] = state;
    }
    return v;
}


std::vector<std::string> CodonAlphabet_Extended::getAlias(const std::string &state) const throw(BadCharException) {
    string locstate = TextTools::toUpper(state);
    if (!isCharInAlphabet(locstate))
        throw BadCharException(locstate, "CodonAlphabet_Extended::getAlias(string): Specified base unknown.");
    vector<string> v;

    if (locstate == "NNN") {
        v.resize(65);
        for (size_t i = 0; i < 65; ++i) {
            v[i] = intToChar(static_cast<int>(i));
        }
    } else {
        v.resize(1);
        v[0] = state;
    }
    return v;
}
