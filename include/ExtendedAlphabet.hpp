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
 * @file ExtendedAlphabet.hpp
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
#ifndef CASTOR_EXTENDEDALPHABET_HPP
#define CASTOR_EXTENDEDALPHABET_HPP


#include <Bpp/Seq/Alphabet/LetterAlphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabetState.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabetState.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Alphabet/WordAlphabet.h>


namespace bpp {

    /**
     * @brief This alphabet is used to deal with DNA_EXTENDED sequences.
     *
     * It supports all 4 nucleotides (A, T, G and C) with their standard denomination.
     * Gaps are coded by '-', unresolved characters are coded by 'X, N, O, 0 or ?'.
     * Extensive support for generic characters (e.g. 'P', 'Y', etc.) is provided.
     */
    class DNA_EXTENDED : public NucleicAlphabet {
    public:
        /**
         * @param exclamationMarkCountsAsGap If yes, '!' characters are replaced by gaps.
         * Otherwise, they are counted as unknown characters.
         */
        DNA_EXTENDED(bool exclamationMarkCountsAsGap = false);

        DNA_EXTENDED(const DNA_EXTENDED &bia) : NucleicAlphabet(bia) {}

        DNA_EXTENDED &operator=(const DNA_EXTENDED &bia) {
            NucleicAlphabet::operator=(bia);
            return *this;
        }

        DNA_EXTENDED *clone() const {
            return new DNA_EXTENDED(*this);
        }

        virtual ~DNA_EXTENDED() {}

    public:
        std::vector<int> getAlias(int state) const throw(BadIntException);

        std::vector<std::string> getAlias(const std::string &state) const throw(BadCharException);

        int getGeneric(const std::vector<int> &states) const throw(BadIntException);

        std::string getGeneric(const std::vector<std::string> &states) const throw(BadCharException);

        std::string getAlphabetType() const { return "DNA_EXTENDED"; }

        // return 4 : A, C, G, T (or U)
        unsigned int getSize() const { return 5; }

        // return 15 : gap isn't included, generic unresolved bases (N, X, ?, O, 0) count for one
        unsigned int getNumberOfTypes() const { return 16; }

        int getUnknownCharacterCode() const { return 15; }

        bool isUnresolved(int state) const { return state > 4; }

        bool isUnresolved(const std::string &state) const { return charToInt(state) > 4; }

        int getGapCharacterCode() const { return 4; }
    };


    /**
     * @brief This alphabet is used to deal with proteins.
     *
     * It supports all 20 amino-acids with their standard denomination.
     * Gaps are coded by '-', unresolved characters are coded by 'X'.
     */

    class ProteicAlphabet_Extended : public ProteicAlphabet {
        /**
         * @name Overloaded methods from AbstractAlphabet
         * @{
         */
    public:
        const ProteicAlphabetState &getState(const std::string &letter) const
        throw(BadCharException) {
            return dynamic_cast<const ProteicAlphabetState &>(
                    AbstractAlphabet::getState(letter)
            );
        }

        const ProteicAlphabetState &getState(int num) const
        throw(BadIntException) {
            return dynamic_cast<const ProteicAlphabetState &>(
                    AbstractAlphabet::getState(num)
            );
        }

    protected:

        const ProteicAlphabetState &getStateAt(size_t pos) const
        throw(IndexOutOfBoundsException) {
            return dynamic_cast<const ProteicAlphabetState &>(
                    AbstractAlphabet::getStateAt(pos)
            );
        }

        ProteicAlphabetState &getStateAt(size_t pos)
        throw(IndexOutOfBoundsException) {
            return dynamic_cast<ProteicAlphabetState &>(
                    AbstractAlphabet::getStateAt(pos)
            );
        }

        /** @} */
    public:
        ProteicAlphabet_Extended();

        // ProteicAlphabet_Extended(const ProteicAlphabet_Extended &bia) : ProteicAlphabet(bia)  {}

        ProteicAlphabet_Extended &operator=(const ProteicAlphabet_Extended &bia) {
            LetterAlphabet::operator=(bia);
            return *this;
        }

        ProteicAlphabet_Extended *clone() const {
            return new ProteicAlphabet_Extended(*this);
        }


        virtual ~ProteicAlphabet_Extended() {}


    public:
        unsigned int getSize() const { return 21; }

        unsigned int getNumberOfTypes() const { return 24; }

        int getUnknownCharacterCode() const { return 23; }

        std::vector<int> getAlias(int state) const throw(BadIntException);

        std::vector<std::string> getAlias(const std::string &state) const throw(BadCharException);

        int getGeneric(const std::vector<int> &states) const throw(BadIntException);

        int getGapCharacterCode() const { return 20; }

        std::string getGeneric(const std::vector<std::string> &states) const throw(BadCharException);

        bool isUnresolved(int state) const { return state > 21; }

        bool isUnresolved(const std::string &state) const { return charToInt(state) > 21; }

        std::string getAlphabetType() const { return "Proteic"; }

    public:

        /**
         * @name Specific methods
         *
         * @{
         */

        /**
         * @brief Get the abbreviation (3 letter code) for a state coded as char.
         *
         * @param aa Char description of the amino-acid to analyse.
         */
        std::string getAbbr(const std::string &aa) const throw(AlphabetException);

        /**
         * @brief Get the abbreviation (3 letter code) for a state coded as int.
         *
         * @param aa Int description of the amino-acid to analyse.
         */
        std::string getAbbr(int aa) const throw(AlphabetException);
        /** @} */

    };


    /**
 * @brief Codon alphabet class.
 * @author Laurent Guéguen, Julien Dutheil
 * 
 * Since codons are made of 3 nucleic bases (RNA or DNA), this class
 * has a NucleicAlphabet field used to check char description. This
 * nucleic alphabet is passed to the constructor. This class also adds
 * some methods specific to codon manipulation.
 */

    class CodonAlphabet_Extended :
            public virtual CoreWordAlphabet,
            public AbstractAlphabet {
    protected:
        const NucleicAlphabet *nAlph_;

    public: // Constructor and destructor.

        /**
         * @brief Builds a new codon alphabet from a nucleic alphabet.
         * 
         * @param alpha The nucleic alphabet to be used.
         */
        CodonAlphabet_Extended(const NucleicAlphabet *alpha) :
                AbstractAlphabet(),
                nAlph_(alpha) {
            build_();
        }

        CodonAlphabet_Extended(const CodonAlphabet_Extended &bia) :
                AbstractAlphabet(bia),
                nAlph_(bia.nAlph_) {}

        CodonAlphabet_Extended &operator=(const CodonAlphabet_Extended &bia) {
            AbstractAlphabet::operator=(bia);
            nAlph_ = bia.nAlph_;

            return *this;
        }

        CodonAlphabet_Extended *clone() const {
            return new CodonAlphabet_Extended(*this);
        }

        virtual ~CodonAlphabet_Extended() {}

        std::string getAlphabetType() const {
            return "Codon(letter=" + nAlph_->getAlphabetType() + ")";
        }

    private:
        /**
         * @name Inner utilitary functions
         *
         * @{
         */
        bool containsUnresolved(const std::string &state) const;

        bool containsGap(const std::string &state) const;

        void build_();

        /** @} */

    public:

        /**
         * @name From AbstractAlphabet
         *
         * @{
         */

        unsigned int getNumberOfTypes() const { return 66; }

        unsigned int getSize() const {
            return 65;
        }

        int getGapCharacterCode() const { return 64; }

        int getUnknownCharacterCode() const {
            return 65;
        }

        bool isUnresolved(int state) const {
            return state >= 65;
        }

        bool isUnresolved(const std::string &state) const {
            return isUnresolved(charToInt(state));
        }

        std::vector<int> getAlias(int state) const throw(BadIntException);

        std::vector<std::string> getAlias(const std::string &state) const throw(BadCharException);

        int getGeneric(const std::vector<int> &states) const throw(BadIntException) {
            return states[0];
        }

        std::string getGeneric(const std::vector<std::string> &states) const throw(BadCharException) {
            return states[0];
        }

        int charToInt(const std::string &state) const throw(BadCharException) {
            if (state.size() != 3)
                throw BadCharException(state, "CodonAlphabet_Extended::charToInt", this);
            if (containsUnresolved(state))
                return static_cast<int>(getSize());
            if (containsGap(state))
                return -1;
            else return AbstractAlphabet::charToInt(state);
        }

        /**
         * @name Codon specific methods
         *
         * @{
         */

        /**
         * @brief Get the int code for a codon given the int code of the three underlying positions.
         *
         * The int code of each position must match the nucleic alphabet specified for this alphabet.
         * @param pos1 Int description for position 1.
         * @param pos2 Int description for position 2.
         * @param pos3 Int description for position 3.
         * @return The int code of the codon.
         */

        int getCodon(int pos1, int pos2, int pos3) const {
            return (nAlph_->isUnresolved(pos1)
                    || nAlph_->isUnresolved(pos2)
                    || nAlph_->isUnresolved(pos3)) ? getUnknownCharacterCode()
                                                   : pos3 + 4 * pos2 + 16 * pos1;
        }

        /**
         * @brief Get the char code for a codon given the char code of the
         * three underlying positions.
         *
         * The char code of each position must match the nucleic alphabet
         * specified for this alphabet.
         *
         * NB: This performs pos1 + pos2 + pos3 after checking for each
         * position validity.
         *
         * @param pos1 Char description for position 1.
         * @param pos2 Char description for position 2.
         * @param pos3 Char description for position 3.
         * @return The Char code of the codon.
         */

        std::string getCodon(const std::string &pos1, const std::string &pos2, const std::string &pos3) const {
            return pos1 + pos2 + pos3;
        }

        /**
         * @brief Get the int code of the first position of a codon given its int description.
         * 
         * @param codon The int description of the codon.
         * @return The int description of the first position of the codon.
         */

        int getFirstPosition(int codon) const {
            return isUnresolved(codon) ? nAlph_->charToInt("N") : codon / 16;
        }

        /**
         * @brief Get the int code of the second position of a codon given its int description.
         * 
         * @param codon The int description of the codon.
         * @return The int description of the second position of the codon.
         */

        int getSecondPosition(int codon) const {
            return isUnresolved(codon) ? nAlph_->charToInt("N") : (codon / 4) % 4;
        }


        /**
         * @brief Get the int code of the third position of a codon given its int description.
         * 
         * @param codon The int description of the codon.
         * @return The int description of the third position of the codon.
         */

        int getThirdPosition(int codon) const {
            return isUnresolved(codon) ? nAlph_->charToInt("N") : codon % 4;
        }

        /**
         * @brief Get the char code of the first position of a codon given its char description.
         * 
         * @param codon The char description of the codon.
         * @return The char description of the first position of the codon.
         */

        std::string getFirstPosition(const std::string &codon) const {
            return codon.substr(0, 1);
        }


        /**
         * @brief Get the char code of the second position of a codon given its char description.
         * 
         * @param codon The char description of the codon.
         * @return The char description of the second position of the codon.
         */

        std::string getSecondPosition(const std::string &codon) const {
            return codon.substr(1, 1);
        }


        /**
         * @brief Get the char code of the third position of a codon given its char description.
         * 
         * @param codon The char description of the codon.
         * @return The char description of the third position of the codon.
         */

        std::string getThirdPosition(const std::string &codon) const {
            return codon.substr(2, 1);
        }


        /**
         * @name From CoreWordAlphabet
         *
         * @{
         */

        unsigned int getLength() const {
            return 3;
        }

        bool hasUniqueAlphabet() const {
            return true;
        }

        const Alphabet *getNAlphabet(size_t n) const {
            return nAlph_;
        }

        int getWord(const Sequence &seq, size_t pos = 0) const {
            if (seq.size() < pos + 3)
                throw IndexOutOfBoundsException("CodonAlphabet_Extended::getWord", pos, 0, seq.size() - 3);
            return getCodon(seq[pos], seq[pos + 1], seq[pos + 2]);
        }

        /**
         * @brief Get the char code for a word given the char code of the
         * underlying positions.
         *
         * The char code of each position must match the corresponding alphabet specified at this position.
         * @param vpos vector description for all the positions.
         * @param pos the start position to match in the vector.
         * @return The string of the word.
         * @throw IndexOutOfBoundsException In case of wrong position.
         */

        std::string getWord(const std::vector<std::string> &vpos, size_t pos = 0) const {
            if (vpos.size() < pos + 3)
                throw IndexOutOfBoundsException("CodonAlphabet_Extended::getWord", pos, 0, vpos.size() - 3);

            return getCodon(vpos[pos], vpos[pos + 1], vpos[pos + 2]);
        }

        int getWord(const std::vector<int> &vpos, size_t pos = 0) const {
            if (vpos.size() < pos + 3)
                throw IndexOutOfBoundsException("CodonAlphabet_Extended::getWord", pos, 0, vpos.size() - 3);

            return getCodon(vpos[pos], vpos[pos + 1], vpos[pos + 2]);
        }


        int getNPosition(int codon, size_t pos) const {
            if (isUnresolved(codon))
                return nAlph_->getUnknownCharacterCode();
            else
                return (pos == 0 ? codon / 16 :
                        (pos == 1 ? (codon / 4) % 4
                                  : codon % 4));
        }

        /**
         * @brief Get the int codes of each position of a word given its int description.
         *
         * @param word The int description of the word.
         * @return The int description of the positions of the codon.
         */

        std::vector<int> getPositions(int word) const {
            if (isUnresolved(word)) {
                int n = nAlph_->getUnknownCharacterCode();
                return std::vector<int>{n, n, n};
            } else
                return std::vector<int>{word / 16, (word / 4) % 4, word % 4};
        }


        /**
         * @brief Get the char code of the Nth position of a codon given its char description.
         * 
         * @param codon The char description of the codon.
         * @param pos the position in the codon (starting at 0)
         * @return The char description of the position of the codon.
         */

        std::string getNPosition(const std::string &codon, size_t pos) const {
            return codon.substr(pos, 1);
        }

        /**
         * @brief Get the char codes of each position of a word given its char description.
         *
         * @param word The char description of the word.
         * @return The char description of the three positions of the word.
         */

        std::vector<std::string> getPositions(const std::string &word) const {
            return std::vector<std::string>{word.substr(0, 1), word.substr(1, 1), word.substr(2, 1)};
        }

        /**
         * @brief Translate a whole sequence from letters alphabet to words alphabet.
         *
         * @param sequence A sequence in letters alphabet.
         * @param pos the start postion (default 0)
         * @return The corresponding sequence in words alphabet.
         * @throw AlphabetMismatchException If the sequence alphabet do not match the source alphabet.
         * @throw Exception                 Other kind of error, depending on the implementation.
         */

        Sequence *translate(const Sequence &sequence, size_t = 0) const;

        /**
         * @brief Translate a whole sequence from words alphabet to letters alphabet.
         *
         * @param sequence A sequence in words alphabet.
         * @return The corresponding sequence in letters alphabet.
         * @throw AlphabetMismatchException If the sequence alphabet do not match the target alphabet.
         * @throw Exception                 Other kind of error, depending on the implementation.
         */

        Sequence *reverse(const Sequence &sequence) const;

        /*
         *
         * @}
         */

        /**
         * @brief Get the number of G+C in codon
         * @param codon The int description of the codon.
         *
         * @return The number of G+C in codon
         */

        int getGCinCodon(int codon) const;

        /**
         * @return The nucleic alphabet associated to this codon alphabet.
         */

        const NucleicAlphabet *getNucleicAlphabet() const {
            return nAlph_;
        }

        /**
         * @name Overloaded AbstractAlphabet methods.
         * @{
         */

        unsigned int getStateCodingSize() const { return 3; }
        /** @} */
    };

}

#endif //CASTOR_EXTENDEDALPHABET_HPP
