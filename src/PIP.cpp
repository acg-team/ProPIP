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
 * @file PIP.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 20 12 2017
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
#include "PIP.hpp"
#include "Utils.hpp"
#include "LambertW.h"
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <boost/algorithm/string.hpp>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <glog/logging.h>


PIP_Nuc::PIP_Nuc(const NucleicAlphabet *alpha, SubstitutionModel *basemodel, const SequenceContainer &data, double lambda, double mu,
                 bool initFreqFromData) :
        AbstractParameterAliasable("PIP."),
        AbstractReversibleNucleotideSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "PIP."),
        lambda_(lambda), mu_(mu), r_(), l_(), k_(), exp1_(), exp2_(), p_(size_, size_) {


    // Setting basemodel to PIP
    submodel_ = basemodel;

    // Inheriting basemodel parameters
    ParameterList parlist = submodel_->getParameters();

    for (int i = 0; i < parlist.size(); i++) {
        addParameter_(new Parameter("PIP." + parlist[i].getName(), parlist[i].getValue(), parlist[i].getConstraint()));
    }

    name_ = submodel_->getName() + "+PIP";
    modelname_ = "PIP." + submodel_->getName();


    addParameter_(new Parameter("PIP.lambda", lambda, &Parameter::R_PLUS));
    addParameter_(new Parameter("PIP.mu", mu, &Parameter::R_PLUS));

    // Compute frequencies from data
    if (initFreqFromData) {

        setFreqFromData(data, 0);

    } else {
        // Add frequency for gap character
        freq_ = submodel_->getFrequencies();
        freq_[alphabet_->getGapCharacterCode()] = 0; // hack for updateMatrices()

    }

    // Check if the sum of the frequencies is 1
    LOG_IF(WARNING, !ComparisonUtils::areLogicallyEqual(bpp::VectorTools::sum(freq_), 1.0)) <<
                                                                                            "The state frequencies sum up to " <<
                                                                                            TextTools::toString(bpp::VectorTools::sum(freq_), 15) <<
                                                                                            " != [1.0]";

    updateMatrices();
}

void PIP_Nuc::updateMatrices() {

    lambda_ = getParameterValue("lambda");
    mu_ = getParameterValue("mu");

    unsigned long eraseCharNum = name_.size();

    for (int i = 0; i < getParameters().size(); i++) {
        //test[i].getName();
        std::string parName = getParameters()[i].getName();
        if (parName.find(modelname_) != std::string::npos) {
            parName.erase(parName.begin(), parName.begin() + eraseCharNum + 1);
            submodel_->setParameterValue(parName, getParameters()[i].getValue());
        }
    }

    // Copy the generator from substitution model + extend it
    const bpp::Matrix<double> &qmatrix = submodel_->getGenerator();

    int cols = qmatrix.getNumberOfColumns();
    int rows = qmatrix.getNumberOfRows();


    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {

            if (i == j) {
                generator_(i, j) = qmatrix(i, j) - mu_;
            } else {

                generator_(i, j) = qmatrix(i, j);
            }

        }
        generator_(i, cols - 1) = mu_;
    }



    // Copy the exchangeability from substitution model + extend it
    const bpp::Matrix<double> &exMatrix = submodel_->getExchangeabilityMatrix();

    // Exchangeability
    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {

            if (i == j) {
                exchangeability_(i, j) = exMatrix(i, j) - mu_;
            } else {
                exchangeability_(i, j) = exMatrix(i, j);
            }

        }

        exchangeability_(i, cols - 1) = mu_;
    }

    // Normalization:
    setDiagonal();
    //normalize();


    // Compute eigen values and vectors:
    if (enableEigenDecomposition()) {
        EigenValue<double> ev(generator_);
        rightEigenVectors_ = ev.getV();
        eigenValues_ = ev.getRealEigenValues();
        iEigenValues_ = ev.getImagEigenValues();
        try {
            MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
            isNonSingular_ = true;
            isDiagonalizable_ = true;
            for (size_t i = 0; i < size_ && isDiagonalizable_; i++) {
                if (abs(iEigenValues_[i]) > NumConstants::TINY())
                    isDiagonalizable_ = false;
            }
        }
        catch (ZeroDivisionException &e) {
            ApplicationTools::displayMessage("Singularity during diagonalization. Taylor series used instead.");

            isNonSingular_ = false;
            isDiagonalizable_ = false;
            MatrixTools::Taylor(generator_, 30, vPowGen_);
        }
    }

    //freq_[alphabet_->getGapCharacterCode()] = 0;

}


void PIP_Nuc::setFreq(std::map<int, double> &freqs) {
    std::vector<double> values;
    for (auto const &element : freqs) {
        values.push_back(element.second);
    }
    freq_ = values;

}

void PIP_Nuc::setFreqFromData(const SequenceContainer &data, double pseudoCount) {
    std::map<int, int> counts;
    SequenceContainerTools::getCounts(data, counts);
    std::map<int, double> freqs;

    int gapkey = data.getAlphabet()->getGapCharacterCode();
    std::map<int, int>::iterator iter = counts.find(gapkey);
    if (iter != counts.end())
        counts.erase(iter);
    else puts("not found");


    std::vector<int> retval;
    for (auto const &element : counts) {
        retval.push_back(element.first);
    }

    double t = 0;
    for (auto &ci : counts)
        t += ci.second;

    t += pseudoCount * (double) counts.size();

    //for (int i = 0; i < static_cast<int>(size_)-1; i++)
    for (auto &key:retval) {
        freqs[key] = (static_cast<double>(counts[key]) + pseudoCount) / t;
    }

    freqs[data.getAlphabet()->getGapCharacterCode()] = 0;
    // Re-compute generator and eigen values:
    setFreq(freqs);

    // Check if the sum of the frequencies is 1
    LOG_IF(WARNING, !ComparisonUtils::areLogicallyEqual(bpp::VectorTools::sum(freq_), 1.0)) <<
                                                                                            "The state frequencies sum up to " <<
                                                                                            TextTools::toString(bpp::VectorTools::sum(freq_), 15) <<
                                                                                            " != [1.0]";
}


PIP_AA::PIP_AA(const ProteicAlphabet *alpha, SubstitutionModel *basemodel, const SequenceContainer &data, double lambda, double mu,
               bool initFreqFromData) :
        AbstractParameterAliasable("PIP."),
        AbstractReversibleProteinSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "PIP."),
        lambda_(lambda), mu_(mu), freqSet_(0) {

    // Setting basemodel to PIP and inherit freq set
    submodel_ = basemodel;

    // Inherit frequency set from basemodel and the associated parameters
    freqSet_ = const_cast<FrequenciesSet *>(submodel_->getFrequenciesSet());

    // Extending namespace
    std::string namespaceModel = submodel_->getName() + "+PIP." + freqSet_->getName();


    std::vector<double> subModelFreqs;
    // Compute frequencies from data
    if (initFreqFromData) {
        setFreqFromData(data, 0);
    }
    //else{
    // Add fixed frequency for gap character
    //     freq_ = freqSet_->getFrequencies();
    //}


    freqSet_->setNamespace(namespaceModel);

    // Set model name and prefix
    name_ = namespaceModel;
    modelname_ = "PIP." + submodel_->getName() + "." + freqSet_->getName();

    // Inheriting basemodel parameters (excluded frequency parameters)
    ParameterList parlist = basemodel->getParameters();

    for (int i = 0; i < parlist.size(); i++) {
        addParameter_(new Parameter("PIP." + parlist[i].getName(), parlist[i].getValue(), parlist[i].getConstraint()));
    }

    // Add PIP parameters
    addParameter_(new Parameter("PIP.lambda", lambda, &Parameter::R_PLUS_STAR));
    addParameter_(new Parameter("PIP.mu", mu, &Parameter::R_PLUS_STAR));



    //Update parameters and re-compute generator and eigen values:
    updateMatrices();
}


void PIP_AA::updateMatrices() {

    lambda_ = getParameterValue("lambda");
    mu_ = getParameterValue("mu");

    // Reset frequency for gap character
    freq_ = freqSet_->getFrequencies();
    freq_.resize(getNumberOfStates());

    // Check if the sum of the frequencies is 1
    LOG_IF(WARNING, !ComparisonUtils::areLogicallyEqual(bpp::VectorTools::sum(freq_), 1.0)) <<
                                                                                            "The state frequencies sum up to " <<
                                                                                            TextTools::toString(bpp::VectorTools::sum(freq_), 15) <<
                                                                                            " != [1.0]";

    unsigned long eraseCharNum = 4; // substring="PIP."

    for (int i = 0; i < getParameters().size(); i++) {
        //test[i].getName();
        std::string parName = getParameters()[i].getName();
        if (parName.find(modelname_) != std::string::npos) {
            parName.erase(parName.begin(), parName.begin() + eraseCharNum + submodel_->getName().size() + 1);
            submodel_->setParameterValue(parName, getParameters()[i].getValue());
        }
    }



    // Copy the generator from substitution model + extend it
    const bpp::Matrix<double> &qmatrix = submodel_->getGenerator();

    int cols = qmatrix.getNumberOfColumns() + 1;
    int rows = qmatrix.getNumberOfRows() + 1;


    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == j) {
                generator_(i, j) = qmatrix(i, j) - mu_;
            } else {
                generator_(i, j) = qmatrix(i, j);
            }
        }
        generator_(i, cols - 1) = mu_;
    }

    // Copy the exchangeability from substitution model + extend it
    const bpp::Matrix<double> &exMatrix = submodel_->getExchangeabilityMatrix();

    // Exchangeability
    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == j) {
                exchangeability_(i, j) = exMatrix(i, j) - mu_;
            } else {
                exchangeability_(i, j) = exMatrix(i, j);
            }
        }
        exchangeability_(i, cols - 1) = mu_;
    }

    // Normalization:
    setDiagonal();
    //normalize();


    // Compute eigen values and vectors:
    if (enableEigenDecomposition()) {
        EigenValue<double> ev(generator_);
        rightEigenVectors_ = ev.getV();
        eigenValues_ = ev.getRealEigenValues();
        iEigenValues_ = ev.getImagEigenValues();
        try {
            MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
            isNonSingular_ = true;
            isDiagonalizable_ = true;
            for (size_t i = 0; i < size_ && isDiagonalizable_; i++) {
                if (abs(iEigenValues_[i]) > NumConstants::TINY())
                    isDiagonalizable_ = false;
            }
        }
        catch (ZeroDivisionException &e) {
            ApplicationTools::displayMessage("Singularity during diagonalization. Taylor series used instead.");

            isNonSingular_ = false;
            isDiagonalizable_ = false;
            MatrixTools::Taylor(generator_, 30, vPowGen_);
        }
    }


}

void PIP_AA::setFreqFromData(const SequenceContainer &data, double pseudoCount) {
    std::map<int, int> counts;
    SequenceContainerTools::getCounts(data, counts);
    std::vector<double> frequencies(getNumberOfStates() - 1);

    int states = static_cast<int>(size_) - 1;
    double t = 0;
    for (int i = 0; i < states; i++) {
        t += (counts[i] + pseudoCount);
    }
    for (size_t i = 0; i < size_ - 1; ++i) frequencies[i] = (static_cast<double>(counts[static_cast<int>(i)]) + pseudoCount) / t;

    freqSet_->setFrequencies(frequencies);
    matchParametersValues(freqSet_->getParameters());

    freq_ = freqSet_->getFrequencies();
    // Reset frequency for gap character

    freq_.resize(getNumberOfStates());
    freq_[data.getAlphabet()->getGapCharacterCode()] = 0;

    LOG_IF(WARNING, !ComparisonUtils::areLogicallyEqual(bpp::VectorTools::sum(freq_), 1.0)) <<
                                                                                            "The state frequencies sum up to " <<
                                                                                            TextTools::toString(bpp::VectorTools::sum(freq_), 15) <<
                                                                                            " != [1.0]";
}


PIP_Codon::PIP_Codon(const CodonAlphabet_Extended *alpha, const GeneticCode *gc, SubstitutionModel *basemodel, const SequenceContainer &data,
                     double lambda, double mu, bool initFreqFromData) :
        AbstractParameterAliasable("PIP."),
        AbstractReversibleSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "PIP."),
        lambda_(lambda), mu_(mu), freqSet_(0) {

    computeFrequencies(false);

    // Setting basemodel to PIP and inherit freq set
    submodel_ = basemodel;

    // Inherit frequency set from basemodel and the associated parameters
    freqSet_ = const_cast<FrequenciesSet *>(submodel_->getFrequenciesSet());

    // Extending namespace
    std::string namespaceModel = submodel_->getName() + "+PIP." + freqSet_->getName();

    std::vector<double> subModelFreqs;
    // Compute frequencies from data
    if (initFreqFromData) {
        setFreqFromData(data, 0);
    } else {
        // Add fixed frequency for gap character
        freq_ = freqSet_->getFrequencies();
    }


    freqSet_->setNamespace(namespaceModel);

    // Set model name and prefix
    name_ = namespaceModel;
    modelname_ = "PIP." + submodel_->getName() + "." + freqSet_->getName();

    // Inheriting basemodel parameters (excluded frequency parameters)
    ParameterList parlist = basemodel->getParameters();

    for (int i = 0; i < parlist.size(); i++) {
        addParameter_(new Parameter("PIP." + parlist[i].getName(), parlist[i].getValue(), parlist[i].getConstraint()));
    }

    // Add PIP parameters
    addParameter_(new Parameter("PIP.lambda", lambda, &Parameter::R_PLUS_STAR));
    addParameter_(new Parameter("PIP.mu", mu, &Parameter::R_PLUS_STAR));


    //Update parameters and re-compute generator and eigen values:
    updateMatrices();

}

PIP_Codon::~PIP_Codon() {}

void PIP_Codon::updateMatrices() {

    lambda_ = getParameterValue("lambda");
    mu_ = getParameterValue("mu");

    // Reset frequency for gap character
    freq_ = freqSet_->getFrequencies();
    freq_.resize(getNumberOfStates());

    // Check if the sum of the frequencies is 1
    LOG_IF(WARNING, !ComparisonUtils::areLogicallyEqual(bpp::VectorTools::sum(freq_), 1.0)) <<
                                                                                            "The state frequencies sum up to " <<
                                                                                            TextTools::toString(bpp::VectorTools::sum(freq_), 15) <<
                                                                                            " != [1.0]";
    unsigned long eraseCharNum = 4; // substring="PIP."

    for (int i = 0; i < getParameters().size(); i++) {
        //test[i].getName();
        std::string parName = getParameters()[i].getName();
        if (parName.find(modelname_) != std::string::npos) {
            parName.erase(parName.begin(), parName.begin() + eraseCharNum + submodel_->getName().size() + 1);
            submodel_->setParameterValue(parName, getParameters()[i].getValue());
        }
    }

    // Copy the generator from substitution model + extend it
    const bpp::Matrix<double> &qmatrix = submodel_->getGenerator();

    int cols = qmatrix.getNumberOfColumns() + 1;
    int rows = qmatrix.getNumberOfRows() + 1;

    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == j) {
                generator_(i, j) = qmatrix(i, j) - mu_;
            } else {
                generator_(i, j) = qmatrix(i, j);
            }
        }
        generator_(i, cols - 1) = mu_;
    }

    // Copy the exchangeability from substitution model + extend it
    const bpp::Matrix<double> &exMatrix = submodel_->getExchangeabilityMatrix();

    // Exchangeability
    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == j) {
                exchangeability_(i, j) = exMatrix(i, j) - mu_;
            } else {
                exchangeability_(i, j) = exMatrix(i, j);
            }
        }
        exchangeability_(i, cols - 1) = mu_;
    }

    // Normalization:
    setDiagonal();
    //normalize();


    // Compute eigen values and vectors:
    if (enableEigenDecomposition()) {
        EigenValue<double> ev(generator_);
        rightEigenVectors_ = ev.getV();
        eigenValues_ = ev.getRealEigenValues();
        iEigenValues_ = ev.getImagEigenValues();
        try {
            MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
            isNonSingular_ = true;
            isDiagonalizable_ = true;
            for (size_t i = 0; i < size_ && isDiagonalizable_; i++) {
                if (abs(iEigenValues_[i]) > NumConstants::TINY())
                    isDiagonalizable_ = false;
            }
        }
        catch (ZeroDivisionException &e) {
            ApplicationTools::displayMessage("Singularity during diagonalization. Taylor series used instead.");

            isNonSingular_ = false;
            isDiagonalizable_ = false;
            MatrixTools::Taylor(generator_, 30, vPowGen_);
        }
    }

}

void PIP_Codon::setFreqFromData(const SequenceContainer &data, double pseudoCount) {
    std::map<int, int> counts;
    SequenceContainerTools::getCounts(data, counts);
    std::vector<double> frequencies(getNumberOfStates() - 1);

    int states = static_cast<int>(size_) - 1;
    double t = 0;
    for (int i = 0; i < states; i++) {
        t += (counts[i] + pseudoCount);
    }
    for (size_t i = 0; i < size_ - 1; ++i) frequencies[i] = (static_cast<double>(counts[static_cast<int>(i)]) + pseudoCount) / t;

    freqSet_->setFrequencies(frequencies);
    matchParametersValues(freqSet_->getParameters());

    freq_ = freqSet_->getFrequencies();
    // Reset frequency for gap character

    freq_.resize(getNumberOfStates());
    freq_[data.getAlphabet()->getGapCharacterCode()] = 0;

    // Check if the sum of the frequencies is 1
    LOG_IF(WARNING, !ComparisonUtils::areLogicallyEqual(bpp::VectorTools::sum(freq_), 1.0)) <<
                                                                                            "The state frequencies sum up to " <<
                                                                                            TextTools::toString(bpp::VectorTools::sum(freq_), 15) <<
                                                                                            " != [1.0]";
}


PIP_Codon::PIP_Codon(const PIP_Codon &pip_codon) :
        AbstractParameterAliasable("PIP."),
        AbstractReversibleSubstitutionModel(pip_codon) {}

PIP_Codon &PIP_Codon::operator=(const PIP_Codon &pip_codon) {
    return *this;
}


double bpp::estimateLambdaFromData(Tree *tree, SequenceContainer *sequences, double proportion) {
    double N = 0;
    double lambda = 0;

    // Get average sequence length
    std::vector<std::string> seqNames = sequences->getSequencesNames();
    for (auto &seqName : seqNames) {
        N += sequences->getSequence(seqName).size();
    }
    N = N / sequences->getNumberOfSequences();

    // Lambda estimates depends on the proportion of columns containing gaps
    lambda = (N * proportion) / tree->getTotalLength();

    return lambda;
}

double bpp::estimateMuFromData(Tree *tree, double proportion) {
    double mu = 0;

    mu = proportion / tree->getTotalLength();

    return mu;
}

/*
double bpp::estimateLambdaFromData(Tree *tree, SiteContainer *alignment) {
    double N = 0;
    double M = 0;

    // Alignment length with gaps
    M = alignment->getNumberOfSites();

    // Compute average sequence length without gaps
    for (auto &seqName : alignment->getSequencesNames()) {
        for (int i = 0; i < alignment->getNumberOfSites(); i++) {
            if (alignment->getSequence(seqName)[i] != alignment->getAlphabet()->getGapCharacterCode()) {
                N++;
            }
        }
    }
    N = N / alignment->getNumberOfSequences();

    return (M - N) / tree->getTotalLength();
}

double bpp::estimateMuFromData(Tree *tree, SiteContainer *alignment) {
    double N = 0;
    double M = 0;

    // Alignment length with gaps
    M = alignment->getNumberOfSites();

    // Compute average sequence length without gaps
    for (auto &seqName : alignment->getSequencesNames()) {
        for (int i = 0; i < alignment->getNumberOfSites(); i++) {
            if (alignment->getSequence(seqName)[i] != alignment->getAlphabet()->getGapCharacterCode()) {
                N++;
            }
        }
    }
    N = N / alignment->getNumberOfSequences();

    return (M - N) / (tree->getTotalLength() * N);
}

*/

double bpp::computeNH(SiteContainer *alignment){
    double nH = 0;
    for (int i = 0; i < alignment->getNumberOfSites(); i++) {
        bool containsGap = false;
        for (auto &seqName : alignment->getSequencesNames()) {
            if (alignment->getSequence(seqName)[i] == alignment->getAlphabet()->getGapCharacterCode()) {
                containsGap = true;
                break;
            }
        }
        if (!containsGap) nH++;
    }

    if (nH <= 0) nH = 1;
    return nH;
}

double bpp::estimateMuFromData(Tree *tree, SiteContainer *alignment) {

    double nH = computeNH(alignment);
    double nG = alignment->getNumberOfSites() - nH;

    double tmp = exp(1) * (nG / nH + 1);
    double t1 = MathUtils::LambertW(0, tmp);
    return (t1 - 1) / tree->getTotalLength();

}

double bpp::estimateLambdaFromData(Tree *tree, SiteContainer *alignment) {
    double nH = computeNH(alignment);
    double nG = alignment->getNumberOfSites() - nH;

    double tmp = MathUtils::LambertW(0,exp(1) * (nG / nH + 1));

    return (nH*exp(tmp-1 ) * (tmp-1))/tree->getTotalLength();

}
