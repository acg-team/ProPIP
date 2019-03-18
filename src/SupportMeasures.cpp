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
 * @file SupportMeasures.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 02 05 2018
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

#include <glog/logging.h>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Seq/Io/Fasta.h>
#include "UnifiedTSHomogeneousTreeLikelihood_PIP.hpp"
#include "SupportMeasures.hpp"
#include "Optimizators.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_Generic.hpp"

using namespace bpp;

Bootstrap::Bootstrap(AbstractHomogeneousTreeLikelihood *tl,
                   const SiteContainer &data,
                   DiscreteDistribution *rDist,
                   tshlib::Utree *utree,
                   UtreeBppUtils::treemap *tm,
                   std::map<std::string, std::string> &params,
                   const std::string& suffix,
                   bool suffixIsOptional,
                   int warn) {

    unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("support.bootstrap.repeats", params, 0, suffix, true, 1);

    auto bestTree = new TreeTemplate<Node>(tl->getTree());
    Tree *tree = bestTree;

    bool indelModels = false;
    if (dynamic_cast<UnifiedTSHomogeneousTreeLikelihood_PIP *>(tl)) {
        indelModels = true;
    }

    AbstractHomogeneousTreeLikelihood *tlRep = nullptr;

    if (nbBS > 0) {
        ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
        bool approx = ApplicationTools::getBooleanParameter("support.bootstrap.approximate", params, true, suffix, true, 2);
        ApplicationTools::displayBooleanResult("Use approximate bootstrap", approx);
        bool bootstrapVerbose = ApplicationTools::getBooleanParameter("support.bootstrap.verbose", params, false, suffix, true, 2);

        const Tree *initTree = &tl->getTree();
        TransitionModel *model = tl->getModel();
        auto *sites = const_cast<SiteContainer *>(tl->getData());

        std::string bsTreesPath = ApplicationTools::getAFilePath("support.bootstrap.output.trees.file", params, false, false, suffix, true);
        std::string bsAlignPath = ApplicationTools::getAFilePath("support.bootstrap.output.alignments.file", params, false, false, suffix, true);

        std::ofstream *out_trees = nullptr;
        std::ofstream *out_msas = nullptr;

        if (bsTreesPath != "none") {
            ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
            out_trees = new std::ofstream(bsTreesPath.c_str(), std::ios::out);
        }
        if (bsAlignPath != "none") {
            ApplicationTools::displayResult("Bootstrap alignments stored in file", bsAlignPath);
        }

        std::string bootstrapMethod = ((approx) ? "Approximate" : "");

        Newick newick;
        ParameterList paramsToIgnore = tl->getSubstitutionModelParameters();
        paramsToIgnore.addParameters(tl->getRateDistributionParameters());

        ApplicationTools::displayTask("Bootstrapping", true);
        std::vector<Tree *> bsTrees(nbBS);

        bpp::Fasta seqWriter;

        for (unsigned int i = 0; i < nbBS; i++) {

           DLOG(INFO) << "[Bootstrap] " <<  bootstrapMethod << " cycle #" << i;

            ApplicationTools::displayGauge(i, nbBS - 1, '=');
            VectorSiteContainer *sample = SiteContainerTools::bootstrapSites(*sites);
            sample->reindexSites();

            // Write sampled alignments
            if (nbBS > 0) seqWriter.writeAlignment(TextUtils::appendToFilePath(bsAlignPath,std::to_string(i+1)), *sample, true);

            if (!approx) {
                tl->getModel()->setFreqFromData(*sample);
            }

            if (indelModels) {
                tlRep = new UnifiedTSHomogeneousTreeLikelihood_PIP(*initTree, *sample, model, rDist, utree, tm, true, params, "", false, false, false);
            } else {
                tlRep = new UnifiedTSHomogeneousTreeLikelihood(*initTree, *sample, model, rDist, utree, tm, true, params, "", false, false, false);
            }

            tlRep->initialize();
            ParameterList parametersRep = tlRep->getParameters();

            if (approx) {
                parametersRep.deleteParameters(paramsToIgnore.getParameterNames());
            }

            tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(Optimizators::optimizeParameters(tlRep,parametersRep, params, "", true, false, 0));

            bsTrees[i] = new TreeTemplate<Node>(tlRep->getTree());
            if (out_trees && i == 0) newick.write(*bsTrees[i], bsTreesPath, true);
            if (out_trees && i > 0) newick.write(*bsTrees[i], bsTreesPath, false);
            delete tlRep;
            delete sample;
        }

        // Close output files
        if (out_trees) out_trees->close();
        if (out_trees) delete out_trees;


        ApplicationTools::displayTaskDone();


        ApplicationTools::displayTask("Compute bootstrap values");
        TreeTools::computeBootstrapValues(*tree, bsTrees);
        ApplicationTools::displayTaskDone();
        for (unsigned int i = 0; i < nbBS; i++) {
            delete bsTrees[i];
        }

        // Write resulting tree:
        PhylogeneticsApplicationTools::writeTree(*bestTree, params);
    }

}

