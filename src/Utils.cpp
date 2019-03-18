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
 * @file utils.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 21 12 2017
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
#include <cfloat>

#include <glog/logging.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Distance/PGMA.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Phyl/Distance/NeighborJoining.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <Bpp/Text/KeyvalTools.h>

#include "Utils.hpp"
#include "Optimizators.hpp"
#include "UnifiedTSHomogeneousTreeLikelihood_Generic.hpp"

using namespace bpp;


void UtreeBppUtils::_traverseTree_b2u(Utree *in_tree, VirtualNode *target, bpp::Tree *refTree, int nodeId, treemap &tm) {


    if (refTree->isLeaf(nodeId)) {
        target->vnode_leaf = true;
    } else {
        target->vnode_leaf = false;
    }

    for (auto &sonId:refTree->getSonsId(nodeId)) {
        auto ichild = new VirtualNode();

        ichild->vnode_id = sonId;
        ichild->vnode_name = refTree->getNodeName(sonId);
        ichild->vnode_branchlength = refTree->getDistanceToFather(sonId);
        // Filling the bidirectional map
        tm.insert(nodeassoc(sonId, ichild->getVnode_id()));

        if (!refTree->isLeaf(sonId)) {

            ichild->vnode_leaf = false;

            target->connectNode(ichild);
            in_tree->addMember(ichild);
            _traverseTree_b2u(in_tree, ichild, refTree, sonId, tm);


        } else {

            // Set all the other directions to null
            ichild->_setNodeLeft(nullptr);
            ichild->_setNodeRight(nullptr);

            // Set the LEAF flag to true
            ichild->vnode_leaf = true;
            in_tree->addMember(ichild);
            target->connectNode(ichild);
        }
    }
}


void UtreeBppUtils::convertTree_b2u(bpp::Tree *in_tree, Utree *out_tree, treemap &tm) {
    int rootId = in_tree->getRootId();

    for (auto &sonId:in_tree->getSonsId(rootId)) {

        auto ichild = new VirtualNode;

        ichild->vnode_id = sonId;
        ichild->vnode_name = in_tree->getNodeName(sonId);
        ichild->vnode_branchlength = in_tree->getDistanceToFather(sonId);

        // Filling the bidirectional map
        tm.insert(nodeassoc(sonId, ichild->getVnode_id()));

        _traverseTree_b2u(out_tree, ichild, in_tree, sonId, tm);

        // Add this node as starting point of the tree
        out_tree->addMember(ichild, true);

    }

    // Collapse multiforcating trees to star tree pointing to the same pseudoroot
    // Pick root node at random within the node-vector

    out_tree->startVNodes.at(0)->_setNodeUp(out_tree->startVNodes.at(1));
    out_tree->startVNodes.at(1)->_setNodeUp(out_tree->startVNodes.at(0));

}


bpp::TreeTemplate<bpp::Node> *UtreeBppUtils::convertTree_u2b(tshlib::Utree *in_tree) {

    auto *RootNode = new bpp::Node;

    RootNode->setName(in_tree->rootnode->getNodeName());
    RootNode->setDistanceToFather(in_tree->rootnode->vnode_branchlength);

    _traverseTree_u2b(RootNode, in_tree->rootnode->getNodeLeft());
    _traverseTree_u2b(RootNode, in_tree->rootnode->getNodeRight());


    RootNode->setId((int) in_tree->listVNodes.size());

    // Set root node on new bpp tree
    auto *tree = new bpp::TreeTemplate<bpp::Node>();

    tree->setRootNode(RootNode);

    return tree;
}


void UtreeBppUtils::_traverseTree_u2b(bpp::Node *target, tshlib::VirtualNode *source) {

    auto *child = new bpp::Node;

    child->setName(source->getNodeName());
    child->setId(source->vnode_id);
    child->setDistanceToFather(source->vnode_branchlength);

    if (!source->isTerminalNode()) {

        _traverseTree_u2b(child, source->getNodeLeft());

        _traverseTree_u2b(child, source->getNodeRight());

    }

    target->addSon(child);


}


void UtreeBppUtils::associateNode2Alignment(bpp::SiteContainer *sites, tshlib::Utree *in_tree) {

    for (auto &node:in_tree->listVNodes) {

        if (node->isTerminalNode()) {

            std::vector<std::string> seqnames = sites->getSequencesNames();

            for (int i = 0; i < seqnames.size(); i++) {

                if (seqnames.at(i).compare(node->vnode_name) == 0) {
                    node->vnode_seqid = i;
                    break;
                }

            }


        }

    }
}


void UtreeBppUtils::associateNode2Alignment(bpp::SequenceContainer *sequences, tshlib::Utree *in_tree) {

    for (auto &node:in_tree->listVNodes) {

        if (node->isTerminalNode()) {

            std::vector<std::string> seqnames = sequences->getSequencesNames();

            for (int i = 0; i < seqnames.size(); i++) {

                if (seqnames.at(i).compare(node->vnode_name) == 0) {
                    node->vnode_seqid = i;
                    break;
                }

            }


        }

    }
}


void UtreeBppUtils::renameInternalNodes(bpp::Tree *in_tree, std::string prefix) {

    // Rename internal nodes with standard Vxx * where xx is a progressive number
    for (auto &nodeId:in_tree->getNodesId()) {

        if (!in_tree->hasNodeName(nodeId)) {

            std::string stringId;
            std::string stringName;

            stringId = std::to_string(nodeId);
            stringName = prefix + stringId;

            in_tree->setNodeName(nodeId, stringName);

        }
    }

}


std::vector<bpp::Node *>
UtreeBppUtils::remapNodeLists(std::vector<int> &inputList, bpp::TreeTemplate<bpp::Node> *tree, UtreeBppUtils::treemap tm) {

    std::vector<bpp::Node *> newList;

    for (auto &vnode:inputList) {

        newList.push_back(tree->getNode(tm.right.at(vnode)));
    }

    return newList;
}


void UtreeBppUtils::updateTree_b2u(bpp::TreeTemplate<bpp::Node> inBTree, tshlib::Utree *inUTree, UtreeBppUtils::treemap &tm) {

    std::vector<tshlib::VirtualNode *> nodelist;

    nodelist = inUTree->listVNodes;

    std::map<int, bpp::Node *> tempMap;

    // reset inBtree
    for (auto &bnode:inBTree.getNodes()) {

        tempMap.insert(std::pair<int, bpp::Node *>(bnode->getId(), bnode));

        bnode->removeSons();
        bnode->removeFather();

    }


    for (auto &vnode:nodelist) {

        std::cerr << "vnode " << vnode->getNodeName();
        if (!vnode->isTerminalNode()) {

            // get corrisponding sons in inBTree
            bpp::Node *leftBNode = tempMap[tm.right.at(vnode->getNodeLeft()->getVnode_id())];
            bpp::Node *rightBNode = tempMap[tm.right.at(vnode->getNodeRight()->getVnode_id())];

            // get corrisponding parent in inBTree
            bpp::Node *pNode = tempMap[tm.right.at(vnode->getVnode_id())];

            leftBNode->setFather(pNode);
            rightBNode->setFather(pNode);

            //Add new sons
            pNode->setSon(0, leftBNode);
            pNode->setSon(1, rightBNode);

            std::cerr << "\t internal";

        } else {

            std::cerr << "\t leaf";



        }
        // in case the current vnode is also the pseudo-root
        if (vnode == vnode->getNodeUp()->getNodeUp()) {
            std::cerr << "\tvnode pseudoroot";


            bpp::Node *leftBNode = tempMap[tm.right.at(vnode->getVnode_id())];
            bpp::Node *rightBNode = tempMap[tm.right.at(vnode->getNodeUp()->getVnode_id())];

            // get corrisponding parent in inBTree
            inBTree.getRootNode()->removeSons();


            leftBNode->setFather(inBTree.getRootNode());
            rightBNode->setFather(inBTree.getRootNode());

            inBTree.getRootNode()->setSon(0, leftBNode);
            inBTree.getRootNode()->setSon(1, rightBNode);

        }


        std::cerr << "\t done\n";

    }

}


void UtreeBppUtils::updateTree_u2b(bpp::Tree *inBTree, tshlib::Utree *inUTree, UtreeBppUtils::treemap &tm) {

}




double MatrixBppUtils::dotProd(const std::vector<double> *x, const std::vector<double> *y) {

    double val;

    val = 0.0;
    for (unsigned long i = 0; i < x->size(); i++) {
        val += (x->at(i) * y->at(i));
    }

    return val;
}


double MatrixBppUtils::dotProd(const bpp::ColMatrix<double> &x, const bpp::ColMatrix<double> &y) {

    double val;

    val = 0.0;
    for (unsigned long i = 0; i < x.getNumberOfRows(); i++) {
        val += (x(i, 0) * y(i, 0));
    }

    return val;
}


std::vector<double> MatrixBppUtils::cwiseProd(std::vector<double> *x, std::vector<double> *y) {

    std::vector<double> val;

    val.resize(x->size());

    for (unsigned long i = 0; i < x->size(); i++) {
        val.at(i) = (x->at(i) * y->at(i));
    }

    return val;
}


double MatrixBppUtils::sumVector(std::vector<double> *x) {

    double val;

    val = 0.0;
    for (unsigned long i = 0; i < x->size(); i++) {
        val += x->at(i);
    }

    return val;
}


std::vector<double> MatrixBppUtils::matrixVectorProd(bpp::RowMatrix<double> &M, std::vector<double> &A) {

    std::vector<double> B;
    B.resize(A.size());

    for (int i = 0; i < A.size(); i++) {
        B[i] = 0;
        for (int j = 0; j < A.size(); j++) {
            B[i] += M(i, j) * A.at(j);

        }
    }

    return B;
}


bpp::DistanceMatrix *InputUtils::parseDistanceMatrix(std::string filepath) {

    std::ifstream inFile;
    inFile.open(filepath);

    int matsize;

    inFile >> matsize;
    std::string character;

    auto outmatrix = new bpp::DistanceMatrix(matsize);
    outmatrix->resize(matsize);
    std::string stringline;

    int x = 0;
    int rownum = 0;
    while (!inFile.eof()) {

        std::getline(inFile, stringline);
        std::stringstream ss(stringline);

        std::string token;
        int y = 0;

        int colnum = 0;
        while (std::getline(ss, token, ' ')) {

            if (!token.empty()) {

                if (y == 0) {
                    std::string seqname = token;
                    (*outmatrix).setName(rownum, seqname);

                } else {
                    double value = std::stod(token);
                    (*outmatrix)(rownum, colnum) = value;
                    (*outmatrix)(colnum, rownum) = value;
                    colnum++;
                }

            }
            //(*outmatrix)(x,y) = std::stod(token);
            y++;


        }
        if (!stringline.empty()) {
            x++;
            rownum++;
        }
    }

    return outmatrix;

}


void OutputUtils::writeOutput2LOG(bpp::AbstractHomogeneousTreeLikelihood *tl) {
    bpp::ParameterList parModel;
    std::ostringstream oss;


    parModel = tl->getSubstitutionModelParameters();
    if (parModel.size() > 0) {
        oss << "model=" << tl->getModel()->getName() << "(";
        for (auto &parameterName:parModel.getParameterNames()) {
            oss << parameterName << "=" << parModel.getParameter(parameterName).getValue() << ",";
        }
        oss << ")";
        DLOG(INFO) << oss.str();
    }
    oss.clear();
    oss.str("");


    parModel = tl->getRateDistributionParameters();
    if (parModel.size() > 0) {
        oss << "rates=" << tl->getRateDistribution()->getName() << "(";
        for (auto &parameterName:parModel.getParameterNames()) {
            oss << parameterName << "=" << parModel.getParameter(parameterName).getValue() << ",";
        }
        oss << ")";
        DLOG(INFO) << oss.str();
    }
    oss.clear();
    oss.str("");

    parModel = tl->getBranchLengthsParameters();
    if (parModel.size() > 0) {
        oss << "branches=" << tl->getBranchLengthsParameters().size() << "(";
        for (auto &parameterName:parModel.getParameterNames()) {
            oss << parameterName << "=" << parModel.getParameter(parameterName).getValue() << ",";;
        }
        oss << ")";
        DLOG(INFO) << oss.str();
    }
    oss.clear();
    oss.str("");
}


std::string OutputUtils::TreeTools::writeTree2String(bpp::Tree *tree) {
    bpp::Newick treeWriter;
    bpp::TreeTemplate<bpp::Node> ttree(*tree);
    std::ostringstream oss;
    treeWriter.write(ttree, oss);
    std::string out = oss.str();
    return out;
}


void OutputUtils::writeOutput2JSON(bpp::AbstractHomogeneousTreeLikelihood *tl,
                                   bpp::SiteContainer *sites,
                                   std::map<std::string, std::string> &params,
                                   const string &suffix,
                                   bool suffixIsOptional,
                                   bool verbose,
                                   int warn) {

    std::string estimate_filename = ApplicationTools::getAFilePath("output.estimates.file", params, false, false, "none", true);

    boost::property_tree::ptree pt; // initial ptree structure for json output

    try {

        bpp::ParameterList parModel;
        std::ostringstream oss;

        pt.put("Input.sites.length", sites->getNumberOfSites());
        pt.put("Input.sites.sequences", sites->getNumberOfSequences());
        pt.put("Input.alphabet.states", sites->getAlphabet()->getNumberOfStates());
        pt.put("Input.alphabet.type", sites->getAlphabet()->getAlphabetType());

        parModel = tl->getSubstitutionModelParameters();
        if (parModel.size() > 0) {

            for (auto &parameterName:parModel.getParameterNames()) {
                pt.put("Model." + parameterName, parModel.getParameter(parameterName).getValue());
            }

        }

        // Add PIP intensity to the list of parameters to output (PIP.Intensity is a pseudo-parameter)
        if (tl->getSubstitutionModel()->getName().find("PIP") != std::string::npos) {
            double pip_intensity = tl->getSubstitutionModelParameters().getParameter("PIP.lambda").getValue() *
                                   tl->getSubstitutionModelParameters().getParameter("PIP.mu").getValue();
            pt.put("Model.PIP.intensity", pip_intensity);
        }

        parModel = tl->getRateDistributionParameters();
        if (parModel.size() > 0) {
            for (auto &parameterName:parModel.getParameterNames()) {
                pt.put("ASVR." + parameterName, parModel.getParameter(parameterName).getValue());
            }

        }

        parModel = tl->getBranchLengthsParameters();
        if (parModel.size() > 0) {
            for (auto &parameterName:parModel.getParameterNames()) {
                pt.put("Tree." + parameterName, parModel.getParameter(parameterName).getValue());
            }
        }

        pt.put("Final.LogLikelihood", tl->getLogLikelihood());

        boost::property_tree::write_json(estimate_filename, pt);

    } catch (std::exception const &e) {
        std::cerr << e.what() << std::endl;
    }


}


void OutputUtils::writeTreeAnnotations2TSV(bpp::Tree *tree, std::string outputfile) {

    std::ofstream outfile;
    outfile.open(outputfile);
    //lkFile << score;

    // Write header
    outfile << "Taxa\t";
    for (auto &attributeName:tree->getBranchPropertyNames(tree->getInnerNodesId().at(1))) {

        outfile << attributeName << "\t";
    }
    outfile << "\n";

    // Write content
    for (auto &nodeID:tree->getNodesId()) {
        outfile << tree->getNodeName(nodeID) << "\t";
        for (auto &attributeName:tree->getBranchPropertyNames(tree->getInnerNodesId().at(1))) {
            if (tree->hasBranchProperty(nodeID, attributeName)) {
                outfile << dynamic_cast<const BppString *>(tree->getBranchProperty(nodeID, attributeName))->toSTL() << "\t";
            } else {
                outfile << "\t";
            }
        }
        outfile << "\n";
    }
    outfile.close();

}


void OutputUtils::writeNexusMetaTree(std::vector<bpp::Tree *> trees, std::map<std::string, std::string> &params, const std::string &suffix,
                                     bool suffixIsOptional, bool verbose, int warn) {

    std::string PAR_output_tree_filename = ApplicationTools::getAFilePath("output.tree.file", params, false, false, "", true, "", 1);
    std::string PAR_output_tree_attributes = ApplicationTools::getStringParameter("output.tree.attributes", params, "", "", true, true);

    bool internaNodeNames = false;
    std::vector<std::string> attributeNames;

    // Split all the attributes to parse for making the tree beautiful
    StringTokenizer st(PAR_output_tree_attributes, ",");
    while (st.hasMoreToken()) {
        try {

            string param = st.nextToken();

            if (param == "InternalNodes") {

                internaNodeNames = true;

            } else if (param != "AllAttributes") {

                attributeNames.push_back(param);

            }

        } catch (ParameterNotFoundException &pnfe) {
            ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");

        }
    }

    std::ofstream out;
    out.open(PAR_output_tree_filename);
    // Checking the existence of specified file, and possibility to open it in write mode
    if (!out) { throw IOException("writeNexusMetaTree::write: failed to write to stream"); }

    out << "#NEXUS" << endl;
    out << endl;
    out << "BEGIN TREES;" << endl;

    //First, we retrieve all leaf names from all trees:
    std::vector<std::string> names;
    for (size_t i = 0; i < trees.size(); i++) {
        names = VectorTools::vectorUnion(names, trees[i]->getLeavesNames());
    }
    //... and create a translation map:
    map<string, size_t> translation;
    size_t code = 1;
    for (size_t i = 0; i < names.size(); i++) {
        translation[names[i]] = code++;
    }

    //Second we translate all leaf names to their corresponding code:
    vector<Tree *> translatedTrees(trees.size());
    for (size_t i = 0; i < trees.size(); i++) {
        vector<int> leavesId = trees[i]->getLeavesId();
        Tree *tree = dynamic_cast<Tree *>(trees[i]->clone());
        for (size_t j = 0; j < leavesId.size(); j++) {
            tree->setNodeName(leavesId[j], TextTools::toString(translation[tree->getNodeName(leavesId[j])]));
        }
        translatedTrees[i] = tree;
    }

    //Third we print the translation command:
    out << "  TRANSLATE";
    size_t count = 0;
    for (map<string, size_t>::iterator it = translation.begin(); it != translation.end(); it++) {
        out << endl << "    " << it->second << "\t" << it->first;
        count++;
        if (count < translation.size())
            out << ",";
    }
    out << ";";

    //Finally we print all tree descriptions:
    for (size_t i = 0; i < trees.size(); i++) {
        out << endl << "  TREE tree" << (i + 1) << " = "
            << OutputUtils::TreeTools::treeToParenthesis(*translatedTrees[i], internaNodeNames, attributeNames);
    }
    out << "END;" << endl;

    //Clean trees:
    for (size_t i = 0; i < translatedTrees.size(); i++) {
        delete translatedTrees[i];
    }

}


void OutputUtils::exportOutput(bpp::AbstractHomogeneousTreeLikelihood *tl,
                               bpp::SiteContainer *sites,
                               std::map<std::string, std::string> &params,
                               const string &suffix,
                               bool suffixIsOptional,
                               bool verbose,
                               int warn) {


    std::string estimate_filename = ApplicationTools::getAFilePath("output.estimates.file", params, false, false, "none", true);
    std::string estimate_format = ApplicationTools::getStringParameter("output.estimates.format", params, "json", "", true, true);

    if (estimate_filename.find("none") == std::string::npos) {

        size_t lastindex = estimate_filename.find_last_of(".");
        std::string rawname = estimate_filename.substr(0, lastindex);

        if (estimate_format.find("json") != std::string::npos) {
            ApplicationTools::displayResult("Output estimates format", TextTools::toString("json"));

            rawname = rawname + ".json";
            params["output.estimates.file"] = rawname;

            // Export estimates in a computer-readable format (JSON)
            writeOutput2JSON(tl, sites, params);


        } else {
            ApplicationTools::displayResult("Output estimates format", TextTools::toString("human-readable text"));

            rawname = rawname + ".log";
            params["output.estimates.file"] = rawname;
            // Output estimates in a text file (TEXT)
            writeOutput2Text(tl, sites, params);
        }

        ApplicationTools::displayResult("Output estimates to file", estimate_filename);

    }

}


void OutputUtils::writeOutput2Text(bpp::AbstractHomogeneousTreeLikelihood *tl,
                                   bpp::SiteContainer *sites,
                                   std::map<std::string, std::string> &params,
                                   const string &suffix,
                                   bool suffixIsOptional,
                                   bool verbose,
                                   int warn) {


    std::string parametersFile = ApplicationTools::getAFilePath("output.estimates.file", params, false, false, "none", true);
    bool withAlias = ApplicationTools::getBooleanParameter("output.estimates.alias", params, true, "", true, 0);

    ApplicationTools::displayResult("Output estimates to file", parametersFile);


    if (parametersFile != "none") {
        StlOutputStream out(new ofstream(parametersFile.c_str(), ios::out));

        int numParametersModel = 0;

        numParametersModel += tl->getTree().getNumberOfNodes() - 1;

        out << "# Log likelihood = ";
        out.setPrecision(20) << (tl->getLogLikelihood());
        out.endLine();
        out << "# Number of sites = ";
        out.setPrecision(20) << sites->getNumberOfSites();
        out.endLine();
        out.endLine();
        out << "# Substitution model parameters:";
        out.endLine();

        tl->getModel()->matchParametersValues(tl->getParameters());
        numParametersModel += tl->getModel()->getNumberOfParameters();
        PhylogeneticsApplicationTools::printParameters(tl->getModel(), out, 1, withAlias);

        out.endLine();
        (out << "# Rate distribution parameters:").endLine();
        tl->getRateDistribution()->matchParametersValues(tl->getParameters());
        numParametersModel += tl->getRateDistribution()->getNumberOfParameters();
        PhylogeneticsApplicationTools::printParameters(tl->getRateDistribution(), out, withAlias);
        out.endLine();
        out << "# Total number of parameters: " << numParametersModel;
        out.endLine();
    }

}

std::string OutputUtils::TreeTools::nodeToParenthesis(const Tree &tree, int nodeId, bool internalNodesNames,
                                                      std::vector<std::string> attributeNames) throw(NodeNotFoundException) {

    if (!tree.hasNode(nodeId))
        throw NodeNotFoundException("OutputUtils::TreeTools::nodeToParenthesis", nodeId);
    ostringstream s;
    if (tree.isLeaf(nodeId)) {
        s << tree.getNodeName(nodeId);
    } else {
        s << "(";
        vector<int> sonsId = tree.getSonsId(nodeId);
        s << OutputUtils::TreeTools::nodeToParenthesis(tree, sonsId[0], internalNodesNames, attributeNames);
        for (size_t i = 1; i < sonsId.size(); i++) {
            s << "," << OutputUtils::TreeTools::nodeToParenthesis(tree, sonsId[i], internalNodesNames, attributeNames);
        }
        s << ")";

        if (internalNodesNames && !tree.getNodeName(nodeId).empty()) {
            s << tree.getNodeName(nodeId);
        }
        //if (bootstrap)
        //{
        //if (tree.hasBranchProperty(nodeId, BOOTSTRAP))
        //    s << (dynamic_cast<const Number<double>*>(tree.getBranchProperty(nodeId, BOOTSTRAP))->getValue());
        //}
        //else
        //{

    }

    // Add metacomments (if not specified elsewhere, then all the attributes associated on the branches)
    if (attributeNames.empty()) {
        attributeNames = tree.getBranchPropertyNames(nodeId);
    }

    if (attributeNames.size() > 0) {
        s << "[&";

        for (int i = 0; i < attributeNames.size(); i++) {
            std::string propertyName = attributeNames.at(i);
            if (tree.hasBranchProperty(nodeId, propertyName))
                s << propertyName << "=" << *(dynamic_cast<const BppString *>(tree.getBranchProperty(nodeId, propertyName)));
            //}
            if (i < attributeNames.size() - 1) {
                s << ",";
            }
        }
        s << "]";
    }

    if (tree.hasDistanceToFather(nodeId))
        s << ":" << tree.getDistanceToFather(nodeId);
    return s.str();
}


std::string OutputUtils::TreeTools::treeToParenthesis(const Tree &tree, bool internalNodesNames, std::vector<std::string> attributeNames) {
    ostringstream s;
    s << "(";
    int rootId = tree.getRootId();
    vector<int> sonsId = tree.getSonsId(rootId);
    if (tree.isLeaf(rootId)) {
        s << tree.getNodeName(rootId);
        for (size_t i = 0; i < sonsId.size(); i++) {
            s << "," << OutputUtils::TreeTools::nodeToParenthesis(tree, sonsId[i], internalNodesNames, attributeNames);
        }
    } else {
        s << OutputUtils::TreeTools::nodeToParenthesis(tree, sonsId[0], internalNodesNames, attributeNames);
        for (size_t i = 1; i < sonsId.size(); i++) {
            s << "," << OutputUtils::TreeTools::nodeToParenthesis(tree, sonsId[i], internalNodesNames, attributeNames);
        }
    }
    s << ")";
    //if (bootstrap)
    //{
    //if (tree.hasBranchProperty(rootId, BOOTSTRAP))
    //    s << (dynamic_cast<const Number<double>*>(tree.getBranchProperty(rootId, BOOTSTRAP))->getValue());
    //}
    //else
    //{

    if (attributeNames.empty()) {
        attributeNames = tree.getBranchPropertyNames(rootId);
    }

    if (attributeNames.size() > 0) {
        s << "[&";

        for (int i = 0; i < attributeNames.size(); i++) {
            std::string propertyName = attributeNames.at(i);
            if (tree.hasBranchProperty(rootId, propertyName))
                s << propertyName << "=" << *(dynamic_cast<const BppString *>(tree.getBranchProperty(rootId, propertyName)));
            //}
            if (i < attributeNames.size() - 1) {
                s << ",";
            }
        }
        s << "]";
    }

    s << ";" << endl;
    return s.str();
}


bpp::DistanceEstimation DistanceUtils::computeDistanceMethod(std::string seqfilename, bpp::Alphabet *alphabet, bpp::GeneticCode *gCode,
                                                             std::map<std::string, std::string> &params) {

    // Create a map containing the required parameters passed by the user
    map<std::string, std::string> parmap;
    parmap = params;
    // Overwrite the model parameter
    parmap["model"] = "JC69";
    // Read the data using the non-extended alphabet
    bpp::Fasta seqReader;
    bpp::SequenceContainer *sequences = seqReader.readAlignment(seqfilename, alphabet);
    bpp::SiteContainer *sites = new bpp::VectorSiteContainer(*sequences);
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    bpp::TransitionModel *model = bpp::PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode, sites, params);

    bpp::DiscreteDistribution *rDist = 0;
    if (model->getNumberOfStates() > model->getAlphabet()->getSize()) {
        //Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();
    } else {
        rDist = bpp::PhylogeneticsApplicationTools::getRateDistribution(params);
    }

    bpp::DistanceEstimation distEstimation(model, rDist, sites, 1, false);

    delete sites;

    return distEstimation;
}


std::string TextUtils::appendToFilePath(std::string inputFilePath, std::string string2append) {

    // Split full filepath first and get filename

    std::vector<std::string> fullpath;
    boost::split(fullpath, inputFilePath, boost::is_any_of("/"));
    std::string filename = fullpath.back();
    fullpath.pop_back(); // remove element
    std::string joinedFullPath = boost::algorithm::join(fullpath, "/");

    // Split filename on dot
    std::vector<std::string> segmentsFileName;
    boost::split(segmentsFileName, filename, boost::is_any_of("."));
    std::reverse(std::begin(segmentsFileName), std::end(segmentsFileName));
    std::string lastelement = segmentsFileName.at(0);
    segmentsFileName.at(0) = string2append;
    std::reverse(std::begin(segmentsFileName), std::end(segmentsFileName));
    segmentsFileName.push_back(lastelement);
    std::string joinedString = boost::algorithm::join(segmentsFileName, ".");

    // recompone string
    std::string fullInitialMSAFilePath = joinedFullPath + "/" + joinedString;

    return fullInitialMSAFilePath;
}


void AlignmentUtils::checkAlignmentConsistency(bpp::SiteContainer &sites) {

    int gapCode = sites.getAlphabet()->getGapCharacterCode();
    int unresolvedCode = sites.getAlphabet()->getUnknownCharacterCode();
    bool nonGapSeen = false;
    bool nonUnkownSeen = false;
    int currentChar;

    for (unsigned long i = 0; i < sites.getNumberOfSites(); i++) {

        nonGapSeen = false;
        nonUnkownSeen = false;

        for (unsigned long s = 0; s < sites.getNumberOfSequences(); s++) {

            currentChar = sites.getSite(i).getValue(s);

            if (currentChar != gapCode) nonGapSeen = true;
            if (currentChar != unresolvedCode) nonUnkownSeen = true;

        }

        LOG_IF(FATAL, !nonGapSeen || !nonUnkownSeen)
        << "Column #" << i + 1 << " of the alignment contains only gaps. Please remove it and try again!";

    }
}


bool ComparisonUtils::areLogicallyEqual(double a, double b) {
    //return std::abs(a - b) < std::numeric_limits<double>::epsilon();
    //return a == b || std::abs(a - b) < std::abs(std::min(a, b)) * std::numeric_limits<double>::epsilon();
    return std::abs(a - b) < 0.00001 * std::max(std::abs(a), std::abs(b));

}

double MathUtils::add_lns(double a_ln,double b_ln){
    //ln(a + b) = ln{exp[ln(a) - ln(b)] + 1} + ln(b)

    double R;
    const double exp_precision =  log(pow(2,(double)DBL_MANT_DIG-1)-1);

    //ApplicationTools::displayResult("Mantissa precision", TextTools::toString(exp_precision, 50));

    if (std::isinf(a_ln) && std::isinf(b_ln)) {
        R = -std::numeric_limits<double>::infinity();
    } else if (std::isinf(a_ln)) {
        R = b_ln;
    } else if (std::isinf(b_ln)) {
        R = a_ln;
    } else if ((abs(a_ln - b_ln) >= exp_precision)) {
        //TODO:check this
        //2^52-1 = 4503599627370495.	log of that is 36.043653389117155867651465390794
        R = max(a_ln, b_ln);
    } else {
        R = log(exp(a_ln - b_ln) + 1) + b_ln;
    }

    return R;

};


