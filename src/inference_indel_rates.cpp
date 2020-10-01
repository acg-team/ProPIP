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
 * @file inference_indel_rates.cpp
 * @author Massimo Maiolo
 * @date 06 08 2020
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
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

/*
* From GSL:
*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

/*
* From SeqLib:
*/
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include <Bpp/Seq/AlphabetIndex/DefaultNucleotideScore.h>
#include <Bpp/Seq/AlphabetIndex/GranthamAAChemicalDistance.h>

#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/Node.h>

#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

#include <Bpp/Phyl/Model/SubstitutionModel.h>

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

#include <Bpp/Numeric/Matrix/Matrix.h>

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include <Bpp/Phyl/Distance/BioNJ.h>

#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include <glog/logging.h>

#include "inference_indel_rates.hpp"

int inference_indel_rates::func_f (const gsl_vector *x, void *params, gsl_vector *f){


    double ntotcol = ((struct data *)params)->tot_num_col;
    double ncol1 = ((struct data *)params)->num_col_1;
    double ncol2 = ((struct data *)params)->num_col_2;
    double bl1 = ((struct data *)params)->bl_1;
    double bl2 = ((struct data *)params)->bl_2;

    double lambda = gsl_vector_get(x, 0);
    double mu = gsl_vector_get(x, 1);


    gsl_vector_set(f, 0, lambda/mu - ntotcol);
    gsl_vector_set(f, 1, lambda/mu * (exp(-mu*bl1)) * (1-exp(-mu*bl2)) + lambda/mu * (1- exp(-mu*bl1)) - ncol2 );
    gsl_vector_set(f, 2, lambda/mu * (exp(-mu*bl2)) * (1-exp(-mu*bl1)) + lambda/mu * (1- exp(-mu*bl2)) - ncol1);

    return GSL_SUCCESS;
}

void inference_indel_rates::callback(const size_t iter, void *params,const gsl_multifit_nlinear_workspace *w){

    gsl_vector * x = gsl_multifit_nlinear_position(w);

    /* print out current location */
    printf("%f %f\n",gsl_vector_get(x, 0),gsl_vector_get(x, 1));
}

void inference_indel_rates::solve_system(gsl_vector *x0, gsl_multifit_nlinear_fdf *fdf,
                                         gsl_multifit_nlinear_parameters *params,double *lambda_0,double *mu_0){

    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    const size_t max_iter = 200;
    const double xtol = 1.0e-8;
    const double gtol = 1.0e-8;
    const double ftol = 1.0e-8;

    const size_t n = fdf->n;
    const size_t p = fdf->p;

    gsl_multifit_nlinear_workspace *work = gsl_multifit_nlinear_alloc(T, params, n, p);
    gsl_vector * f = gsl_multifit_nlinear_residual(work);
    gsl_vector * x = gsl_multifit_nlinear_position(work);

    int info;
    double chisq0, chisq, rcond;

    /* initialize solver */
    gsl_multifit_nlinear_init(x0, fdf, work);

    /* store initial cost */
    gsl_blas_ddot(f, f, &chisq0);

    /* iterate until convergence */
    //gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,callback, NULL, &info, work);
    gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,NULL, NULL, &info, work);

    /* store final cost */
    gsl_blas_ddot(f, f, &chisq);

    /* store cond(J(x)) */
    gsl_multifit_nlinear_rcond(&rcond, work);

    /*
    // print summary
    fprintf(stderr, "NITER         = %zu\n", gsl_multifit_nlinear_niter(work));
    fprintf(stderr, "NFEV          = %zu\n", fdf->nevalf);
    fprintf(stderr, "NJEV          = %zu\n", fdf->nevaldf);
    fprintf(stderr, "NAEV          = %zu\n", fdf->nevalfvv);
    fprintf(stderr, "initial cost  = %.12e\n", chisq0);
    fprintf(stderr, "final cost    = %.12e\n", chisq);
    fprintf(stderr, "final x       = (%.12e, %.12e)\n",gsl_vector_get(x, 0), gsl_vector_get(x, 1));
    fprintf(stderr, "final cond(J) = %.12e\n", 1.0 / rcond);
    */

    *lambda_0 = gsl_vector_get(x, 0);
    *mu_0 = gsl_vector_get(x, 1);

    gsl_multifit_nlinear_free(work);
}

void inference_indel_rates::infere_indel_rates_from_sequences(std::string PAR_input_sequences,
                                                              std::string PAR_Alphabet,
                                                              bool PAR_alignment,
                                                              bool PAR_model_indels,
                                                              std::map<std::string, std::string> &params,
                                                              bpp::Tree *tree,
                                                              double *lambda_from_pairs,
                                                              double *mu_from_pairs,
                                                              const bpp::GeneticCode*  gCode,
                                                              std::map<std::string, std::string> modelMap){

    *lambda_from_pairs=0.0;
    *mu_from_pairs=0.0;

    int id = 0;
    int Id1 = 0;
    int Id2 = 0;
    int nseq = 0;
    int gapIndex = 0;
    int localRootId = 0;
    int s1 = 0;
    int s2 = 0;
    double countGapCol_1=0.0;
    double countGapCol_2=0.0;
    double countNoGapCol=0.0;
    double tot_num_col = 0.0;
    double count_pairs=0.0;
    double bl1 = 0.0;
    double bl2 = 0.0;
    double lambda_0 = 0.0;
    double mu_0 = 0.0;
    bpp::Fasta seqReader;
    bpp::SequenceContainer *sequencesCanonical = nullptr;
    std::vector<std::string> seqnames;
    bpp::VectorSiteContainer *allSites = nullptr;
    bpp::SiteContainer *sites = nullptr;

    const bpp::Alphabet* alpha;

    double gap_penalty = -10.0;
    const size_t n = 3; // num. of equations
    const size_t p = 2; // num. of unknown

    if (PAR_Alphabet.find("DNA") != std::string::npos ) {
        alpha = &bpp::AlphabetTools::DNA_ALPHABET;
    } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
        alpha = &bpp::AlphabetTools::PROTEIN_ALPHABET;
    } else{

    }

    if(PAR_alignment){
        sequencesCanonical = seqReader.readSequences(PAR_input_sequences, alpha);
        seqnames = sequencesCanonical->getSequencesNames();
        nseq=seqnames.size();
    }else {
        allSites = bpp::SequenceApplicationTools::getSiteContainer(alpha, params);
        sites = bpp::SequenceApplicationTools::getSitesToAnalyse(*allSites, params, "", true, !PAR_model_indels,true, 1);
        seqnames = sites->getSequencesNames();
        nseq=sites->getNumberOfSequences();
    }

    count_pairs=0.0;
    for(int i=0;i<nseq-1;i++){
        for(int j=i+1;j<nseq;j++){

            bpp::SiteContainer* alignedSeq;

            if(PAR_alignment){

                if (PAR_Alphabet.find("DNA") != std::string::npos) {
                    alignedSeq = bpp::SiteContainerTools::alignNW(sequencesCanonical->getSequence(seqnames.at(i)),
                                                                  sequencesCanonical->getSequence(seqnames.at(j)),
                                                                  bpp::DefaultNucleotideScore(
                                                                          &bpp::AlphabetTools::DNA_ALPHABET),
                                                                  gap_penalty);
                } else if (PAR_Alphabet.find("Protein") != std::string::npos) {
                    alignedSeq = bpp::SiteContainerTools::alignNW(sequencesCanonical->getSequence(seqnames.at(i)),
                                                                  sequencesCanonical->getSequence(seqnames.at(j)),
                                                                  bpp::GranthamAAChemicalDistance(),
                                                                  gap_penalty);
                } else {

                }

            }else {

                bpp::VectorSequenceContainer *msa = new bpp::VectorSequenceContainer(alpha);

                msa->addSequence(*(new bpp::BasicSequence(sites->getName(i), sites->getSequence(i).toString(), alpha)), true);
                msa->addSequence(*(new bpp::BasicSequence(sites->getName(j), sites->getSequence(j).toString(), alpha)), true);

                alignedSeq = new bpp::VectorSiteContainer(*msa);

            }

            gapIndex = alignedSeq->getAlphabet()->getGapCharacterCode();

            countGapCol_1=0.0;
            countGapCol_2=0.0;
            countNoGapCol=0.0;
            for(int k=0;k<alignedSeq->getNumberOfSites();k++){
                bpp::Site s = alignedSeq->getSite(k);

                s1 = s.getValue(0);
                s2 = s.getValue(1);
                if(s1==gapIndex & s1==s2){
                    //DLOG(WARNING) << "full gap column!";
                }else if(s1==gapIndex){
                    countGapCol_1+=1.0;
                }else if(s2==gapIndex){
                    countGapCol_2+=1.0;
                }else{
                    countNoGapCol+=1.0;
                }
            }

            tot_num_col=(countNoGapCol+countGapCol_1+countGapCol_2);

            Id1 = tree->getLeafId(seqnames.at(i));
            Id2 = tree->getLeafId(seqnames.at(j));

            std::vector< int > nodeIds;
            nodeIds.push_back(Id1);
            nodeIds.push_back(Id2);

            localRootId = bpp::TreeTools::getLastCommonAncestor(*tree,nodeIds);

            id=Id1;
            bl1=0;
            while(id!=localRootId){
                bl1+=tree->getDistanceToFather(id);
                id=tree->getFatherId(id);
            }

            id=Id2;
            bl2=0;
            while(id!=localRootId){
                bl2+=tree->getDistanceToFather(id);
                id=tree->getFatherId(id);
            }

            gsl_vector *f = gsl_vector_alloc(n);
            gsl_vector *x = gsl_vector_alloc(p);
            gsl_multifit_nlinear_fdf fdf;
            gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

            inference_indel_rates:: data d= {tot_num_col,
                                             countGapCol_1,
                                             countGapCol_2,
                                             bl1,
                                             bl2};

            /* define function to be minimized */
            fdf.f = inference_indel_rates::func_f;
            fdf.df = NULL;
            fdf.fvv = NULL;
            fdf.n = n;
            fdf.p = p;
            fdf.params = &d;

            /* starting point */
            gsl_vector_set(x, 0, tot_num_col/10.0); //lambda
            gsl_vector_set(x, 1, 1.0/10.0); //mu

            lambda_0 = 0.0;
            mu_0 = 0.0;

            // This selects the Levenberg-Marquardt algorithm.
            fdf_params.trs = gsl_multifit_nlinear_trs_lm;
            inference_indel_rates::solve_system(x, &fdf, &fdf_params,&lambda_0,&mu_0);

            gsl_vector_free(f);
            gsl_vector_free(x);

            *lambda_from_pairs+=lambda_0;
            *mu_from_pairs+=mu_0;
            count_pairs+=1.0;
        }
    }

    *lambda_from_pairs=(*lambda_from_pairs/count_pairs);
    *mu_from_pairs=(*mu_from_pairs/count_pairs);

}
