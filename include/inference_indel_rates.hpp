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
 * @file inference_indel_rates.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 06 02 2018
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
 * @see For more information visit:
 */
#ifndef CASTOR_INFERENCE_INDEL_RATES_HPP
#define CASTOR_INFERENCE_INDEL_RATES_HPP

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

namespace inference_indel_rates{

    struct data {
        double tot_num_col;
        double num_col_1;
        double num_col_2;
        double bl_1;
        double bl_2;
    };

    int func_f(const gsl_vector *x, void *params, gsl_vector *f);

    void solve_system(gsl_vector *x0, gsl_multifit_nlinear_fdf *fdf,
            gsl_multifit_nlinear_parameters *params,double *lambda_0,double *mu_0);

    void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);

    void infere_indel_rates_from_sequences(std::string PAR_input_sequences,
                                           std::string PAR_Alphabet,
                                           bpp::Tree *tree,
                                           double *lambda_from_pairs,
                                           double *mu_from_pairs);

}

#endif //CASTOR_INFERENCE_INDEL_RATES_HPP
