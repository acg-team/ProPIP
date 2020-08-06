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
 * @file DistanceFactoryML.cpp
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

#include <cmath>
#include <string>
#include <iostream>

#include "DistanceFactory.hpp"
#include "DistanceFactoryML.hpp"

#define MAXITER 20
#define EPSILON 1e-5

using namespace DistanceFactoryPrographMSA;

DistanceFactoryPrographMSA::DistanceFactoryML::DistanceFactoryML(bpp::SubstitutionModel *smodel,int ALPHABET_DIM) {

    this->smodel_=smodel;
    this->ALPHABET_DIM_=ALPHABET_DIM;
    this->CountMatrix_.resize(ALPHABET_DIM,ALPHABET_DIM);

    if(ALPHABET_DIM==4){
        DIST_MAX=DIST_MAX_DNA;
        VAR_MAX=VAR_MAX_DNA;
        VAR_MIN=VAR_MIN_DNA;
    }else if(ALPHABET_DIM==20){
        DIST_MAX=DIST_MAX_AA;
        VAR_MAX=VAR_MAX_AA;
        VAR_MIN=VAR_MIN_AA;
    }else{
        DIST_MAX=DIST_MAX_Codon;
        VAR_MAX=VAR_MAX_Codon;
        VAR_MIN=VAR_MIN_Codon;
    }

}

DistanceFactoryPrographMSA::DistanceFactoryML::~DistanceFactoryML() {
}

distvar_t DistanceFactoryML::computeMLDist(const CountMatrix &counts,
                                           unsigned int gaps,
                                           double seqlen,
                                           double dist0,
                                           double var0,
                                           bool mldist_gap_flag,
                                           double indel_rate) {

    double dist_min = 0;
    double dist_max = INFINITY;

    double dist = dist0, var = var0;

    double delta = 1;
    unsigned int iteration = 0;

    while (std::abs(delta) > EPSILON) {
        if (iteration > MAXITER) {
#ifdef DEBUG
            std::cerr << "maximum number of iterations reached: " << dist << " / " << dist0 << std::endl;
#endif
            if (dist_max == INFINITY) {
                dist = DIST_MAX;
                var = VAR_MAX;
            } else {
                dist = dist0;
                var = var0;
            }

            break;
        }

        bpp::RowMatrix<double> pr;
        pr=this->smodel_->getPij_t(dist);
        bpp::RowMatrix<double> q;
        q=this->smodel_->getGenerator();


        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> p;
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> qe;
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> pp;
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> ppp;

        p.resize(pr.getNumberOfRows(),pr.getNumberOfColumns());
        qe.resize(pr.getNumberOfRows(),pr.getNumberOfColumns());
        pp.resize(q.getNumberOfRows(),q.getNumberOfColumns());
        ppp.resize(pr.getNumberOfRows(),pr.getNumberOfColumns());

        for(int i=0;i<pr.getNumberOfRows();i++){
            for(int j=0;j<pr.getNumberOfColumns();j++){
                p(i,j)=pr(i,j);
                qe(i,j)=q(i,j);
            }
        }

        pp=qe*p;
        ppp=qe*pp;

        double f, ff;

        if (mldist_gap_flag) {
            double grate = indel_rate * seqlen * dist;
            //double g0 = std::pow(grate,gaps) * std::exp(-grate); // XXX / gaps!
            double g = (-grate + gaps) / dist;
            double gg = -gaps / (dist * dist);

            //lf = (p.array().log() * counts.cast<double>().array()).sum() + gaps * std::log(g0);
            f = (counts.template cast<double>().array() * pp.array() / p.array()).sum() + g;
            ff = ((counts.template cast<double>().array() * (ppp.array() * p.array() - pp.array().square())) /
                  p.array().square()).sum() + gg;
        } else {
            //lf = (p.array().log() * counts.cast<double>().array()).sum();
            f = (counts.template cast<double>().array() * pp.array() / p.array()).sum();
            ff = ((counts.template cast<double>().array() * (ppp.array() * p.array() - pp.array().square())) /
                  p.array().square()).sum();
        }

        /* variance estimated by Fisher information */
        var = -1.0 / ff;

        /* compute next estimate of distance by either newton or bisection */
        if (f > 0) {
            dist_min = std::max(dist_min, dist);
        } else {
            dist_max = std::min(dist_max, dist);
        }

        double new_dist = dist - f / ff;
        if (!(new_dist < dist_max && new_dist > dist_min)) {
            double upper = (dist_max == INFINITY) ? dist * 3 : dist_max;
            double lower = dist_min;
            new_dist = (upper + lower) / 2.0;
        }
        delta = 1.0 - new_dist / dist;
        dist = new_dist;

        ++iteration;
    }

    return distvar_t(dist, var);
}

distvar_t DistanceFactoryML::computeDistance(const CountMatrix &counts,
                                             unsigned int gaps,
                                             double seqlen,
                                             bool mldist_flag,
                                             bool mldist_gap_flag,
                                             double cutoff_dist,
                                             double indel_rate) {

    double ident = counts.diagonal().sum();
    double total = counts.sum();

    double dist0 = 1.0 - ident / total;
    double dist;
    double var;

    if (mldist_flag || mldist_gap_flag) {
        if (total == 0 || dist0 > 0.85) {
            dist = DIST_MAX;
            dist0 = DIST_MAX;
            var = VAR_MAX;
        } else {
            dist = dist0 = -std::log(1.0 - dist0 - 0.2 * dist0 * dist0);
            var = dist / total;
        }

        if (total > 0 && ident != total) {
            distvar_t dv = this->computeMLDist(counts, gaps, seqlen, dist, var, mldist_gap_flag, indel_rate);
            dist = dv.dist;
            var = dv.var;
        }
    } else {
        if (total == 0) {
            dist = 1.0;
            dist0 = 1.0;
            var = VAR_MAX;
        } else {
            dist = dist0;
            var = dist0 / total;
        }
    }

    if (!(dist < DIST_MAX)) {
        dist = DIST_MAX;
        var = VAR_MAX;
    }

    if (dist > cutoff_dist) {
        dist = cutoff_dist;
    }

    if (var < VAR_MIN) {
        var = VAR_MIN;
    }

    if (!(var < VAR_MAX)) {
        var = VAR_MAX;
    }

    return distvar_t(dist, var);
}

#undef EPSILON
#undef MAXITER