#pragma once
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include "../surface/SVI.hpp"
#include "../surface/VolSurface.hpp"
#include "RoughBergomiCalibrator.hpp"
#include "RoughBergomiModel.hpp"

class RoughVolDiagnostics {
public:
    struct SmilePoint {
        double strike;
        double impliedVol;
    };

    struct TermPoint {
        double expiry;
        double impliedVol;
    };

    struct SviSliceFit {
        double expiry;
        SVIParams params;
        double rmse;
        double maxAbsError;
    };

    struct SviSurfaceFit {
        std::vector<SviSliceFit> slices;
        double overallRmse;
    };

    struct VariancePoint {
        int paths;
        double rmse;
        double rmseStdDev;
    };

    struct CalibrationNoisePoint {
        std::string samplerLabel;
        int pathsPerEvaluation;
        double parameterRmseMean;
        double parameterRmseStdDev;
        double objectiveRmseMean;
        double objectiveRmseStdDev;
        double xi0StdDev;
        double etaStdDev;
        double rhoStdDev;
        double hurstStdDev;
    };

    static std::vector<SmilePoint> smileSlice(
        const VolSurface& surface,
        double expiry,
        const std::vector<double>& strikes)
    {
        std::vector<SmilePoint> out;
        out.reserve(strikes.size());
        for (double K : strikes) {
            out.push_back({K, surface.impliedVol(expiry, K)});
        }
        return out;
    }

    static std::vector<TermPoint> termStructureSlice(
        const VolSurface& surface,
        double strike)
    {
        std::vector<TermPoint> out;
        const std::vector<double>& expiries = surface.getExpiries();
        out.reserve(expiries.size());
        for (double T : expiries) {
            out.push_back({T, surface.impliedVol(T, strike)});
        }
        return out;
    }

    static double smileLeftSkewMetric(const std::vector<SmilePoint>& smile) {
        if (smile.size() < 2)
            throw std::runtime_error("RoughVolDiagnostics: need at least two smile points");
        return smile.front().impliedVol - smile.back().impliedVol;
    }

    static SviSurfaceFit fitSviSurface(
        const VolSurface& surface,
        double spot,
        double rate)
    {
        const std::vector<double>& strikes = surface.getStrikes();
        const std::vector<double>& expiries = surface.getExpiries();

        SviSurfaceFit fit;
        fit.overallRmse = 0.0;
        int totalCount = 0;
        double totalSq = 0.0;

        for (double T : expiries) {
            const double F = spot * std::exp(rate * T);
            std::vector<double> marketVols;
            marketVols.reserve(strikes.size());
            for (double K : strikes)
                marketVols.push_back(surface.impliedVol(T, K));

            SVIParams params = SVI::fit(strikes, marketVols, F, T);

            double sq = 0.0;
            double maxAbs = 0.0;
            for (std::size_t i = 0; i < strikes.size(); ++i) {
                const double sviVol = SVI::impliedVol(params, strikes[i], F, T);
                const double diff = sviVol - marketVols[i];
                sq += diff * diff;
                maxAbs = std::max(maxAbs, std::abs(diff));
            }

            const double rmse = std::sqrt(sq / strikes.size());
            fit.slices.push_back({T, params, rmse, maxAbs});
            totalSq += sq;
            totalCount += static_cast<int>(strikes.size());
        }

        fit.overallRmse = std::sqrt(totalSq / std::max(1, totalCount));
        return fit;
    }

    static std::vector<VariancePoint> varianceAnalysis(
        const RoughBergomiModel::Parameters& params,
        const RoughBergomiModel::SurfaceSpec& baseSpec,
        const std::vector<int>& pathCounts,
        int referencePaths,
        int outerRuns)
    {
        if (outerRuns <= 0)
            throw std::runtime_error("RoughVolDiagnostics: outerRuns must be positive");

        RoughBergomiModel::SurfaceSpec referenceSpec = baseSpec;
        referenceSpec.paths = referencePaths;
        const VolSurface referenceSurface =
            averagedBenchmarkSurface(params, referenceSpec, std::max(6, outerRuns));

        std::vector<VariancePoint> diagnostics;
        diagnostics.reserve(pathCounts.size());

        for (int nPaths : pathCounts) {
            std::vector<double> rmses;
            rmses.reserve(outerRuns);

            RoughBergomiModel::SurfaceSpec testSpec = baseSpec;
            testSpec.paths = nPaths;

            for (int run = 0; run < outerRuns; ++run) {
                VolSurface trial = RoughBergomiModel::generateVolSurface(
                    params,
                    testSpec,
                    100 + run * 17 + nPaths);
                rmses.push_back(surfaceRmse(trial, referenceSurface));
            }

            diagnostics.push_back({
                nPaths,
                mean(rmses),
                stdDev(rmses)
            });
        }

        return diagnostics;
    }

    static void writeSmileCsv(
        const std::string& path,
        const std::vector<SmilePoint>& smile)
    {
        std::ofstream out(path.c_str());
        out << "strike,implied_vol\n";
        for (std::size_t i = 0; i < smile.size(); ++i)
            out << smile[i].strike << "," << smile[i].impliedVol << "\n";
    }

    static void writeTermCsv(
        const std::string& path,
        const std::vector<TermPoint>& term)
    {
        std::ofstream out(path.c_str());
        out << "expiry,implied_vol\n";
        for (std::size_t i = 0; i < term.size(); ++i)
            out << term[i].expiry << "," << term[i].impliedVol << "\n";
    }

    static void writeSviCsv(
        const std::string& path,
        const VolSurface& surface,
        const SviSurfaceFit& fit,
        double spot,
        double rate)
    {
        std::ofstream out(path.c_str());
        out << "expiry,strike,rough_iv,svi_iv,error\n";
        const std::vector<double>& strikes = surface.getStrikes();

        for (std::size_t si = 0; si < fit.slices.size(); ++si) {
            const double T = fit.slices[si].expiry;
            const double F = spot * std::exp(rate * T);
            for (std::size_t ki = 0; ki < strikes.size(); ++ki) {
                const double roughIv = surface.impliedVol(T, strikes[ki]);
                const double sviIv = SVI::impliedVol(fit.slices[si].params, strikes[ki], F, T);
                out << T << "," << strikes[ki] << "," << roughIv << "," << sviIv
                    << "," << (sviIv - roughIv) << "\n";
            }
        }
    }

    static void writeVarianceCsv(
        const std::string& path,
        const std::vector<VariancePoint>& variance)
    {
        std::ofstream out(path.c_str());
        out << "paths,rmse,rmse_stddev\n";
        for (std::size_t i = 0; i < variance.size(); ++i) {
            out << variance[i].paths << "," << variance[i].rmse
                << "," << variance[i].rmseStdDev << "\n";
        }
    }

    static std::vector<CalibrationNoisePoint> calibrationNoiseStudy(
        const RoughBergomiModel::Parameters& trueParams,
        const VolSurface& targetSurface,
        const RoughBergomiCalibrator::CalibrationSpec& baseSpec,
        const std::vector<int>& pathCounts,
        int outerRuns)
    {
        if (outerRuns <= 0)
            throw std::runtime_error("RoughVolDiagnostics: outerRuns must be positive");

        struct SamplerChoice {
            RoughBergomiModel::SamplingMethod method;
            const char* label;
            bool useOwen;
        };
        const SamplerChoice samplers[] = {
            {RoughBergomiModel::MC, "MC", false},
            {RoughBergomiModel::QMC, "QMC", false},
            {RoughBergomiModel::RQMC, "RQMC", true},
        };

        std::vector<CalibrationNoisePoint> out;
        for (const auto& sampler : samplers) {
            for (int paths : pathCounts) {
                std::vector<double> paramErrors;
                std::vector<double> objectiveErrors;
                std::vector<double> xi0s;
                std::vector<double> etas;
                std::vector<double> rhos;
                std::vector<double> hursts;

                for (int run = 0; run < outerRuns; ++run) {
                    RoughBergomiCalibrator::CalibrationSpec spec = baseSpec;
                    spec.pathsPerEvaluation = paths;
                    spec.objectiveSamplingMethod = sampler.method;
                    spec.objectiveUseOwenScrambling = sampler.useOwen;
                    spec.randomizeObjectiveSeed = (sampler.method != RoughBergomiModel::QMC);
                    spec.objectiveBaseSeed = 1000ULL + static_cast<std::uint64_t>(paths * 37 + run * 101);

                    try {
                        const RoughBergomiCalibrator::Result fit =
                            RoughBergomiCalibrator::calibrate(targetSurface, spec);

                        const double paramRmse =
                            std::sqrt(0.25 * (
                                sq(fit.parameters.xi0 - trueParams.xi0) +
                                sq(fit.parameters.eta - trueParams.eta) +
                                sq(fit.parameters.rho - trueParams.rho) +
                                sq(fit.parameters.hurst - trueParams.hurst)));

                        paramErrors.push_back(paramRmse);
                        objectiveErrors.push_back(fit.rmse);
                        xi0s.push_back(fit.parameters.xi0);
                        etas.push_back(fit.parameters.eta);
                        rhos.push_back(fit.parameters.rho);
                        hursts.push_back(fit.parameters.hurst);
                    } catch (...) {
                        paramErrors.push_back(1.0e6);
                        objectiveErrors.push_back(1.0e6);
                        xi0s.push_back(baseSpec.initialGuess.xi0);
                        etas.push_back(baseSpec.initialGuess.eta);
                        rhos.push_back(baseSpec.initialGuess.rho);
                        hursts.push_back(baseSpec.initialGuess.hurst);
                    }
                }

                out.push_back({
                    sampler.label,
                    paths,
                    mean(paramErrors),
                    stdDev(paramErrors),
                    mean(objectiveErrors),
                    stdDev(objectiveErrors),
                    stdDev(xi0s),
                    stdDev(etas),
                    stdDev(rhos),
                    stdDev(hursts)
                });
            }
        }

        return out;
    }

    static void writeCalibrationNoiseCsv(
        const std::string& path,
        const std::vector<CalibrationNoisePoint>& points)
    {
        std::ofstream out(path.c_str());
        out << "sampler,paths,parameter_rmse_mean,parameter_rmse_stddev,objective_rmse_mean,objective_rmse_stddev,xi0_stddev,eta_stddev,rho_stddev,hurst_stddev\n";
        for (const auto& p : points) {
            out << p.samplerLabel << "," << p.pathsPerEvaluation << ","
                << p.parameterRmseMean << "," << p.parameterRmseStdDev << ","
                << p.objectiveRmseMean << "," << p.objectiveRmseStdDev << ","
                << p.xi0StdDev << "," << p.etaStdDev << ","
                << p.rhoStdDev << "," << p.hurstStdDev << "\n";
        }
    }

private:
    static double sq(double x) { return x * x; }

    static double surfaceRmse(const VolSurface& a, const VolSurface& b) {
        const std::vector<double>& strikes = a.getStrikes();
        const std::vector<double>& expiries = a.getExpiries();
        double sq = 0.0;
        int count = 0;
        for (std::size_t i = 0; i < expiries.size(); ++i) {
            for (std::size_t j = 0; j < strikes.size(); ++j) {
                const double diff = a.impliedVol(expiries[i], strikes[j])
                                  - b.impliedVol(expiries[i], strikes[j]);
                sq += diff * diff;
                ++count;
            }
        }
        return std::sqrt(sq / std::max(1, count));
    }

    static double mean(const std::vector<double>& values) {
        return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    }

    static double stdDev(const std::vector<double>& values) {
        if (values.size() < 2)
            return 0.0;
        const double mu = mean(values);
        double sq = 0.0;
        for (std::size_t i = 0; i < values.size(); ++i) {
            const double d = values[i] - mu;
            sq += d * d;
        }
        return std::sqrt(sq / (values.size() - 1));
    }

    static VolSurface averagedBenchmarkSurface(
        const RoughBergomiModel::Parameters& params,
        const RoughBergomiModel::SurfaceSpec& spec,
        int runs)
    {
        const std::vector<double>& strikes = spec.strikes;
        const std::vector<double>& expiries = spec.expiries;
        std::vector<std::vector<double>> avg(
            expiries.size(),
            std::vector<double>(strikes.size(), 0.0));

        for (int run = 0; run < runs; ++run) {
            VolSurface surface = RoughBergomiModel::generateVolSurface(
                params,
                spec,
                777 + run * 101);

            for (std::size_t i = 0; i < expiries.size(); ++i) {
                for (std::size_t j = 0; j < strikes.size(); ++j) {
                    avg[i][j] += surface.impliedVol(expiries[i], strikes[j]);
                }
            }
        }

        for (std::size_t i = 0; i < expiries.size(); ++i)
            for (std::size_t j = 0; j < strikes.size(); ++j)
                avg[i][j] /= runs;

        return VolSurface(strikes, expiries, avg);
    }
};
