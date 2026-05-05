#pragma once
#include <vector>
#include "RoughVolDiagnostics.hpp"
#include "FractionalBrownianMotion.hpp"
#include "FractionalSDE.hpp"
#include "RoughBergomiCalibrator.hpp"
#include "RoughBergomiModel.hpp"

class RoughVolatilityResearch {
public:
    using RoughBergomiParameters = RoughBergomiModel::Parameters;
    using RoughSurfaceSpec = RoughBergomiModel::SurfaceSpec;
    using CalibrationSpec = RoughBergomiCalibrator::CalibrationSpec;
    using CalibrationResult = RoughBergomiCalibrator::Result;
    using SmilePoint = RoughVolDiagnostics::SmilePoint;
    using TermPoint = RoughVolDiagnostics::TermPoint;
    using SviSurfaceFit = RoughVolDiagnostics::SviSurfaceFit;
    using VariancePoint = RoughVolDiagnostics::VariancePoint;
    using CalibrationNoisePoint = RoughVolDiagnostics::CalibrationNoisePoint;

    static std::vector<double> sampleFractionalBrownianMotion(
        double T, int steps, double hurst, Random& rng)
    {
        return FractionalBrownianMotion::samplePath(T, steps, hurst, rng);
    }

    static FractionalSDE::Path simulateFractionalSDE(
        double x0,
        double T,
        const std::vector<double>& fBmPath,
        const std::function<double(double, double)>& drift,
        const std::function<double(double, double)>& diffusion)
    {
        return FractionalSDE::simulateEuler(x0, T, fBmPath, drift, diffusion);
    }

    static VolSurface generateRoughBergomiSurface(
        const RoughBergomiParameters& params,
        const RoughSurfaceSpec& spec,
        int seed = 42,
        RoughBergomiModel::SamplingMethod method = RoughBergomiModel::MC,
        bool useOwenScrambling = true)
    {
        return RoughBergomiModel::generateVolSurface(params, spec, seed, method, useOwenScrambling);
    }

    static CalibrationResult calibrateRoughBergomi(
        const VolSurface& target,
        const CalibrationSpec& spec)
    {
        return RoughBergomiCalibrator::calibrate(target, spec);
    }

    static std::vector<SmilePoint> smileSlice(
        const VolSurface& surface,
        double expiry,
        const std::vector<double>& strikes)
    {
        return RoughVolDiagnostics::smileSlice(surface, expiry, strikes);
    }

    static std::vector<TermPoint> termStructureSlice(
        const VolSurface& surface,
        double strike)
    {
        return RoughVolDiagnostics::termStructureSlice(surface, strike);
    }

    static double smileLeftSkewMetric(const std::vector<SmilePoint>& smile) {
        return RoughVolDiagnostics::smileLeftSkewMetric(smile);
    }

    static SviSurfaceFit fitSviSurface(
        const VolSurface& surface,
        double spot,
        double rate)
    {
        return RoughVolDiagnostics::fitSviSurface(surface, spot, rate);
    }

    static std::vector<VariancePoint> varianceAnalysis(
        const RoughBergomiParameters& params,
        const RoughSurfaceSpec& spec,
        const std::vector<int>& pathCounts,
        int referencePaths,
        int outerRuns)
    {
        return RoughVolDiagnostics::varianceAnalysis(params, spec, pathCounts, referencePaths, outerRuns);
    }

    static std::vector<CalibrationNoisePoint> calibrationNoiseStudy(
        const RoughBergomiParameters& trueParams,
        const VolSurface& targetSurface,
        const CalibrationSpec& spec,
        const std::vector<int>& pathCounts,
        int outerRuns)
    {
        return RoughVolDiagnostics::calibrationNoiseStudy(trueParams, targetSurface, spec, pathCounts, outerRuns);
    }
};
