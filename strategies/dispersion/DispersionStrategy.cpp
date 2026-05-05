#include "DispersionStrategy.hpp"

DispersionStrategy::RunResult DispersionStrategy::runDispersion(
    double r, double T, int nPaths, std::uint64_t seed)
{
    const int nAssets = static_cast<int>(config.spotVec.size());
    const double dt   = 1.0 / 252;

    // ── 1. Compute dispersion signal ──────────────────────────────────────
    // Average pairwise correlation from the CorrelatedGBM
    double modelCorr = 0.0;
    {
        const auto& corrMatrix = config.process.correlation();
        int pairs = 0;
        for (int i = 0; i < nAssets; ++i)
            for (int j = i + 1; j < nAssets; ++j) {
                modelCorr += corrMatrix[i][j];
                ++pairs;
            }
        if (pairs > 0) modelCorr /= pairs;
    }

    const auto metrics = DispersionAnalyzer::analyze(
        config.indexImpliedVol,
        config.constituentIVs,
        config.basketWeights,
        modelCorr,
        config.entryThreshold);

    if (!metrics.dispersionOpportunity) {
        return {BacktestResult{0,0,0,0,0,0, std::vector<double>(nPaths, 0.0)},
                metrics, false, 0.0, 0.0};
    }

    // ── 2. Simulate correlated constituent paths via MultiAssetQMCSimulator
    // Result: terminalPrices[path][asset] — we need full paths for hedging.
    // Use RQMC for efficiency.
    const auto simResult = MultiAssetQMCSimulator::simulate(
        config.process, config.spotVec, T, dt, nPaths,
        MultiAssetQMCSimulator::RQMC, config.dirFile,
        /*useBrownianBridge=*/true, /*useOwenScrambling=*/true,
        /*rqmcReplicates=*/8, seed);

    // Reconstruct full paths per asset from terminal prices is not enough —
    // we need step-by-step paths to do delta hedging.
    // Re-simulate using PathSimulator logic with CorrelatedGBM.
    // Since MultiAssetQMCSimulator only returns terminals, we simulate
    // full paths directly using the multi-dim Euler-Maruyama integrator.
    const int nSteps = static_cast<int>(T / dt);

    // Full multi-asset paths: fullPaths[path][step][asset]
    std::vector<std::vector<std::vector<double>>> fullPaths(
        nPaths, std::vector<std::vector<double>>(
            nSteps + 1, std::vector<double>(nAssets)));

    {
        Random rng(static_cast<unsigned int>(seed));
        for (int i = 0; i < nPaths; ++i) {
            fullPaths[i][0] = config.spotVec;
            std::vector<double> x = config.spotVec;
            double t = 0.0;
            for (int j = 0; j < nSteps; ++j) {
                x = MultiDimensionalEulerMaruyama::step(
                    config.process, x, t, dt, rng);
                fullPaths[i][j+1] = x;
                t += dt;
            }
        }
    }

    // ── 3. Build basket paths and constituent paths ───────────────────────
    // basketPaths[i][j] = sum_k w_k * fullPaths[i][j][k]
    std::vector<std::vector<double>> basketPaths(
        nPaths, std::vector<double>(nSteps + 1));
    std::vector<std::vector<std::vector<double>>> constituentPaths(
        nAssets, std::vector<std::vector<double>>(
            nPaths, std::vector<double>(nSteps + 1)));

    for (int i = 0; i < nPaths; ++i) {
        for (int j = 0; j <= nSteps; ++j) {
            double basket = 0.0;
            for (int k = 0; k < nAssets; ++k) {
                basket += config.basketWeights[k] * fullPaths[i][j][k];
                constituentPaths[k][i][j] = fullPaths[i][j][k];
            }
            basketPaths[i][j] = basket;
        }
    }

    // ── 4. Run legs ───────────────────────────────────────────────────────
    // SHORT the basket (index) straddle — positive when basket vol < indexIV
    DeltaHedgedStraddleStrategy indexLeg(
        config.indexStrike, config.indexImpliedVol,
        config.hedgeInterval, config.transactionCost);
    BacktestResult indexRes = indexLeg.run(basketPaths, r, T);

    // LONG constituent straddles — positive when constituent vol > constituentIV
    // Net edge: we sell index vol at indexIV and buy constituent vol at constituentIV.
    // When realised constituent vol > constituentIV and realised basket vol < indexIV,
    // both legs profit.
    std::vector<double> totalPnL(nPaths, 0.0);
    for (int i = 0; i < nPaths; ++i)
        totalPnL[i] = indexRes.pnlPath[i];

    for (int k = 0; k < nAssets; ++k) {
        DeltaHedgedStraddleStrategy constLeg(
            config.constituentStrikes[k],
            config.constituentIVs[k],
            config.hedgeInterval,
            config.transactionCost);
        BacktestResult cRes = constLeg.run(constituentPaths[k], r, T);

        // LONG = negate the short-straddle PnL, weighted
        const double w = config.basketWeights[k];
        for (int i = 0; i < nPaths; ++i)
            totalPnL[i] -= w * cRes.pnlPath[i];
    }

    // ── 5. Realised vol diagnostics ───────────────────────────────────────
    const double rvIndex = realisedVol(basketPaths, T);
    double rvConstSum = 0.0;
    for (int k = 0; k < nAssets; ++k)
        rvConstSum += realisedVol(constituentPaths[k], T);
    const double rvConst = rvConstSum / nAssets;

    // ── 6. Statistics ─────────────────────────────────────────────────────
    const double mean = std::accumulate(totalPnL.begin(), totalPnL.end(), 0.0) / nPaths;
    double sq = 0.0;
    for (double p : totalPnL) sq += (p - mean) * (p - mean);
    const double stdDev = std::sqrt(sq / std::max(1, nPaths - 1));
    const double sharpe = (stdDev > 1e-12) ? mean / stdDev : 0.0;
    const double maxDD  = -*std::min_element(totalPnL.begin(), totalPnL.end());

    std::vector<double> sorted = totalPnL;
    std::sort(sorted.begin(), sorted.end());
    const int vi = std::max(0, static_cast<int>(std::floor(0.05 * nPaths)) - 1);
    const double var95 = -sorted[vi];
    double esSum = 0.0;
    const int esc = std::max(1, vi + 1);
    for (int i = 0; i < esc; ++i) esSum += sorted[i];

    BacktestResult agg{mean, stdDev, sharpe, maxDD, var95, -esSum/esc, sorted};
    return {agg, metrics, true, rvIndex, rvConst};
}