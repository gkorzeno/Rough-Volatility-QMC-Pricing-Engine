#include "TermStructureStrategy.hpp"
#include <numeric>

BacktestResult TermStructureStrategy::run(
    const std::vector<std::vector<double>>& paths,
    double r, double T)
{
    return runDetailed(paths, r, T).aggregate;
}

TermStructureStrategy::RunResult TermStructureStrategy::runDetailed(
    const std::vector<std::vector<double>>& paths,
    double r, double T)
{
    const int nPaths = static_cast<int>(paths.size());
    RunResult result;

    // ── 1. Compute skew metrics per expiry ────────────────────────────────
    for (double expiry : config.surface.getExpiries())
        result.skewByExpiry.push_back(
            SkewAnalyzer::computeSkew(config.surface, expiry, config.spot, config.rate));

    // ── 2. Generate signals ───────────────────────────────────────────────
    result.signals = SkewAnalyzer::generateSignals(
        config.surface, config.spot, config.rate,
        config.skewThreshold, config.termSpreadThreshold);

    if (result.signals.empty()) {
        result.aggregate = {0,0,0,0,0,0, std::vector<double>(nPaths, 0.0)};
        return result;
    }

    const int nSignals = std::min(config.maxSignals,
                                  static_cast<int>(result.signals.size()));

    // ── 3. Execute each signal as a delta-hedged pair ─────────────────────
    std::vector<double> totalPnL(nPaths, 0.0);

    for (int s = 0; s < nSignals; ++s) {
        const auto& sig = result.signals[s];

        // Short leg
        DeltaHedgedStraddleStrategy shortLeg(
            sig.shortStrike, sig.shortIV,
            config.hedgeInterval, config.transactionCost);
        BacktestResult shortRes = shortLeg.run(paths, r, T);

        // Long leg (negate short straddle PnL)
        DeltaHedgedStraddleStrategy longLeg(
            sig.longStrike, sig.longIV,
            config.hedgeInterval, config.transactionCost);
        BacktestResult longRes = longLeg.run(paths, r, T);

        for (int i = 0; i < nPaths; ++i)
            totalPnL[i] += shortRes.pnlPath[i] - longRes.pnlPath[i];
    }

    // Normalise
    if (nSignals > 0) {
        const double scale = 1.0 / nSignals;
        for (double& p : totalPnL) p *= scale;
    }

    // ── 4. Statistics ─────────────────────────────────────────────────────
    const double mean = std::accumulate(totalPnL.begin(), totalPnL.end(), 0.0) / nPaths;
    double sq = 0.0;
    for (double p : totalPnL) sq += (p - mean) * (p - mean);
    const double stdDev = std::sqrt(sq / std::max(1, nPaths - 1));
    const double sharpe = (stdDev > 1e-12) ? mean / stdDev : 0.0;
    const double maxDD  = -*std::min_element(totalPnL.begin(), totalPnL.end());

    std::vector<double> sorted = totalPnL;
    std::sort(sorted.begin(), sorted.end());
    const int vi = std::max(0, static_cast<int>(std::floor(0.05*nPaths))-1);
    const double var95 = -sorted[vi];
    double esSum = 0.0; const int esc = std::max(1, vi+1);
    for (int i=0;i<esc;++i) esSum+=sorted[i];

    result.aggregate = {mean, stdDev, sharpe, maxDD, var95, -esSum/esc, sorted};
    return result;
}