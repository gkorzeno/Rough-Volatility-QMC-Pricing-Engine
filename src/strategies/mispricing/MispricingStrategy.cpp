#include "MispricingStrategy.hpp"
#include <numeric>

std::vector<double> MispricingStrategy::runLeg(
    const std::vector<std::vector<double>>& paths,
    double strike, double impliedVol,
    double r, double T, bool isShort) const
{
    DeltaHedgedStraddleStrategy leg(
        strike, impliedVol, config.hedgeInterval, config.transactionCost);

    BacktestResult res = leg.run(paths, r, T);

    // For a long leg, negate the PnL (buying instead of selling)
    if (!isShort)
        for (double& p : res.pnlPath) p = -p;

    return res.pnlPath;
}

MispricingStrategy::RunResult MispricingStrategy::runDetailed(
    const std::vector<std::vector<double>>& paths,
    double r, double T)
{
    // ── 1. Scan for mispricings ───────────────────────────────────────────
    auto mispricings = MispricingScanner::scan(
        config.market, config.sviSurface,
        config.rbergomiSurface, config.localvolSurface,
        config.minMispricing);

    auto trades = MispricingScanner::buildRVTrades(mispricings, config.maxLegs);

    if (trades.empty()) {
        const int n = static_cast<int>(paths.size());
        return {BacktestResult{0,0,0,0,0,0, std::vector<double>(n,0)},
                {}, mispricings};
    }

    // ── 2. Run each trade leg, sum PnL ────────────────────────────────────
    const int nPaths = static_cast<int>(paths.size());
    std::vector<double> totalPnL(nPaths, 0.0);

    for (const auto& trade : trades) {
        // Long leg (buy undervalued)
        auto longPnL = runLeg(paths, trade.longStrike, trade.longIV, r, T, false);
        // Short leg (sell overvalued)
        auto shortPnL = runLeg(paths, trade.shortStrike, trade.shortIV, r, T, true);

        for (int i = 0; i < nPaths; ++i)
            totalPnL[i] += longPnL[i] + shortPnL[i];
    }

    // Normalise by number of trades so scale is per-pair, not per-leg-count
    if (!trades.empty()) {
        const double scale = 1.0 / trades.size();
        for (double& p : totalPnL) p *= scale;
    }

    // ── 3. Compute aggregate statistics ───────────────────────────────────
    const double mean = std::accumulate(totalPnL.begin(), totalPnL.end(), 0.0) / nPaths;
    double sq = 0.0;
    for (double p : totalPnL) sq += (p - mean) * (p - mean);
    const double stdDev = std::sqrt(sq / std::max(1, nPaths - 1));
    const double sharpe = (stdDev > 1e-12) ? mean / stdDev : 0.0;
    const double maxDD  = -*std::min_element(totalPnL.begin(), totalPnL.end());

    std::vector<double> sorted = totalPnL;
    std::sort(sorted.begin(), sorted.end());
    const int    varIdx = std::max(0, static_cast<int>(std::floor(0.05 * nPaths)) - 1);
    const double var95  = -sorted[varIdx];
    double esSum = 0.0;
    const int esCount = std::max(1, varIdx + 1);
    for (int i = 0; i < esCount; ++i) esSum += sorted[i];
    const double es95 = -esSum / esCount;

    BacktestResult agg{mean, stdDev, sharpe, maxDD, var95, es95, sorted};
    return {agg, trades, mispricings};
}

BacktestResult MispricingStrategy::run(
    const std::vector<std::vector<double>>& paths,
    double r, double T)
{
    return runDetailed(paths, r, T).aggregate;
}