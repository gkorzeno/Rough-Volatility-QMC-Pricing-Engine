#include "VRPStrategy.hpp"

BacktestResult VRPStrategy::run(
    const std::vector<std::vector<double>>& paths,
    double r, double T)
{
    if (paths.empty()) return {};
    const double spot = paths[0][0];

    // ── 1. Compute signal ─────────────────────────────────────────────────
    const Signal sig = computeSignal(spot, r, T);

    // If the signal says don't trade, return flat result
    if (!sig.shouldTrade) {
        const int n = static_cast<int>(paths.size());
        return BacktestResult{0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              std::vector<double>(n, 0.0)};
    }

    // ── 2. Delegate to delta-hedged straddle ──────────────────────────────
    // Size position proportionally to the VRP signal strength
    // (larger VRP = larger edge = larger position, capped at 1x).
    const double positionScale = std::min(1.0, sig.vrp / 0.05);

    DeltaHedgedStraddleStrategy hedged(
        config.strike,
        config.impliedVol,
        config.hedgeInterval,
        config.transactionCost);

    BacktestResult result = hedged.run(paths, r, T);

    // Scale all PnL statistics by position size
    result.meanPnL    *= positionScale;
    result.stdPnL     *= positionScale;
    result.maxDrawdown *= positionScale;
    result.var95      *= positionScale;
    result.es95       *= positionScale;
    for (double& p : result.pnlPath) p *= positionScale;
    // Sharpe is scale-invariant

    return result;
}