#include "ShortStraddleStrategy.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

BacktestResult ShortStraddleStrategy::run(
    const std::vector<std::vector<double>>& paths,
    double r,
    double T) {

    const int nPaths = static_cast<int>(paths.size());
    if (nPaths == 0) return {};

    // ── Premium collected at inception ────────────────────────────────────
    // S0 = paths[i][0] for any i (same starting price across all paths)
    const double S0 = paths[0][0];

    // Short straddle: sell call + sell put at strike K
    const double callPremium = BlackScholes::callPrice(S0, K, r, sigma, T);
    const double putPremium  = BlackScholes::putPrice(S0, K, r, sigma, T);
    const double premium     = callPremium + putPremium - tc; // net after transaction cost

    // ── PnL per path ─────────────────────────────────────────────────────
    // PnL = premium collected - payoff of straddle at expiry - transaction cost
    // (we already subtracted tc from premium above)
    std::vector<double> pnl(nPaths);
    for (int i = 0; i < nPaths; ++i) {
        const double ST       = paths[i].back();
        const double callPayoff = std::max(ST - K, 0.0);
        const double putPayoff  = std::max(K - ST, 0.0);
        // Discount payoffs back to t=0 to compare fairly with premium
        const double disc     = std::exp(-r * T);
        pnl[i] = premium - disc * (callPayoff + putPayoff);
    }

    // ── Statistics ───────────────────────────────────────────────────────
    const double meanPnL = std::accumulate(pnl.begin(), pnl.end(), 0.0) / nPaths;

    double sq = 0.0;
    for (double p : pnl) sq += (p - meanPnL) * (p - meanPnL);
    const double stdPnL = std::sqrt(sq / std::max(1, nPaths - 1));

    const double sharpe = (stdPnL > 1e-12) ? meanPnL / stdPnL : 0.0;

    // Max drawdown on cumulative PnL path (sequential order of paths)
    std::vector<double> cumPnL(nPaths);
    cumPnL[0] = pnl[0];
    for (int i = 1; i < nPaths; ++i)
        cumPnL[i] = cumPnL[i - 1] + pnl[i];

    double peak       = cumPnL[0];
    double maxDD      = 0.0;
    for (int i = 1; i < nPaths; ++i) {
        if (cumPnL[i] > peak) peak = cumPnL[i];
        const double dd = peak - cumPnL[i];
        if (dd > maxDD) maxDD = dd;
    }

    // VaR and ES at 95% confidence (left tail of PnL distribution)
    std::vector<double> sorted = pnl;
    std::sort(sorted.begin(), sorted.end());

    const int varIdx = static_cast<int>(std::floor(0.05 * nPaths));
    const double var95 = -sorted[std::max(0, varIdx)]; // sign convention: VaR is positive loss

    double esSum = 0.0;
    const int esCount = std::max(1, varIdx);
    for (int i = 0; i < esCount; ++i)
        esSum += sorted[i];
    const double es95 = -esSum / esCount;

    return BacktestResult{
        meanPnL,
        stdPnL,
        sharpe,
        maxDD,
        var95,
        es95,
        cumPnL   // pnlPath = cumulative PnL across paths
    };
}