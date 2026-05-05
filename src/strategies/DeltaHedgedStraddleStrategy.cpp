#include "DeltaHedgedStraddleStrategy.hpp"
#include <numeric>
#include <algorithm>

BacktestResult DeltaHedgedStraddleStrategy::run(
    const std::vector<std::vector<double>>& paths,
    double r,
    double T)
{
    const int nPaths = static_cast<int>(paths.size());
    if (nPaths == 0) return {};
    if (paths[0].empty()) return {};

    const double S0     = paths[0][0];
    const int    nSteps = static_cast<int>(paths[0].size()) - 1;
    const double dtPath = T / nSteps;

    const int stepsPerHedge = std::max(1,
        static_cast<int>(std::round(dtHedge / dtPath)));

    // ── Inception premium ─────────────────────────────────────────────────
    // This is income. We subtract one round-trip transaction cost.
    const double premium = straddlePrice(S0, K, r, sigma, T) - tc;

    std::vector<double> pnl(nPaths, 0.0);

    for (int i = 0; i < nPaths; ++i) {
        const auto& path = paths[i];

        // ── t=0: set up the self-financing portfolio ───────────────────────
        // We hold `shares` of the underlying, funded by the cash account.
        // Short straddle delta = -(2N(d1) - 1).
        double shares = -straddleDelta(S0, K, r, sigma, T);

        // Cash account starts with premium received, minus cost of
        // purchasing the initial hedge position.
        double cash = premium - shares * S0;

        // ── Step through the path ─────────────────────────────────────────
        for (int j = 1; j <= nSteps; ++j) {
            const double S_curr = path[j];
            const double tau    = T - j * dtPath;

            // Accrue interest on cash balance over this step.
            cash *= std::exp(r * dtPath);

            // Rebalance at hedge intervals (but not at the final step —
            // that is handled by the expiry settlement block below).
            if (j % stepsPerHedge == 0 && tau > 1e-10) {
                const double newShares = -straddleDelta(S_curr, K, r, sigma, tau);
                const double trade     = newShares - shares;

                // Buying `trade` shares costs `trade * S_curr` from cash.
                // Negative trade = selling = cash inflow.
                cash   -= trade * S_curr;

                // Proportional transaction cost on the notional traded.
                cash   -= tc * std::abs(trade * S_curr);

                shares  = newShares;
            }
        }

        // ── Expiry settlement ─────────────────────────────────────────────
        const double ST = path[nSteps];

        // 1. Liquidate the residual share position at ST.
        cash += shares * ST;

        // 2. Pay the straddle holder (we are short).
        const double callPayoff = std::max(ST - K, 0.0);
        const double putPayoff  = std::max(K  - ST, 0.0);
        cash -= (callPayoff + putPayoff);

        // 3. `cash` IS the PnL — it already contains premium, hedge gains/losses,
        //    interest accrual, and the final settlement. Nothing else to add.
        pnl[i] = cash;
    }

    // ── Statistics ────────────────────────────────────────────────────────
    const double meanPnL = std::accumulate(pnl.begin(), pnl.end(), 0.0) / nPaths;

    double sq = 0.0;
    for (double p : pnl) sq += (p - meanPnL) * (p - meanPnL);
    const double stdPnL = std::sqrt(sq / std::max(1, nPaths - 1));
    const double sharpe = (stdPnL > 1e-12) ? meanPnL / stdPnL : 0.0;

    // Max drawdown: worst peak-to-trough on sorted running-mean equity,
    // not on the raw sequential pnl vector (which has no time ordering).
    // Report cross-sectional worst-case instead: max single-path loss.
    const double maxDD = -*std::min_element(pnl.begin(), pnl.end());

    // VaR and ES at 95%
    std::vector<double> sorted = pnl;
    std::sort(sorted.begin(), sorted.end());
    const int    varIdx = std::max(0, static_cast<int>(std::floor(0.05 * nPaths)) - 1);
    const double var95  = -sorted[varIdx];

    double esSum = 0.0;
    const int esCount = std::max(1, varIdx + 1);
    for (int i = 0; i < esCount; ++i)
        esSum += sorted[i];
    const double es95 = -esSum / esCount;

    // PnL path = sorted distribution (useful for plotting the histogram)
    return BacktestResult{meanPnL, stdPnL, sharpe, maxDD, var95, es95, sorted};
}