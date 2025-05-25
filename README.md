# Wait-Less Online Decision-Making

## 1. **Abstract** ##

Online linear programming (OLP) has found broad applications in revenue management and resource allocation. State-of-the-art OLP algorithms achieve low regret by repeatedly solving linear programming (LP) subproblems that incorporate updated resource information. However, LP-based methods are computationally expensive and often inefficient for large-scale applications. By contrast, recent first-order OLP algorithms are more computationally efficient but typically suffer from weaker regret guarantees. To address these shortcomings, we propose a new algorithm that combines the strengths of LP-based and first-order OLP algorithms. Our algorithm re-solves the LP subproblems periodically at a predefined frequency $f$ and uses the latest dual prices to guide online decision-making. In parallel, a first-order method runs during each interval between LP re-solves and smooths resource consumption. Our algorithm achieves $\mathcal{O}(\log (T/f) + \sqrt{f})$ regret and delivers a "wait-less" online decision-making process that balances computational efficiency and regret guarantees. Extensive experiments demonstrate at least 10-fold improvements in regret over first-order methods and 100-fold improvements in runtime over LP-based methods. 

## 2. **Experiments** ##

We conduct extensive experiments to evaluate our algorithm's performance and validate our theoretical results. In the first part, we evaluate our main algorithms across different choices of re-solving frequency. In the second part, we compare our algorithms with LP-based and first-order methods in terms of regret and running time. All implementations are in MATLAB. We organize those files as follows: 

- main: include the main structure of experiments and different algorithms
- olptwopath_freq: include Algorithm 1 in our paper, frequently solving LP with subgradient fine-tune only in the first and last batches
- olptwopath_freq2: include Algorithm 2 in our paper, frequently solving LP with subgradient fine-tune throughout the whole horizon
- olp_infrequent: include an algorithm from other papers for comparison, infrequently solving LP
- olpgurobi: include the solver used for the LP-based method
- olptwopath_grad: include the first-order method
- olpgetdata: include data generation for the customer's bidding price (reward), resource consumption, and total resource
