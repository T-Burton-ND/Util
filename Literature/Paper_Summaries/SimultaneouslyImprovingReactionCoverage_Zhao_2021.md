---
Paper Title: Simultaneously improving reaction coverage and computational cost in automated reaction prediction tasks
DOI: https://doi.org/10.1038/s43588-021-00101-3
First Author: Qiyuan Zhao
Last Author: Brett M Savoie
Institution: Purdue university
My Title:  Introducing YARP
---

## 🧩 Takeaways
1. 
2. 
3. 
---
## 🔬 Relevance Summary
(Explain how this paper connects to your work — e.g., implications for YARP’s CRN pruning, potential metrics to test, methods to adapt, etc.)


---
## Chat-GPT Summary
This Nature Computational Science article introduces **Yet Another Reaction Program (YARP)**—a cost-efficient, automated reaction discovery platform that improves both **reaction coverage** and **transition-state convergence**.  
By uniting graph-based enumeration (via the b2f2 ERS), joint reactant–product geometry optimization, and semi-empirical (GFN2-xTB) pre-screening, YARP achieves:

- ~**100× reduction** in DFT cost compared with prior GSM-based approaches,  
- **>85% TS success rates**, and  
- **Broad mechanistic diversity** including both concerted and multistep reactions.

Benchmarks against the Zimmerman dataset, Grambow’s 3-hydroperoxypropanal network, and Yang’s Diels–Alder systems show that YARP **recovers all known reactions**, identifies **many new kinetically favorable products**, and provides detailed energetic maps (∆G‡, ∆Hf) across entire networks.  

The discussion highlights that simple ERSs (b2f2) capture most physical reactivity while remaining well-conditioned for TS convergence—suggesting that deeper, kinetically prioritized searches (as later realized in **YAKS**) can scale linearly rather than exponentially.  
Overall, the paper establishes YARP as a **scalable, physically grounded framework** for automated reaction prediction and a precursor to modern kinetic-exploration tools like **YAKS**.



