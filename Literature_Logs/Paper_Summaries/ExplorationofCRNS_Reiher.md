---
Paper Title: The Exploration of Chemical Reaction Networks
Url: https://www.annualreviews.org/content/journals/10.1146/annurev-physchem-071119-040123
First Author: Jan Unsleber
Last Author: Markus Reiher
Institution: ETH Zurich
My Title:  A review of the state of CRN exploration, 2020
Tags: [YAKS, YARP, CRN,]
---

# 🧭 3×5 Literature Sprint (20–30 min)

## 1️⃣ Warm-Up (5–7 min)
**Goal:** Get the “shape” of the paper.

- [X] Read Abstract  
- [X] Read 1st + last paragraph of Introduction  
- [X] Skim all Figure captions  
- [X] Read Conclusion paragraph  
- [X] Ask yourself:
  - *What question is this paper answering?*
Simultaneously addressing 2 key questions: what is the current state of CRN exploration, and where does that 
leave us positioned for upcoming developments? In a way, this paper directly predicts the emergence of **YAKS** and
other CRN exploration algorithims that are more than purely synthetic chemically focused. 
  - *What approach or tool are they using?*
They try and establish a benchmark system--or at least propose benchmark systems--and key parameters and metrics
for measuring software packages and solutions against each other.
  - *What results do I expect?*
I can see the development of **YAKS** to meet some of the metrics, particularly on the **YARP** side, but it has
yet to even attempt to reach other metrics (user availability and automation, for example). 
📝 **Quick Impression (1-3 sentences):**
> The limitations of current methods outlined in this paper seem to persist, but are rapidly being addressed only 5 years later. The remaining limitation is probably more systemic to academia then it is intrinsic to the problem.

---

## 2️⃣ Deep Dive (10–15 min)
**Goal:** Extract core contributions and translate to your domain.

- [X] Read one key section (Methods *or* Results) carefully (I read section 4, **Algorithim Features**) 
- [X] Note 3 takeaways or findings  
- [X] Highlight one anchor figure, table, or equation  
- [X] Translate ideas into your research language (e.g., YARP, CRNs, GED metrics)

### 🧩 Takeaways
1. Breadth and depth are necessary to resolve side-chain reactions and multi-step reaction kinetics fully *breadth* corresponds to *graph fidelity*, *depth* to *node fidelity.*
2. Can we explore high barrier reactions in low level only to get more reaction coverage? Is it worth it? Could we do *breadth with ML barriers* and *depth with more accurate quantum chemical barriers*? quote-un-quote, **situational analysis?**
3. I **need** to get my hands on the latest **YAKS** code and starting prepping it. That needs to be my end goal, is situational *YAKS* that can utilize ML barriers and low-level *YARP* for extreme breadth and maintain the high level *YARP/DFT* analysis for depth and situational analysis.
4. In CRN exploration as in ML, uncertainty quantification is vital and publishable.
5. both **YARP** and **YAKS** have to be scalable, computationally tractable, and extremely user-friendly if they are ever, ever going to used by anyone outside the group. How can I help achieve that? how do I turn that goal into a research outcome?
6. Data needs to be resuable and transferable. How can I help make that happen? Can we convert RGD-1 into a YARP/YAKS friendly library/cache that can be accessed to accelerate explorations? In the long term, can *toy-reactions* be stored in a cache as a starting point for lager molecule exploration?

📊 **Anchor figure / concept:**  

> Figure 1 - 3 types of exploration algorithim. Forward open-end, backward open-start, and start-to-end. YAKS/YARP are currently type 1. We want to have the option for type 3. 

🧠 **Key terms or ideas to revisit:**  

> Certain fascets are required for a useful algorithim, and many are laid out here, including uer availability and automation, uncertainty quantification, and the need for transferability.

---

## 3️⃣ Reflection (5–8 min)
**Goal:** Lock in understanding through context and synthesis.

### 🔬 Relevance Summary
(How this connects to your work — e.g., implications for YARP’s CRN pruning, potential metrics to test, or methods to adapt.)
> The work towards double ended metrics will help us expand from type 1 into type 3 exploration metics. Simultaneously, a major goal that must be balanced with research is the publishing of both **YARP** and **YAKS** in user-friendly *enough* formats that other groups can begin to use them and improve upon them. Data sharing and transferability can be secondary, but the reliable algorithim that dynamically provides accurate results with quantifiable uncertainty is paramount.

**If I applied this idea, what problem would it solve?**  
> Applying all these ideas will take a significant effort and years of work, but to actually contribute originally to the field this might be my best path forward. Dont I pride myself on patience and on the ability to achieve goals that others give up on?

**How is it similar/different from another paper I’ve read?**  
> Kind of my first paper of this kind, where we are tackling broad needs in the field rather than applicaitons. Would be great to look around and see if there is a more recent follow up. 

---

## 🗂️ Review Checklist
- [X] Saved in `/Paper_Notes/`  
- [X] Added to `reading_log.csv`  
- [X] Tagged for later synthesis (#metrics, #CRN, #automation, etc.)
- [X] Scheduled for re-read in 2 weeks  

---

✅ *End of Sprint — close the tab, you’re done.*
